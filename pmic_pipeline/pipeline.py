#!/usr/bin/env python3

"""
PMIC Pipeline: Antimicrobial resistance analysis pipeline.

This pipeline analyzes NGS data for antimicrobial resistance genes and virulence factors.
Supports Illumina paired-end reads and Oxford Nanopore reads.
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def plot_gene(data_gene_file, png_file_plot, title):
    """
    Plot gene coverage and identity from abricate output.

    Args:
        data_gene_file (str): Path to the CSV file with gene data.
        png_file_plot (str): Path to save the plot.
        title (str): Title for the plot.
    """
    df = pd.read_csv(data_gene_file, sep='\t', index_col='GENE')
    # Filter for coverage >= 60 and identity >= 95
    df = df.loc[(df["%COVERAGE"] >= 60) & (df["%IDENTITY"] >= 95)]
    df_filter = df[['%COVERAGE', '%IDENTITY']]
    if df_filter.shape[0] > 1:
        df_filter.plot(kind='bar', xlabel='GENE', ylabel='IDENTITY & COVERAGE', rot=20, figsize=(20, 4), title=title)
        plt.savefig(png_file_plot)
        plt.close()  # Close to free memory
    else:
        logger.info("No genes meet the filter criteria for coverage and identity.")


def plot_gene_illumina(data_gene_file, png_file_plot, title):
    """
    Plot gene coverage and identity for Illumina data.

    Args:
        data_gene_file (str): Path to the TSV file with gene data.
        png_file_plot (str): Path to save the plot.
        title (str): Title for the plot.
    """
    df = pd.read_csv(data_gene_file, sep='\t', index_col='GENE')
    df = df.loc[(df["%COVERAGE"] >= 60) & (df["%IDENTITY"] >= 95)]
    df_filter = df[['%COVERAGE', '%IDENTITY']]
    if df_filter.shape[0] > 1:
        df_filter.plot(kind='bar', xlabel='GENE', ylabel='IDENTITY & COVERAGE', rot=20, figsize=(20, 4), title=title)
        plt.savefig(png_file_plot)
        plt.close()
    else:
        logger.info("No genes meet the filter criteria for coverage and identity.")


def download_file(url, dest):
    """
    Download a file from a URL using wget.

    Args:
        url (str): Download URL.
        dest (str): Local file path to save.
    """
    run_command(['wget', '-c', url, '-O', dest])


def download_mash_refseq(dest_dir):
    """
    Fetch a prebuilt Mash RefSeq database sketch from NCBI FTP.

    The script downloads the sketch file and places it under `dest_dir`.
    Users can then point `mash dist` to this file when running the pipeline.
    """
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    # URL for a public Mash sketch; adjust if NCBI changes location
    url = 'https://ftp.ncbi.nlm.nih.gov/refseq/release/mash/refseq.genomes.k21s1000.msh'
    out_file = dest / 'refseq.genomes.k21s1000.msh'
    logger.info(f"Downloading Mash RefSeq sketch to {out_file}")
    download_file(url, str(out_file))
    logger.info("Mash RefSeq database download complete.")


def download_plasmid_db(dest_dir):
    """
    Download a plasmid sequence database from NCBI.

    Currently fetches the RefSeq plasmid FASTA file. The resulting file can be
    formatted into a BLAST/abricate database as needed by the pipeline.
    """
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    url = 'https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.1.genomic.fna.gz'
    out_file = dest / 'plasmid.fna.gz'
    logger.info(f"Downloading plasmid database to {out_file}")
    download_file(url, str(out_file))
    logger.info("Plasmid database download complete.")


def run_command(cmd, cwd=None, check=True, stdout=None):
    """
    Run a shell command using subprocess.

    Args:
        cmd (list): Command as list of strings.
        cwd (str): Working directory.
        check (bool): Raise exception on non-zero exit.
        stdout: File object to redirect stdout to, or None to capture.

    Returns:
        subprocess.CompletedProcess: Result of the command.
    """
    logger.info(f"Running command: {' '.join(cmd)}")
    try:
        if stdout is not None:
            result = subprocess.run(cmd, cwd=cwd, check=check, stdout=stdout, stderr=subprocess.PIPE, text=True)
        else:
            result = subprocess.run(cmd, cwd=cwd, check=check, capture_output=True, text=True)
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        logger.error(f"Error: {e.stderr}")
        raise


def download_reference(accession, output_file):
    """
    Download a reference sequence from NCBI.

    Args:
        accession (str): Accession number.
        output_file (str): Path to save the FASTA file.
    """
    logger.info(f"Downloading reference for accession: {accession}")
    Entrez.email = 'ndaom403@gmail.com'  # Should be configurable
    with Entrez.efetch(db='nuccore', id=accession, rettype='fasta') as handle:
        with open(output_file, 'w') as f:
            f.write(handle.read())


def parse_mash_accession(line):
    """
    Extract the accession part from a mash output line containing
    "genomic.fna.gz:". The text after the colon may include an optional
    bracketed sequence count, e.g.

        genomic.fna.gz:[1870 seqs] NC_004354.4
        genomic.fna.gz:NC_012967.1

    This function returns the accession (e.g. "NC_004354.4" or
    "NC_012967.1"). Raises ValueError if the pattern cannot be matched.
    """
    m = re.search(r'genomic\.fna\.gz:(?:\[[^\]]+\]\s*)?(?P<acc>\S+)', line)
    if not m:
        raise ValueError(f"Unable to extract accession from mash line: {line!r}")
    return m.group('acc')


def run_mapping(reference, reads1, reads2, output_prefix, threads=12, rg_id=None, rg_sm=None):
    """
    Run BWA mapping pipeline.

    Args:
        reference (str): Path to reference FASTA.
        reads1 (str): Path to forward reads.
        reads2 (str): Path to reverse reads.
        output_prefix (str): Prefix for output files.
        threads (int): Number of threads.
        rg_id (str): optional read group ID to inject into alignments.
        rg_sm (str): optional sample name for the read group.
    """
    # Index reference if necessary
    run_command(['bwa', 'index', reference])

    # Map reads, attach a read group so Picard will not crash later
    sam_file = f"{output_prefix}.sam"
    mem_cmd = ['bwa', 'mem', '-t', str(threads)]
    if rg_id or rg_sm:
        parts = []
        if rg_id:
            parts.append(f'ID:{rg_id}')
        if rg_sm:
            parts.append(f'SM:{rg_sm}')
        rg_string = '@RG\\t' + '\\t'.join(parts)
        mem_cmd += ['-R', rg_string]
    mem_cmd += [reference, reads1, reads2]
    run_command(mem_cmd, stdout=open(sam_file, 'w'))

    # Convert to BAM
    bam_file = f"{output_prefix}.bam"
    run_command(['samtools', 'view', '-O', 'BAM', sam_file], stdout=open(bam_file, 'w'))

    # Sort BAM
    sorted_bam = f"{output_prefix}.sorted.bam"
    run_command(['samtools', 'sort', '-T', 'temp', '-O', 'bam', '-o', sorted_bam, bam_file])

    # Index sorted BAM
    run_command(['samtools', 'index', sorted_bam])

    # Mark duplicates using new Picard syntax and require read groups
    markdup_bam = f"{output_prefix}_markdup.bam"
    metrics_file = f"{output_prefix}_metrics.txt"
    run_command(['picard', 'MarkDuplicates', '-I', sorted_bam, '-O', markdup_bam, '-M', metrics_file])

    # Index markdup BAM
    run_command(['samtools', 'index', markdup_bam])

    # Stats
    stats_file = f"{output_prefix}_bamstats.txt"
    run_command(['samtools', 'stats', markdup_bam], stdout=open(stats_file, 'w'))

    # Clean up
    os.remove(sam_file)
    os.remove(bam_file)

    return markdup_bam


def run_variant_calling(reference, bam_file, output_prefix):
    """
    Run variant calling with bcftools.

    Args:
        reference (str): Path to reference FASTA.
        bam_file (str): Path to BAM file.
        output_prefix (str): Prefix for output files.
    """
    # Mpileup
    bcf_file = f"{output_prefix}_variant.bcf"
    #run_command(['samtools', 'mpileup', '-t', 'DP', '-Bug', '-m', '4', '-f', reference, bam_file], stdout=open(bcf_file, 'w'))
    run_command(['bcftools', 'mpileup', '--annotate', 'FORMAT/DP', '-B', '-m', '4', '-f', reference, bam_file, '-O', 'b', '-o', bcf_file])

    # Call variants
    vcf_file = f"{output_prefix}_variant.vcf"
    run_command(['bcftools', 'call', '-mv', '-O', 'v', '-o', vcf_file, bcf_file])

    # Filter variants
    filtered_vcf = f"{output_prefix}_variant_filtered.vcf"
    run_command(['bcftools', 'filter', '-i', 'type="snp" && QUAL>=50 && FORMAT/DP>5 && MQ>=30', '-g10', '-G10', vcf_file, '-o', filtered_vcf])

    # Compress and index
    filtered_vcf_gz = f"{filtered_vcf}.gz"
    run_command(['bcftools', 'view', '-O', 'z', '-o', filtered_vcf_gz, filtered_vcf])
    run_command(['bcftools', 'index', filtered_vcf_gz])

    # Consensus
    consensus_fasta = f"{output_prefix}_consensus.fasta"
    run_command(['bcftools', 'consensus', '-f', reference, filtered_vcf_gz], stdout=open(consensus_fasta, 'w'))

    return consensus_fasta


def run_abricate(fasta_file, db, output_file):
    """
    Run abricate for gene detection.

    Args:
        fasta_file (str): Path to FASTA file.
        db (str): Database name.
        output_file (str): Path to output CSV.
    """
    run_command(['abricate', '-db', db, fasta_file], stdout=open(output_file, 'w'))

def get_plasmid_db_path():
    """Return the filesystem path to the plasmid BLAST database.

    Uses the BLASTDB environment variable if set, otherwise defaults to
    "$HOME/blastdb/plsdb". This ensures any user can run the pipeline
    without hardcoding a specific home directory.
    """
    base = os.environ.get('BLASTDB', os.path.expanduser('~/blastdb'))
    return str(Path(base) / 'plsdb')


def process_illumina_sample(sample_name, forward, reverse, output_dir, threads=12):
    """
    Process a single Illumina sample.

    Args:
        sample_name (str): Sample name.
        forward (str): Path to forward reads.
        reverse (str): Path to reverse reads.
        output_dir (str): Output directory.
        threads (int): Number of threads.
    """
    output_path = Path(output_dir) / sample_name
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Processing Illumina sample: {sample_name}")

    # Assembly with SPAdes
    run_command(['spades.py', '-t', str(threads), '--only-assembler', '--plasmid', '-1', forward, '-2', reverse, '-o', str(output_path)])

    # Extract contigs and identify plasmids
    scaffolds_fasta = output_path / "scaffolds.fasta"
    if not scaffolds_fasta.exists():
        logger.error(f"Scaffolds file not found: {scaffolds_fasta}")
        return

    # Get component numbers (plasmid identifiers)
    result = run_command(['grep', 'component', str(scaffolds_fasta)])
    components = [line.split('_')[7] for line in result.stdout.strip().split('\n') if line]
    components = list(set(components))

    logger.info(f"Found {len(components)} plasmids")

    # Get sequence IDs
    result = run_command(['grep', '>', str(scaffolds_fasta)])
    id_list = [line[1:] for line in result.stdout.strip().split('\n') if line]

    num_plasmid = [seqid[-2:] for seqid in id_list]
    num_plasmid = list(set(num_plasmid))

    for num in num_plasmid:
        sequence_ids = [seqid for seqid in id_list if seqid.endswith(num)]

        # Extract plasmid contigs
        contig_fasta = output_path / f"{sample_name}_contig_plasmid{num}.fasta"
        with open(contig_fasta, 'w') as f:
            for seqid in sequence_ids:
                run_command(['perl', '-ne', f'if(/^>(\\S+)/){{$c=grep{{/^{seqid}$/}}qw({" ".join(sequence_ids)})}}print if $c', str(scaffolds_fasta)], stdout=f)

        # BLAST against plasmid database
        plasmid_blast = output_path / f"{sample_name}_plasmid_like{num}.csv"
        run_command(['blastn', '-query', str(contig_fasta), '-db', get_plasmid_db_path(), '-outfmt', '10 qseqid sseqid qstart qend pident bitscore sstart send qlen slen', '-out', str(plasmid_blast)])

        if plasmid_blast.stat().st_size > 0:
            # Get reference accession (2nd field in BLAST CSV: sseqid)
            result = run_command(['awk', '-F,', '{print $2}', str(plasmid_blast)])
            accessions = result.stdout.strip().split('\n')
            # Skip header if present, take the first result's accession
            ref_accession = accessions[1] if len(accessions) > 1 else accessions[0]

            # Download reference
            ref_fasta = output_path / f"{sample_name}_plasmid_reference{num}.fasta"
            download_reference(ref_accession, str(ref_fasta))

            # Mapping (include read-group tags so Picard MarkDuplicates can run)
            bam_prefix = output_path / f"{sample_name}_plasmid{num}"
            markdup_bam = run_mapping(str(ref_fasta), forward, reverse, str(bam_prefix), threads,
                                      rg_id=sample_name + f"_plasmid{num}", rg_sm=sample_name)

            # Variant calling
            consensus = run_variant_calling(str(ref_fasta), markdup_bam, str(bam_prefix))

            # AMR and VF detection
            amr_csv = output_path / f"{sample_name}_plasmid_AMR{num}.csv"
            vf_csv = output_path / f"{sample_name}_plasmid_VF{num}.csv"
            run_abricate(str(consensus), 'card', str(amr_csv))
            run_abricate(str(consensus), 'vfdb', str(vf_csv))

            # Plot if genes found
            if amr_csv.stat().st_size > 122:
                plot_png = output_path / f"{sample_name}_plasmid_AMR{num}.png"
                plot_gene_illumina(str(amr_csv), str(plot_png), 'AMR of Plasmid')
            if vf_csv.stat().st_size > 122:
                plot_png = output_path / f"{sample_name}_plasmid_VF{num}.png"
                plot_gene_illumina(str(vf_csv), str(plot_png), 'VF of Plasmid')
        else:
            logger.info(f"No plasmid found for {contig_fasta}")

    # Chromosome analysis
    # Mash for species identification
    combined_reads = output_path / "dna_samples.fastq.gz"
    run_command(['cat', forward, reverse], stdout=open(str(combined_reads), 'w'))
    run_command(['mash', 'sketch', '-m', '2', str(combined_reads)])
    distance_csv = output_path / f"{sample_name}_chr_distance.csv"
    run_command(['mash', 'dist', '-C', 'refseq.genomes.k21s1000.msh', f"{combined_reads}.msh"], stdout=open(str(distance_csv), 'w'))
    os.remove(str(combined_reads))
    os.remove(f"{combined_reads}.msh")

    # Select best reference
    df_chr = pd.read_csv(str(distance_csv), sep="\t", names=['Query_ID', 'Query_samples', 'Query_dist', 'Median_multiplicity', 'Shared_hashes'])
    tri_distance_chr = df_chr.sort_values(by=['Query_dist']).head()
    tri_distance_csv = output_path / f"{sample_name}_chr_distance_trier.csv"
    tri_distance_chr.to_csv(str(tri_distance_csv))

    # Get accession
    result = run_command(['awk', '-F","', '{print $1}', str(tri_distance_csv)])
    accessions = result.stdout.strip().split('\n')
    # find the first line containing the genome indicator
    line = next((l for l in accessions[1:] if 'genomic.fna.gz:' in l), None)
    if line is None:
        raise ValueError("No line containing 'genomic.fna.gz:' found in the first column")
    ref_accession = parse_mash_accession(line)
    logger.debug(f"Extracted accession {ref_accession} from mash line: {line}")

    # Download reference
    chr_ref_fasta = output_path / f"{sample_name}_chr_reference.fasta"
    download_reference(ref_accession, str(chr_ref_fasta))

    # Mapping for chromosome
    chr_bam_prefix = output_path / f"{sample_name}_chr"
    chr_markdup_bam = run_mapping(str(chr_ref_fasta), forward, reverse, str(chr_bam_prefix), threads,
                                  rg_id=sample_name + "_chr", rg_sm=sample_name)

    # Variant calling
    chr_consensus = run_variant_calling(str(chr_ref_fasta), chr_markdup_bam, str(chr_bam_prefix))

    # AMR and VF
    chr_amr_csv = output_path / f"{sample_name}_chr_AMR.csv"
    chr_vf_csv = output_path / f"{sample_name}_chr_VF.csv"
    run_abricate(str(chr_consensus), 'card', str(chr_amr_csv))
    run_abricate(str(chr_consensus), 'vfdb', str(chr_vf_csv))

    # Plot
    if chr_amr_csv.stat().st_size > 122:
        plot_png = output_path / f"{sample_name}_chr_AMR.png"
        plot_gene_illumina(str(chr_amr_csv), str(plot_png), 'AMR of Chromosome')
    if chr_vf_csv.stat().st_size > 122:
        plot_png = output_path / f"{sample_name}_chr_VF.png"
        plot_gene_illumina(str(chr_vf_csv), str(plot_png), 'VF of Chromosome')

    logger.info(f"Finished processing sample: {sample_name}")


def process_nanopore_sample(sample_name, reads, output_dir, threads=12):
    """
    Process a single Nanopore sample.

    Args:
        sample_name (str): Sample name.
        reads (str): Path to reads file.
        output_dir (str): Output directory.
        threads (int): Number of threads.
    """
    output_path = Path(output_dir) / sample_name
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Processing Nanopore sample: {sample_name}")

    # Assembly with Flye
    run_command(['flye', '--nano-raw', reads, '--threads', str(threads), '-o', str(output_path)])

    assembly_fasta = output_path / "assembly.fasta"
    if not assembly_fasta.exists():
        logger.error(f"Assembly file not found: {assembly_fasta}")
        return

    # Mash for species identification
    combined_reads = output_path / "dna_samples.fastq.gz"
    run_command(['cat', reads], stdout=open(str(combined_reads), 'w'))
    run_command(['mash', 'sketch', str(combined_reads)])
    distance_csv = output_path / f"{sample_name}_chr_distance.csv"
    run_command(['mash', 'dist', '-C', 'refseq.genomes.k21s1000.msh', f"{combined_reads}.msh"], stdout=open(str(distance_csv), 'w'))
    os.remove(str(combined_reads))
    os.remove(f"{combined_reads}.msh")

    # Select best reference
    df_chr = pd.read_csv(str(distance_csv), sep="\t", names=['Query_ID', 'Query_samples', 'Query_dist', 'Median_multiplicity', 'Shared_hashes'])
    tri_distance_chr = df_chr.sort_values(by=['Query_dist']).head()
    tri_distance_csv = output_path / f"{sample_name}_chr_distance_trier.csv"
    tri_distance_chr.to_csv(str(tri_distance_csv))

    result = run_command(['awk', '-F\t', '{print $1}', str(tri_distance_csv)])
    accessions = result.stdout.strip().split('\n')
    line = next((l for l in accessions[1:] if 'genomic.fna.gz:' in l), None)
    if line is None:
        raise ValueError("No line containing 'genomic.fna.gz:' found in the first column")
    ref_accession = parse_mash_accession(line)

    # Download reference
    chr_ref_fasta = output_path / f"{sample_name}_chr_reference.fasta"
    download_reference(ref_accession, str(chr_ref_fasta))

    # AMR and VF detection on assembly
    amr_csv = output_path / f"{sample_name}_AMR.csv"
    vf_csv = output_path / f"{sample_name}_VF.csv"
    plasmid_csv = output_path / f"{sample_name}_plasmid_like.csv"
    run_abricate(str(assembly_fasta), 'card', str(amr_csv))
    run_abricate(str(assembly_fasta), 'vfdb', str(vf_csv))
    run_abricate(str(assembly_fasta), 'plasmidfinder', str(plasmid_csv))

    # Plot
    if amr_csv.stat().st_size > 122:
        plot_png = output_path / f"{sample_name}_AMR.png"
        plot_gene(str(amr_csv), str(plot_png), 'AMR of Assembly')
    if vf_csv.stat().st_size > 122:
        plot_png = output_path / f"{sample_name}_VF.png"
        plot_gene(str(vf_csv), str(plot_png), 'VF of Assembly')

    # Plasmid analysis
    result = run_command(['awk', '-F,', '{print $2}', str(plasmid_csv)])
    contigs_plasmid = [line for line in result.stdout.strip().split('\n') if line and line != 'SEQUENCE']
    contigs_plasmid = list(set(contigs_plasmid))

    for contig in contigs_plasmid:
        logger.info(f"Processing plasmid contig: {contig}")

        # Extract contig
        contig_fasta = output_path / f"{sample_name}_{contig}.fasta"
        with open(contig_fasta, 'w') as f:
            run_command(['perl', '-ne', f'if(/^>(\\S+)/){{$c=grep{{/^{contig}$/}}qw({contig})}}print if $c', str(assembly_fasta)], stdout=f)

        # BLAST
        blast_csv = output_path / f"{sample_name}_{contig}_plasmid_ref.csv"
        run_command(['blastn', '-query', str(contig_fasta), '-db', '~/blastdb/plsdb', '-outfmt', '10 qseqid sseqid qstart qend pident bitscore start send qlen slen', '-out', str(blast_csv)])

        # AMR and VF on contig
        contig_amr = output_path / f"{sample_name}_{contig}_AMR_plasmid.csv"
        contig_vf = output_path / f"{sample_name}_{contig}_VF_plasmid.csv"
        run_abricate(str(contig_fasta), 'card', str(contig_amr))
        run_abricate(str(contig_fasta), 'vfdb', str(contig_vf))

        # Plot
        if contig_amr.stat().st_size > 122:
            plot_png = output_path / f"{sample_name}_{contig}_AMR_plasmid.png"
            plot_gene(str(contig_amr), str(plot_png), 'AMR of Plasmid Contig')
        if contig_vf.stat().st_size > 122:
            plot_png = output_path / f"{sample_name}_{contig}_VF_plasmid.png"
            plot_gene(str(contig_vf), str(plot_png), 'VF of Plasmid Contig')

        # Get reference for contig
        df = pd.read_csv(str(blast_csv), sep=",", names=["qseqid", "sseqid", "qstart", "qend", "pident", "bitscore", "start", "send", "qlen", "slen"])
        tri_distance = df.sort_values(by=["pident", "bitscore"], ascending=False).head()
        tri_csv = output_path / f"{sample_name}_{contig}_blast_trie.csv"
        tri_distance.to_csv(str(tri_csv))

        result = run_command(['awk', '-F,', '{print $2}', str(tri_csv)])
        accessions = result.stdout.strip().split('\n')
        # Take the top hit's accession (skip header)
        ref_accession = accessions[1] if len(accessions) > 1 else accessions[0]

        # Download reference
        ref_fasta = output_path / f"{sample_name}_{contig}_reference.fasta"
        download_reference(ref_accession, str(ref_fasta))

    logger.info(f"Finished processing sample: {sample_name}")

def setup_databases():
    """
    Download required databases for PMIC pipeline.
    """
    from pathlib import Path
    import os

    blastdb_dir = Path(os.path.expanduser('/home/genomic/blastdb'))
    blastdb_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading Mash RefSeq database")

    mash_url = "https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh"
    
    run_command([
        "wget",
        "-O",
        "refseq.genomes.k21s1000.msh",
        mash_url
    ])

    logger.info("Downloading NCBI plasmid database")

    plasmid_url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.1.genomic.fna.gz"

    plasmid_fna_gz = blastdb_dir / "plsdb.fasta.gz"

    run_command([
        "wget",
        "-O",
        str(plasmid_fna_gz),
        plasmid_url
    ])

    plasmid_fna = blastdb_dir / "plsdb.fasta"

    run_command([
        "gunzip",
        str(plasmid_fna_gz)
    ])

    logger.info("Creating BLAST database")

    run_command([
        "makeblastdb",
        "-in",
        str(plasmid_fna),
        "-dbtype",
        "nucl",
        "-out",
        str(blastdb_dir / "plsdb")
    ])

    logger.info("Databases installed successfully")

def main():
    parser = argparse.ArgumentParser(
        description="PMIC pipeline: Antimicrobial resistance analysis from NGS data"
    )

    parser.add_argument(
        "-in", "--input",
        required=False,
        default=None,
        help="Input directory containing samples"
    )

    parser.add_argument(
        "-out", "--output",
        required=False,
        default=None,
        help="Output directory"
    )

    parser.add_argument(
        "--illumina",
        action="store_true",
        help="Illumina paired-end reads"
    )

    parser.add_argument(
        "--nanopore",
        action="store_true",
        help="Oxford Nanopore reads"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=12,
        help="Number of threads to use"
    )

    parser.add_argument(
        "--download-databases",
        action="store_true",
        help="Download required databases (Mash RefSeq and NCBI plasmid)"
    )

    args = parser.parse_args()

    if args.download_databases:
        setup_databases()
        sys.exit(0)

    # -in and -out are required for normal pipeline runs
    if not args.input or not args.output:
        parser.error("Arguments -in/--input and -out/--output are required unless --download-databases is used")

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.illumina:
        logger.info("Running Illumina pipeline")
        import re
        # Match files ending in _R1.fastq.gz or _R2.fastq.gz
        pattern = re.compile(r'^(.+)_R([12])\.fastq\.gz$')
        samples = [f for f in input_dir.iterdir() if f.is_file()]
        sample_dict = {}
        for sample in samples:
            m = pattern.match(sample.name)
            if not m:
                logger.warning(f"Skipping unrecognized file: {sample.name}")
                continue
            name, read_num = m.group(1), m.group(2)
            if name not in sample_dict:
                sample_dict[name] = {}
            if read_num == '1':
                sample_dict[name]['forward'] = str(sample)
            else:
                sample_dict[name]['reverse'] = str(sample)

        for sample_name, reads in sample_dict.items():
            if 'forward' in reads and 'reverse' in reads:
                process_illumina_sample(sample_name, reads['forward'], reads['reverse'], str(output_dir), args.threads)

    elif args.nanopore:
        logger.info("Running Nanopore pipeline")
        samples = [f for f in input_dir.iterdir() if f.is_file()]
        for sample_file in samples:
            # use the filename without its extension(s) as sample name
            # e.g. "sample_name.fastq.gz" -> "sample_name"
            sample_name = sample_file.stem
            process_nanopore_sample(sample_name, str(sample_file), str(output_dir), args.threads)
    else:
        parser.error("Please specify --illumina or --nanopore")

    logger.info("Pipeline completed successfully")


if __name__ == "__main__":
    main()