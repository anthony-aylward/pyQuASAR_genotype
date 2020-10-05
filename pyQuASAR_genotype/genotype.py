#!/usr/bin/env python3
#===============================================================================
# genotype.py
#===============================================================================

# Imports ======================================================================

import json
import os
import os.path
import pandas as pd
import pyhg19
import pyQuASAR
import seqalign
import subprocess
import tempfile
import wasp_map

from argparse import ArgumentParser
from functools import partial
from itertools import groupby
from multiprocessing import Pool




# Functions ====================================================================

def prepare_quasar_input(
    input_file_path: str,
    bam_dir: str,
    intermediate_dir: str,
    reference_genome_path: str,
    mapping_quality: int,
    blacklist_path: str,
    snps_path: str,
    processes: int,
    memory: int,
    paired_end: bool,
    skip_preprocessing: bool,
    write_bam: bool,
    algorithm_switch_bp: int = 70,
    library_name=None,
    algorithm=None,
    temp_dir=None
) -> str:
    """Format data into input files for QuASAR
    
    Parameters
    ----------
    input_file_path : str
        Path to an input file
    bam_dir : str
        Directory to write BAM files
    intermediate_dir : str
        Directory to write intermediate pileup / bed files
    reference_genome_path : str
        Path to reference genome
    mapping_quality : int
        Minimum quality score for filtering alignment
    blacklist_path : str
        Path to ENCODE mappability blacklist
    snps_path : str
        Path to file containing SNPs to genotype
    processes : int
        Number of processes
    memory : int
        Memory limit
    paired_end : bool
        Indicator for paired_end reads
    skip_preprocessing : bool
        Indicator to skip preprocessing steps
    write_bam : bool
        Indicator to write a BAM file to disk
    algorithm_switch_bp : int
        Read length threshold for switching to `bwa mem`
    algorithm : str or None
        Force use of either `aln` or `mem` algorithm, if supplied
    temp_dir
        directory to use for temporary files
    
    Returns
    -------
    str
        Path to a QuASAR input file
    """

    input_file_path = (
        input_file_path.split(',') if ',' in input_file_path
        else input_file_path
    )
    bam_prefix, intermediate_prefix = (
        os.path.join(
            directory,
            (
                library_name if library_name else os.path.basename(
                    input_file_path if isinstance(input_file_path, str)
                    else input_file_path[0]
                )
                .replace('.bam', '')
                .replace('.fastq', '')
                .replace('.gz', '')
            )
        )
        for directory in (bam_dir, intermediate_dir)
    )
    with open('{}.align.log'.format(bam_prefix), 'w') as log:
        sa = seqalign.SequenceAlignment(
            input_file=input_file_path,
            mapping_quality=mapping_quality,
            processes=processes,
            log=log,
            aligner=seqalign.BWA(
                reference_genome_path=reference_genome_path,
                trim_qual=15,
                algorithm_switch_bp=algorithm_switch_bp,
                algorithm=algorithm
            ),
            dedupper=wasp_map.RmDup(
                processes=processes,
                paired_end=paired_end,
                temp_dir=temp_dir
            ),
            temp_dir=temp_dir
        )
        if not skip_preprocessing:
            sa.apply_quality_filter()
            sa.remove_supplementary_alignments()
            sa.remove_blacklisted_reads(blacklist_path=blacklist_path)
            sa.samtools_sort(memory_limit=memory)
            sa.samtools_index()
            sa.restrict_chromosomes(*(tuple(range(1, 23)) + ('X',)))
            sa.samtools_index()
            if write_bam:
                sa.write('{}.filt.bam'.format(bam_prefix))
            sa.remove_duplicates()
            sa.samtools_sort()
        compressed_pileup_bed_path = f'{intermediate_prefix}.pileup.bed.gz'
        pyQuASAR.write_compressed_pileup_bed(
            sa.samtools_mpileup(
                positions=snps_path,
                reference_genome=reference_genome_path
            ),
            compressed_pileup_bed_path,
            snps_bed_path=snps_path
        )
    pyQuASAR.bed_to_quasar(compressed_pileup_bed_path)
    quasar_input_file_path = '{}.quasar.in.gz'.format(intermediate_prefix)
    return (
        quasar_input_file_path if os.path.isfile(quasar_input_file_path)
        else None
    )


def prepare_quasar_input_from_metadata(
    input_file_path: str,
    library_name: str,
    bam_dir: str,
    intermediate_dir: str,
    reference_genome_path: str,
    mapping_quality: int,
    blacklist_path: str,
    snps_path: str,
    processes: int,
    memory: int,
    paired_end: bool,
    skip_preprocessing: bool,
    write_bam: bool,
    algorithm_switch_bp: int = 70,
    algorithm=None,
    temp_dir=None
) -> str:
    return prepare_quasar_input(
        input_file_path=input_file_path,
        bam_dir=bam_dir,
        intermediate_dir=intermediate_dir,
        reference_genome_path=reference_genome_path,
        mapping_quality=mapping_quality,
        blacklist_path=blacklist_path,
        snps_path=snps_path,
        processes=processes,
        memory=memory,
        paired_end=paired_end,
        skip_preprocessing=skip_preprocessing,
        write_bam=write_bam,
        algorithm_switch_bp=algorithm_switch_bp,
        library_name=library_name,
        algorithm=library_name,
        temp_dir=temp_dir
    )


def generate_collated_metadata(metadata_dict):
    for e, experiment in metadata_dict.items():
        for l, library in experiment['libraries'].items():
            if isinstance(library, str):
                yield (
                    os.path.join(experiment['dir'], f"{experiment['assay']}.{library}.fastq.gz"),
                    '.'.join((experiment['assay'], l))
                )
            elif isinstance(library, list) and len(library) == 2:
                yield (
                    ','.join(os.path.join(experiment['dir'], f"{experiment['assay']}.{lib}.fastq.gz") for lib in library),
                    '.'.join((experiment['assay'], l))
                )


def count_input_paths(metadata_item):
    return len(metadata_item[0].split(','))


def collate_metadata(metadata_dict):
    meta = sorted(
        (k, tuple(g)) for k, g in groupby(
            generate_collated_metadata(metadata_dict), key=count_input_paths
        )
    )
    meta_se = tuple(g for k, g in meta if k == 1)
    meta_pe = tuple(g for k, g in meta if k == 2)
    return meta_se, meta_pe



def get_genotypes(
    single_end: list,
    paired_end: list,
    metadata: str,
    bam_dir: str,
    intermediate_dir: str,
    reference_genome_path: str,
    mapping_quality: int,
    blacklist_path: str,
    snps_path: str,
    processes: int,
    memory: int,
    skip_preprocessing: bool = False,
    write_bam: bool = False,
    algorithm_switch_bp: int = 70,
    algorithm=None,
    temp_dir=None
):
    """Obtain genotypes from sequencing data using QuASAR
    
    Parameters
    ----------
    single_end : list
        List of single-end input files
    paired_end : list
        List of paired-end input files
    metadata : dict
        Dict of input file metadata
    bam_dir : str
        Directory to write BAM files
    intermediate_dir : str
        Directory to write intermediate pileup / bed files
    reference_genome_path : str
        Path to reference genome
    mapping_quality : int
        Minimum quality score for filtering alignment
    blacklist_path : str
        Path to ENCODE mappability blacklist
    snps_path : str
        Path to file containing SNPs to genotype
    processes : int
        Number of processes
    memory : int
        Memory limit
    skip_preprocessing : bool
        Indicator to skip preprocessing steps
    write_bam : bool
        Indicator to write a BAM file to disk
    algorithm_switch_bp : int
        Read length threshold for switching to `bwa mem`
    algorithm : str or None
        Force use of either `aln` or `mem` algorithm, if supplied
    temp_dir
        directory to use for temporary files
    """
    
    n_single_end = len(single_end)
    n_paired_end = len(paired_end)
    if not metadata:
        metadata_dict = {}
    else:
        with open(metadata, 'r') as f:
            metadata_dict = json.load(f)
    n_metadata = sum(len(x['libraries']) for x in metadata_dict.values())

    def prepare_quasar_input_params(temp_dir_name, n, pe=False):
        return {
            'bam_dir': bam_dir if bam_dir else temp_dir_name,
            'intermediate_dir': (
                intermediate_dir if intermediate_dir
                else temp_dir_name
            ),
            'reference_genome_path': reference_genome_path,
            'mapping_quality': mapping_quality,
            'blacklist_path': blacklist_path,
            'snps_path': snps_path,
            'processes': max(1, int(processes / n)),
            'memory': memory / min(processes, n),
            'paired_end': pe,
            'skip_preprocessing': skip_preprocessing,
            'write_bam': write_bam,
            'algorithm_switch_bp': algorithm_switch_bp,
            'algorithm': algorithm,
            'temp_dir': temp_dir
        }
    
    with tempfile.TemporaryDirectory(dir=temp_dir) as temp_dir_name:
        with Pool(processes=min(processes, max(n_single_end, n_paired_end, n_metadata))) as pool:
            if n_single_end > 0:
                single_end_quasar_input_paths = pool.map(
                    partial(
                        prepare_quasar_input,
                        **prepare_quasar_input_params(temp_dir_name, n_single_end, pe=False)
                    ),
                    single_end
                )
            else:
                single_end_quasar_input_paths = []
            
            if n_paired_end > 0:
                paired_end_quasar_input_paths = pool.map(
                    partial(
                        prepare_quasar_input,
                        **prepare_quasar_input_params(temp_dir_name, n_paired_end, pe=True)
                    ),
                    paired_end
                )
            else:
                paired_end_quasar_input_paths = []
            
            if n_metadata > 0:
                meta_se, meta_pe = collate_metadata(metadata_dict)
                if len(meta_se) > 0:
                    metadata_quasar_input_paths_se = pool.starmap(
                        partial(
                            prepare_quasar_input_from_metadata,
                            **prepare_quasar_input_params(temp_dir_name, len(meta_se), pe=False)
                        ),
                        meta_se[0]
                    )
                else:
                    metadata_quasar_input_paths_se = []
                if len(meta_pe) > 0:
                    metadata_quasar_input_paths_pe = pool.starmap(
                        partial(
                            prepare_quasar_input_from_metadata,
                            **prepare_quasar_input_params(temp_dir_name, len(meta_pe), pe=True)
                        ),
                        meta_pe[0]
                    )
                else:
                    metadata_quasar_input_paths_pe = []
            else: 
                metadata_quasar_input_paths_se, metadata_quasar_input_paths_pe = [], []
        return pyQuASAR.genotype(
            *filter(
                None,
                single_end_quasar_input_paths
                + paired_end_quasar_input_paths
                + metadata_quasar_input_paths_se
                + metadata_quasar_input_paths_pe
            )
        )


def parse_arguments():
    parser = ArgumentParser(
        description=(
            'Infer genotypes from ChIP-seq or ATAC-seq data using QuASAR'
        )
    )
    io_group = parser.add_argument_group('I/O arguments')
    io_group.add_argument(
        'single_end',
        metavar='<path/to/single_end_data.{fa/fq/bam}>',
        nargs='*',
        help='Paths to single-end FASTQ or BAM files'
    )
    io_group.add_argument(
        '--paired-end',
        metavar='<path/to/paired_end_data.{fa/fq/bam}>',
        nargs='+',
        default=[],
        help='Paths to paired-end FASTQ or BAM files. Paired FASTQ files should be joined by a comma.'
    )
    io_group.add_argument(
        '--metadata',
        metavar='<path/to/metadata.json>',
        help='path to JSON file containing input metadata'
    )
    io_group.add_argument(
        '--bam-dir',
        metavar='<path/to/bam_dir/>',
        help='directory in which to place BAM files'
    )
    io_group.add_argument(
        '--inter-dir',
        metavar='<path/to/inter_dir/>',
        help='prefix for intermediate files'
    )
    io_group.add_argument(
        '--vcf-chr',
        metavar='<output/vcf/prefix>',
        help='Prefix for output VCFs split by chromosome'
    )
    io_group.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress printing of VCF file to standard output'
    )
  
    align_group = parser.add_argument_group('Alignment arguments')
    align_group.add_argument(
        '--processes',
        metavar='<int>',
        type=int,
        default=1,
        help='Number of processes to use [1]'
    )
    align_group.add_argument(
        '--memory',
        metavar='<int>',
        type=int,
        default=8,
        help='Maximum memory usage in GB [8]'
    )
    align_group.add_argument(
        '--quality',
        metavar='<int>',
        type=int,
        default=10,
        help='Mapping quality cutoff for samtools [10]'
    )
    align_group.add_argument(
        '--reference',
        metavar='<path/to/reference_genome.fa>',
        default=pyhg19.PATH,
        help=f'Path to reference genome prepared for BWA [{pyhg19.PATH}]'
    )
    align_group.add_argument(
        '--blacklist',
        metavar='<path/to/blacklist.bed>',
        default=pyhg19.BLACKLIST,
        help=f'Path to ENCODE blacklist file [{pyhg19.BLACKLIST}]'
    )
    align_group.add_argument(
        '--write-bam',
        action='store_true',
        help='Write bam files to disk'
    )
    align_group.add_argument(
        '--algorithm-switch-bp',
        metavar='<int>',
        default=70,
        help='Read length threshold for switching to `bwa mem` [70]'
    )
    align_group.add_argument(
        '--algorithm',
        choices={'aln', 'mem'},
        default=None,
        help='Force use of either `bwa aln` or bwa mem`'
    )
    
    quasar_group = parser.add_argument_group('QuASAR arguments')
    quasar_group.add_argument(
        '--model-fitting-snps',
        metavar='<path/to/snps_file.bed>',
        default=pyQuASAR.SNPS_BED_PATH,
        help=f'BED file containing 1KGP SNPs [{pyQuASAR.SNPS_BED_PATH}]'
    )
    quasar_group.add_argument(
        '--query-snps',
        metavar='<path/to/snps_file.bed>',
        help='BED file containing query SNPs'
    )
    quasar_group.add_argument(
        '--skip-preprocessing',
        action='store_true',
        help='skip preprocessing steps'
    )
    
    vcf_group = parser.add_argument_group('VCF arguments')
    vcf_group.add_argument(
        '--sample',
        metavar='<sample_id>',
        default='SAMPLE',
        help='Name for the sample [SAMPLE]'
    )
    vcf_group.add_argument(
        '--threshold',
        metavar='<float>',
        default=0.99,
        type=float,
        help='Probability threshold for genotype calls [0.99]'
    )
    vcf_group.add_argument(
        '--het-only',
        action='store_true',
        help='Output heterozygous variants only'
    )

    config_group = parser.add_argument_group('configuration arguments')
    config_group.add_argument(
        '--tmp-dir',
        metavar='<temp/file/dir/>',
        help='directory to use for temporary files'
    )
    args = parser.parse_args()
    if len(args.single_end) + len(args.paired_end) == 0:
        if not args.metadata:
            raise RuntimeError('No input files provided')
    return args


def main():
    args = parse_arguments()
    if args.query_snps:
        df_fitting = pd.read_table(args.model_fitting_snps, header=None)
        fitting_coord_set = {(chrom, pos) for chrom, pos in zip(df_fitting.iloc[:,0], df.fitting.iloc[:,2])}
        df_query = pd.read_table(args.query_snps, header=None)
        df_query_filt = pd.DataFrame(r for i, r in df_query.iterrows() if (r[0], r[2]) not in fitting_coord_set)
        df_combined = pd.concat((df_fitting, df_query_filt)).sort_values([0, 2])
        with tempfile.NamedTemporaryFile(dir=args.tmp_dir) as temp:
            temp_file_name = temp.name
        df_combined.to_csv(temp_file_name, sep='\t', index=False, header=False)
        snps_path = temp_file_name
    else:
        snps_path = args.model_fitting_snps
    vcf = pyQuASAR.genotype_to_vcf(
        get_genotypes(
            args.single_end,
            args.paired_end,
            args.metadata,
            args.bam_dir,
            args.inter_dir,
            args.reference,
            args.quality,
            args.blacklist,
            snps_path,
            args.processes,
            args.memory,
            skip_preprocessing=args.skip_preprocessing,
            write_bam=args.write_bam,
            algorithm_switch_bp=args.algorithm_switch_bp,
            algorithm=args.algorithm,
            temp_dir=args.tmp_dir
        ),
        sample_name=args.sample,
        snps_bed_path=snps_path,
        threshold=args.threshold,
        het_only=args.het_only,
        temp_file_dir=args.tmp_dir
    )
    if args.vcf_chr and not args.quiet:
        vcf = tuple(vcf)
    if args.vcf_chr:
        pyQuASAR.write_split_vcf(vcf, args.vcf_chr)
    os.remove(temp_file_name)
    if not args.quiet:
        for line in vcf:
            print(line)




# Execute ======================================================================

if __name__ == '__main__':
    main()
