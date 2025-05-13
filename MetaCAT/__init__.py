import os
import platform
from argparse import ArgumentParser, RawTextHelpFormatter
from uuid import uuid4

from . import abundance
from . import abundance_test
from . import benchmark_gt
from . import benchmark_rw
from . import checkm2
from . import cluster
from . import cluster2mapping
from . import coverage
from . import gtdbtk
from . import index_bam
from . import mapping2cluster
from . import mwas
from . import plot_mwas
from . import representative
from . import seed
from . import variant
from .c import findBinary
from .colors import isValidColor
from .io_test import findExecutablePath, isReadable, isWriteable


def __init__():
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'MetaCAT: Metagenome Clustering and Association Tool',
        epilog = 'congcong_liu@icloud.com.'
    )
    parser.add_argument(
        '-v', '--version', action = 'version', version = '%(prog)s 1.0.0'
    )

    subparsers = parser.add_subparsers(
        title = 'command', dest = 'command', required = True,
        description = 'Try "%(prog)s [command] -h|--help" for full help.',
    )

    # coverage parser #
    coverage_parser = subparsers.add_parser(
        'coverage', formatter_class = RawTextHelpFormatter,
        help = 'Generate COVERAGE file from bam files.',
        usage = '%(prog)s [options] -b <BAMs> -o <COVERAGE>'
    )
    coverage_parser.add_argument(
        '-b', '--bam', type = str, nargs = '+', required = True, metavar = '<str>',
        help = 'Path to the sorted bam files.\nAll bam files must have the same header.'
    )
    coverage_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the coverage file.'
    )
    coverage_parser.add_argument(
        '-ti', '--threads-index', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads for indexing files.\nWorks only with multiple bam files.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    coverage_parser.add_argument(
        '-tc', '--threads-count', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads for counting depth.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    coverage_parser.add_argument(
        '-m', '--min-sequence-length', default = 100, type = int, required = False, metavar = '<int>',
        help = 'Calculations applies only to the sequences with length being greater than or equal to this value.\nThe value should be greater than 2 * "trim" + 1.\nDefault: 100.'
    )
    coverage_parser.add_argument(
        '-l', '--long', default = False, action = 'store_true', required = False,
        help = 'Compute coverage for long-read bam files.\nDefault: False.'
    )
    coverage_parser.add_argument(
        '--trim', default = 0, type = int, required = False, metavar = '<int>',
        help = 'Ignore the regions at both ends of a sequence when do calculations.\nThe value should be a non-negative integer.\nDefault: 0.'
    )
    coverage_parser.add_argument(
        '--mapq', default = 0, type = int, required = False, metavar = '<int>',
        help = 'MAPQ threshold.\nThe value should be a non-negative integer.\nDefault: 0.'
    )
    coverage_parser.add_argument(
        '--identity', default = 0.90, type = float, required = False, metavar = '<float>',
        help = 'Identity threshold.\nThe value should be from 0 to 1.\nDefault: 0.90.'
    )

    # seed parser #
    seed_parser = subparsers.add_parser(
        'seed', formatter_class = RawTextHelpFormatter,
        help = 'Generate SEED file from an assembly file.',
        usage = '%(prog)s [options] -f <ASSEMBLY> -o <SEED>'
    )
    seed_parser.add_argument(
        '-f', '--fasta', type = str, required = True, metavar = '<str>',
        help = 'Path to the fasta formatted assembly file.'
    )
    seed_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the seed file.'
    )
    seed_parser.add_argument(
        '-t', '--threads', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    seed_parser.add_argument(
        '--fraggenescan', type = str, required = False, metavar = '<str>',
        help = 'Path to the "FragGeneScan".\nDefault: internal FragGeneScan.'
    )
    seed_parser.add_argument(
        '--hmmpress', type = str, required = False, metavar = '<str>',
        help = 'Path to the "hmmpress".\nDefault: internal hmmpress.'
    )
    seed_parser.add_argument(
        '--hmmsearch', type = str, required = False, metavar = '<str>',
        help = 'Path to the "hmmsearch".\nDefault: internal hmmsearch.'
    )

    # cluster parser #
    cluster_parser = subparsers.add_parser(
        'cluster', formatter_class = RawTextHelpFormatter,
        help = 'Cluster sequences based on ASSEMBLY, COVERAGE and SEED files.',
        usage = '%(prog)s [options] -f <ASSEMBLY> -c <COVERAGE> -s <SEED> -o <METACAT>'
    )
    cluster_parser.add_argument(
        '-f', '--fasta', type = str, required = True, metavar = '<str>',
        help = 'Path to the fasta formatted assembly file.'
    )
    cluster_parser.add_argument(
        '-c', '--coverage', type = str, required = True, metavar = '<str>',
        help = 'Path to the coverage file.'
    )
    cluster_parser.add_argument(
        '-s', '--seed', type = str, required = True, metavar = '<str>',
        help = 'Path to the seed file.'
    )
    cluster_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Prefix of the output files.\n\"*.*.fasta\" the fasta formatted cluster file.\n\"*.mapping\" the mapping of sequences to clusters.'
    )
    cluster_parser.add_argument(
        '-m', '--min-sequence-length', type = int, required = False, metavar = '<int>',
        help = 'Sequences with the length being greater than or equal to this value can be involved in the clustering algorithm.\nThe value should be a positive integer.\nDefault: 1500 (#samples = 1), 500 (#samples > 1).'
    )
    cluster_parser.add_argument(
        '-t', '--threads', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    cluster_parser.add_argument(
        '--no-fasta-clusters', default = False, action = 'store_true', required = False,
        help = 'Do not output fasta formatted cluster files.\nDefault: False.'
    )
    cluster_parser.add_argument(
        '--min-pca-variance-ratio', default = 0.90, type = float, required = False, metavar = '<float>',
        help = 'Minimum percentage of variance explained in PCA in SWDPGMM.\nThe value should be from 0 to 1.\nDefault: 0.90.'
    )
    cluster_parser.add_argument(
        '--min-pca-components', default = 20, type = int, required = False, metavar = '<int>',
        help = 'Minimum number of components retained in PCA in SWDPGMM.\nThe value should be a positive integer.\nDefault: 20.'
    )
    cluster_parser.add_argument(
        '--swdpgmm-engine', default = 'auto', type = str, choices = ('auto', 'cpu', 'gpu'), required = False,
        help = 'Device used to run SWDPGMM.\nDefault: auto (SWDPGMM will automatically run on the GPU if it is available).'
    )
    cluster_parser.add_argument(
        '--min-swdpgmm-clusters', default = 100, type = int, required = False, metavar = '<int>',
        help = 'SWDPGMM is enabled if the estimated number of genomes is greater than or equal to this value.\nThe value should be an integer.\nDefault: 100.'
    )
    cluster_parser.add_argument(
        '--min-swdpgmm-sequences', default = 500000, type = int, required = False, metavar = '<int>',
        help = 'SWDPGMM is enabled if the number of sequences that meet the length threshold is greater than or equal to this value.\nThe value should be an integer.\nDefault: 500000.'
    )
    cluster_parser.add_argument(
        '--max-swdpgmm-iterations', default = 30, type = int, required = False, metavar = '<int>',
        help = 'Maximum number of SWDPGMM iterations to perform.\nThe value should be a positive integer.\nDefault: 30.'
    )
    cluster_parser.add_argument(
        '--kmer-frequence-weight', default = [0.10, 0.30, 0.50, 0.70, 0.90], type = float, nargs = '+', required = False, metavar = '<float>',
        help = 'Weights of kmer frequency probabilistic model for constructing the affinity graph.\nThe values should be from 0 to 1.\nDefault: 0.10 0.30 0.50 0.70 0.90.'
    )
    cluster_parser.add_argument(
        '--min-kmer-frequence-probability', default = 0.50, type = float, required = False, metavar = '<float>',
        help = 'Probability of kmer frequence probabilistic model lower than this value will be set to 0.\nThe value should be from 0 to 1.\nDefault: 0.50.'
    )
    cluster_parser.add_argument(
        '--coverage-weight', default = [0.90, 0.70, 0.50, 0.30, 0.10], type = float, nargs = '+', required = False, metavar = '<float>',
        help = 'Weight of read coverage probabilistic model for constructing the affinity graph.\nThe values should be from 0 to 1.\nDefault: 0.90 0.70 0.50 0.30 0.10.'
    )
    cluster_parser.add_argument(
        '--min-coverage-probability', default = 0.50, type = float, required = False, metavar = '<float>',
        help = 'Probability of coverage probabilistic model lower than this value will be set to 0.\nThe value should be from 0 to 1.\nDefault: 0.50.'
    )
    cluster_parser.add_argument(
        '--neighbors', default = 10, type = int, required = False, metavar = '<int>',
        help = 'Number of neighbors in affinity model.\nThe values should be positive integer.\nDefault: 10.'
    )
    cluster_parser.add_argument(
        '--min-cluster', default = 100000, type = int, required = False, metavar = '<int>',
        help = 'Minimum size of a cluster to output.\nThe value should be a positive integer.\nDefault: 100000 bp.'
    )
    cluster_parser.add_argument(
        '--max-seeds', default = 30000, type = int, required = False, metavar = '<int>',
        help = 'Maximum number of seed sequences involved in the seed model.\nThe value should be a positive integer.\nDefault: 30000.'
    )
    cluster_parser.add_argument(
        '--random-number', default = 0, type = int, required = False, metavar = '<int>',
        help = 'Random number generator seeded by the given integer.\nDefault: 0.'
    )

    # cluster2mapping parser #
    cluster2mapping_parser = subparsers.add_parser(
        'cluster2mapping', formatter_class = RawTextHelpFormatter,
        help = 'Generate mapping file from fasta formatted cluster files.',
        usage = '%(prog)s -f <FASTAs> -o <MAPPING>'
    )
    cluster2mapping_parser.add_argument(
        '-f', '--fasta', type = str, nargs = '+', required = True, metavar = '<str>',
        help = 'Path to the fasta formatted cluster files.'
    )
    cluster2mapping_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )

    # mapping2cluster parser #
    mapping2cluster_parser = subparsers.add_parser(
        'mapping2cluster', formatter_class = RawTextHelpFormatter,
        help = 'Generate fasta formatted cluster files from assembly and mapping file.',
        usage = '%(prog)s [options] -f <ASSEMBLY> -m <MAPPING>'
    )
    mapping2cluster_parser.add_argument(
        '-f', '--fasta', type = str, required = True, metavar = '<str>',
        help = 'Path to the fasta formatted assembly file.'
    )
    mapping2cluster_parser.add_argument(
        '-m', '--mapping', type = str, required = True, metavar = '<str>',
        help = 'Path to the mapping file.\nA header line "Sequence ID<tab>Cluster ID" should be present.'
    )
    mapping2cluster_parser.add_argument(
        '-o', '--output', default = os.getcwd(), type = str, required = False, metavar = '<str>',
        help = f'Directory to output fasta formatted cluster files.\nDefault: {os.getcwd()}.'
    )
    mapping2cluster_parser.add_argument(
        '--checkm2', type = str, required = False, metavar = '<str>',
        help = 'Path to the file generated by MetaCAT\'s checkm2, only available if "--cluster" is disable.\nA header line "Name<tab>Completeness<tab>Contamination ..." should be present.'
    )
    mapping2cluster_parser.add_argument(
        '--contamination', default = 0.10, type = float, required = False, metavar = '<float>',
        help = 'Contamination threshold, only available if "--checkm2" is enabled.\nThe value should be from 0 to 1.\nDefault: 0.10.'
    )
    mapping2cluster_parser.add_argument(
        '--completeness', default = 0.70, type = float, required = False, metavar = '<float>',
        help = 'Completeness threshold, only available if "--checkm2" is enabled.\nThe value should be from 0 to 1.\nDefault: 0.70.'
    )
    mapping2cluster_parser.add_argument(
        '--min-cluster', default = 1, type = int, required = False, metavar = '<int>',
        help = 'Minimum size of a cluster to output.\nThe value should be a positive integer.\nDefault: 1.'
    )
    mapping2cluster_parser.add_argument(
        '--cluster', type = str, nargs = '+', required = False, metavar = '<str>',
        help = 'Output only the selected clusters.\nDefault: all clusters.'
    )

    # benchmarkGT parser #
    benchmarkGT_parser = subparsers.add_parser(
        'benchmarkGT', formatter_class = RawTextHelpFormatter,
        help = 'Benchmark with ground truth.',
        usage = '%(prog)s [options] -gt <GROUND-TRUTH> -m <MAPPINGs> -o <BENCHMARK>'
    )
    benchmarkGT_parser.add_argument(
        '-gt', '--ground-truth', type = str, required = True, metavar = '<str>',
        help = 'Ground truth file is tab-delimited, contains 3 columns and a single header line as follows.\nSequence ID<tab>Genome ID<tab>Length.'
    )
    benchmarkGT_parser.add_argument(
        '-m', '--mapping', type = str, nargs = '+', required = True, metavar = '<str>',
        help = 'Path to the mapping files.\nA header line "Sequence ID<tab>Cluster ID" should be present.'
    )
    benchmarkGT_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Prefix of the output files.'
    )
    benchmarkGT_parser.add_argument(
        '--dataset', type = str, required = False, metavar = '<str>',
        help = 'Title of dataset.\nDefault: value of "--ground-truth".'
    )
    benchmarkGT_parser.add_argument(
        '--label', type = str, nargs = '+', required = False, metavar = '<str>',
        help = 'Labels of the mapping files.\nIf "--label" is selected, the number of values should be equal to the number of mapping files.'
    )
    benchmarkGT_parser.add_argument(
        '--precision', default = [0.95, 0.90], type = float, nargs = '+', required = False, metavar = '<float>',
        help = 'Precision thresholds to plot.\nThe value should be from 0 to 1.\nDefault: 0.95 0.90.'
    )
    benchmarkGT_parser.add_argument(
        '--recall', default = [0.95, 0.90, 0.80, 0.70], type = float, nargs = '+', required = False, metavar = '<float>',
        help = 'Recall thresholds to plot.\nThe value should be from 0 to 1.\nDefault: 0.95 0.90 0.80 0.70.'
    )
    benchmarkGT_parser.add_argument(
        '--color', default = 'Reds_r', type = str, required = False, metavar = '<str>',
        help = 'The color string.\nDefault: Reds_r.'
    )
    benchmarkGT_parser.add_argument(
        '--color-min', default = 0.20, type = float, required = False, metavar = '<float>',
        help = 'Minimum color.\nThe value should be from 0 to 1.\nDefault: 0.20.'
    )
    benchmarkGT_parser.add_argument(
        '--color-max', default = 0.75, type = float, required = False, metavar = '<float>',
        help = 'Maximum color.\nThe value should be from 0 to 1.\nDefault: 0.75.'
    )

    # benchmarkRW parser #
    benchmarkRW_parser = subparsers.add_parser(
        'benchmarkRW', formatter_class = RawTextHelpFormatter,
        help = 'Benchmark for real-world datasets.',
        usage = '%(prog)s [options] --combine -i <CHECKM2> -o <BENCHMARK>'
    )
    benchmarkRW_parser.add_argument(
        '-i', '--input', type = str, required = True, metavar = '<str>',
        help = 'Path to the file generated by MetaCAT\'s checkm2.\nA header line "Name<tab>Completeness<tab>Contamination ..." should be present.\nThe format of the 1st column must be "Dataset.Program.ClusterID".\n"Dataset" and "ClusterID" fields must contain no dots.'
    )
    benchmarkRW_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Prefix of the output files.'
    )
    benchmarkRW_parser.add_argument(
        '-c', '--combine', default = False, action = 'store_true', required = False,
        help = 'Combine all datasets.\nDefault: False.'
    )
    benchmarkRW_parser.add_argument(
        '--combined-dataset', default = 'Total Datasets', type = str, required = False, metavar = '<str>',
        help = 'Title of combined dataset.\nDefault: Total Datasets.'
    )
    benchmarkRW_parser.add_argument(
        '--contamination', default = [0.05, 0.10], type = float, nargs = '+', required = False, metavar = '<float>',
        help = 'Contamination thresholds to plot.\nThe value should be from 0 to 1.\nDefault: 0.05 0.10.'
    )
    benchmarkRW_parser.add_argument(
        '--completeness', default = [0.95, 0.90, 0.80, 0.70], type = float, nargs = '+', required = False, metavar = '<float>',
        help = 'Completeness thresholds to plot.\nThe value should be from 0 to 1.\nDefault: 0.95 0.90 0.80 0.70.'
    )
    benchmarkRW_parser.add_argument(
        '--rows', type = int, required = False, metavar = '<int>',
        help = 'Number of rows of the subplot grid.\nThe value should be a positive integer.\nDefault: auto.'
    )
    benchmarkRW_parser.add_argument(
        '--columns', type = int, required = False, metavar = '<int>',
        help = 'Number of columns of the subplot grid.\nThe value should be a positive integer.\nDefault: auto.'
    )
    benchmarkRW_parser.add_argument(
        '--program', type = str, nargs = '+', required = False, metavar = '<str>',
        help = 'Order of programs in the figure.'
    )
    benchmarkRW_parser.add_argument(
        '--label', type = str, nargs = '+', required = False, metavar = '<str>',
        help = 'Labels of the programs.\nIf "--label" is selected, the number of values should be equal to the number of programs.'
    )
    benchmarkRW_parser.add_argument(
        '--color', default = 'Reds_r', type = str, required = False, metavar = '<str>',
        help = 'The color string.\nDefault: Reds_r.'
    )
    benchmarkRW_parser.add_argument(
        '--color-min', default = 0.20, type = float, required = False, metavar = '<float>',
        help = 'Minimum color.\nThe value should be from 0 to 1.\nDefault: 0.20.'
    )
    benchmarkRW_parser.add_argument(
        '--color-max', default = 0.75, type = float, required = False, metavar = '<float>',
        help = 'Maximum color.\nThe value should be from 0 to 1.\nDefault: 0.75.'
    )

    # indexBam parser #
    indexBam_parser = subparsers.add_parser(
        'indexBam', formatter_class = RawTextHelpFormatter,
        help = 'Index bam files.',
        usage = '%(prog)s [options] -i <BAMs>'
    )
    indexBam_parser.add_argument(
        '-i', '--input', nargs = '+', type = str, required = True, metavar = '<str>',
        help = 'Path to the sorted bam files.\nAll bam files must have the same header.'
    )
    indexBam_parser.add_argument(
        '-t', '--threads', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads.\nWorks only with multiple bam files.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )

    # representative parser #
    representative_parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'Select representative.',
        epilog = 'congcong_liu@icloud.com.'
    )
    representative_parser = subparsers.add_parser(
        'representative', formatter_class = RawTextHelpFormatter,
        help = 'Select representatives from fasta formatted cluster files.',
        usage = '%(prog)s [options] -f <FASTAs> -c <CHECKM2> -g <GTDBTK> -o <REPRESENTATIVE>'
    )
    representative_parser.add_argument(
        '-f', '--fasta', nargs = '+', type = str, required = True, metavar = '<str>',
        help = 'Path to the fasta files.'
    )
    representative_parser.add_argument(
        '-c', '--checkm2', type = str, required = True, metavar = '<str>',
        help = 'Path to the file generated by MetaCAT\'s checkm2.\nA header line "Name<tab>Completeness<tab>Contamination<tab>...<tab>...<tab>...<tab>Contig_N50 ..." should be present.'
    )
    representative_parser.add_argument(
        '-g', '--gtdbtk', type = str, required = True, metavar = '<str>',
        help = 'Path to the file generated by MetaCAT\'s gtdbtk.\nA header line "user_genome<tab>classification ..." should be present.'
    )
    representative_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Prefix of the output files.\n\"*.annotation\": the classifications of all clusters.\n\"*.assembly\": the combined representative assembly.\n\"*.mapping\": the mapping of sequences to clusters.'
    )
    representative_parser.add_argument(
        '-t', '--threads', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    representative_parser.add_argument(
        '--mash', type = str, required = False, metavar = '<str>',
        help = 'Path to the "mash".\nDefault: internal mash.'
    )
    representative_parser.add_argument(
        '--distance', default = 0.05, type = float, required = False, metavar = '<float>',
        help = 'Maximum mash distance.\nThe value should be from 0 to 1.\nDefault: 0.05.'
    )
    representative_parser.add_argument(
        '--kmer-size', default = 21, type = int, required = False, metavar = '<int>',
        help = 'K-mer size used in Mash.\nThe value should be a positive integer.\nDefault: 21.'
    )
    representative_parser.add_argument(
        '--sketch-size', default = 1000, type = int, required = False, metavar = '<int>',
        help = 'Number of non-redundant min-hashes used in Mash.\nThe value should be a positive integer.\nDefault: 1000.'
    )
    representative_parser.add_argument(
        '--contamination', default = 0.10, type = float, required = False, metavar = '<float>',
        help = 'Contamination threshold.\nThe value should be from 0 to 1.\nDefault: 0.10.'
    )
    representative_parser.add_argument(
        '--completeness', default = 0.70, type = float, required = False, metavar = '<float>',
        help = 'Completeness threshold.\nThe value should be from 0 to 1.\nDefault: 0.70.'
    )
    representative_parser.add_argument(
        '--rename', default = False, action = 'store_true', required = False,
        help = 'Sequences are stored in the assembly file from 1 to N.\nDefault: False.'
    )

    # abundance parser #
    abundance_parser = subparsers.add_parser(
        'abundance', formatter_class = RawTextHelpFormatter,
        help = 'Generate abundance table from metagenomic data.',
        usage = '%(prog)s [options] -a <ANNOTATION> -b <BAMs> -m <MAPPINGs> -o <ABUNDANCE>'
    )
    abundance_parser.add_argument(
        '-a', '--annotation', type = str, required = True, metavar = '<str>',
        help = 'Path to the annotation file generated by MetaCAT\'s representative.\nA header line "Cluster ID<tab>Classification" should be present.'
    )
    abundance_parser.add_argument(
        '-b', '--bam', type = str, nargs = '+', required = True, metavar = '<str>',
        help = 'Path to the sorted bam files.'
    )
    abundance_parser.add_argument(
        '-m', '--mapping', type = str, nargs = '+', required = True, metavar = '<str>',
        help = 'Path to the mapping files.\nA header line "Sequence ID<tab>Cluster ID" should be present.\nThe number of values should be equal to the number of bam files, if the number of values is not one.\nOtherwise all bam files must have the same header.'
    )
    abundance_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    abundance_parser.add_argument(
        '-ti', '--threads-index', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads for indexing files.\nWorks only with multiple bam files.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    abundance_parser.add_argument(
        '-tc', '--threads-count', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads for counting depth.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    abundance_parser.add_argument(
        '-l', '--long', default = False, action = 'store_true', required = False,
        help = 'Compute abundance for long-read bam files.\nDefault: False.'
    )
    abundance_parser.add_argument(
        '--trim', default = 0, type = int, required = False, metavar = '<int>',
        help = 'Ignore the regions at both ends of a sequence when do calculations.\nThe value should be a non-negative integer.\nDefault: 0.'
    )
    abundance_parser.add_argument(
        '--mapq', default = 0, type = int, required = False, metavar = '<int>',
        help = 'MAPQ threshold.\nThe value should be a non-negative integer.\nDefault: 0.'
    )
    abundance_parser.add_argument(
        '--identity', default = 0.90, type = float, required = False, metavar = '<float>',
        help = 'Identity threshold.\nThe value should be from 0 to 1.\nDefault: 0.90.'
    )
    abundance_parser.add_argument(
        '--min-abundance', default = 0.10, type = float, required = False, metavar = '<float>',
        help = 'Absolute abundance less than this value will be set to 0.\nThe value should be from 0 to 1.\nDefault: 0.10.'
    )

    # abundanceTest parser #
    abundanceTest_parser = subparsers.add_parser(
        'abundanceTest', formatter_class = RawTextHelpFormatter,
        help = 'Significance test for abundance.',
        usage = '%(prog)s [options] -a <ABUNDANCE> -g <GROUP> -o <TEST>'
    )
    abundanceTest_parser.add_argument(
        '-a', '--abundance', type = str, required = True, metavar = '<str>',
        help = 'Path to the abundance file.\nA header line "Abundance<tab>ID1<tab>ID2 ..." should be present.'
    )
    abundanceTest_parser.add_argument(
        '-g', '--group', type = str, required = True, metavar = '<str>',
        help = 'Path to the group file.\nA header line "ID<tab>Group" should be present.'
    )
    abundanceTest_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Prefix of the output files.\n\"*.G1-G2.*: significance test for group G1 and G2.\"'
    )
    abundanceTest_parser.add_argument(
        '--coverage', default = 0.50, type = float, required = False, metavar = '<float>',
        help = 'Minimum coverage in each of groups.\nThe value should be from 0 to 1.\nDefault: 0.50.'
    )
    abundanceTest_parser.add_argument(
        '--rows', type = int, required = False, metavar = '<int>',
        help = 'Number of rows of the subplot grid.\nThe value should be a positive integer.\nDefault: auto.'
    )
    abundanceTest_parser.add_argument(
        '--columns', type = int, required = False, metavar = '<int>',
        help = 'Number of columns of the subplot grid.\nThe value should be a positive integer.\nDefault: auto.'
    )
    abundanceTest_parser.add_argument(
        '--classification', nargs = '+', default = ['s'], type = str, choices = ('s', 'g', 'f', 'o', 'c', 'p', 'd'), required = False, metavar = '<str>',
        help = 'Classification levels used for comparison.\nThe values can be s, g, f, o, c, p, d.\nDefault: s.'
    )
    abundanceTest_parser.add_argument(
        '--comparison', nargs = 2, type = str, required = False, metavar = '<str>',
        help = 'Paired groups for comparison.\nDefault: all paired groups.'
    )
    abundanceTest_parser.add_argument(
        '--alpha', default = 0.05, type = float, required = False, metavar = '<float>',
        help = 'The probability of making the wrong decision when the null hypothesis is true.\nDefault: 0.05.'
    )
    abundanceTest_parser.add_argument(
        '--multiple-test', default = 'bonferroni', type = str, choices = ('bonferroni', 'benjamini-hochberg', 'none'), required = False, metavar = '<str>',
        help = 'Method for multiple hypothesis testing.\nThe values can be bonferroni, benjamini-hochberg, none.\nDefault: bonferroni.'
    )

    # variant parser #
    variant_parser = subparsers.add_parser(
        'variant', formatter_class = RawTextHelpFormatter,
        help = 'Call SNPs from metagenomic data.',
        usage = '%(prog)s [options] -i <BAMs> -o <VARIANT>'
    )
    variant_parser.add_argument(
        '-i', '--input', nargs = '+', type = str, required = True, metavar = '<str>',
        help = 'Path to the sorted bam files.\nAll bam files must have the same header.'
    )
    variant_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    variant_parser.add_argument(
        '-tc', '--threads-call', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads for calling variants.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    variant_parser.add_argument(
        '-ti', '--threads-index', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads for indexing files.\nWorks only with multiple bam files.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    variant_parser.add_argument(
        '--coverage', default = 0.50, type = float, required = False, metavar = '<float>',
        help = 'Minimum population coverage of a SNP.\nDefault: 0.50.'
    )
    variant_parser.add_argument(
        '--mapq', default = 0, type = int, required = False, metavar = '<int>',
        help = 'MAPQ threshold.\nThe value should be a non-negative integer.\nDefault: 0.'
    )
    variant_parser.add_argument(
        '--identity', default = 0.90, type = float, required = False, metavar = '<float>',
        help = 'Identity threshold.\nThe value should be from 0 to 1.\nDefault: 0.90.'
    )
    variant_parser.add_argument(
        '--depth', default = 1, type = int, required = False, metavar = '<int>',
        help = 'Read depth for a sample less than this value will be set to 0.\nThe value should be a positive integer.\nDefault: 1.'
    )

    # mwas parser #
    mwas_parser = subparsers.add_parser(
        'mwas', formatter_class = RawTextHelpFormatter,
        help = 'MWAS for metagenomic data.',
        usage = '%(prog)s [options] -p <PHENOTYPE> -m <MAPPING> -ad <ABUNDANCE> -at <ANNOTATION> -v <VARIANT> -o <MWAS>'
    )
    mwas_parser.add_argument(
        '-p', '--phenotype', type = str, required = True, metavar = '<str>',
        help = 'Path to the phenotype file.\nA header line "ID<tab>Phenotype1<tab>Phenotype2 ..." should be present.\nEach phenotype is a numeric value.'
    )
    mwas_parser.add_argument(
        '-m', '--mapping', type = str, required = True, metavar = '<str>',
        help = 'Path to the mapping file.\nA header line "Sequence ID<tab>Cluster ID" should be present.'
    )
    mwas_parser.add_argument(
        '-ad', '--abundance', type = str, required = True, metavar = '<str>',
        help = 'Path to the abundance file.\nA header line "Abundance<tab>ID1<tab>ID2 ..." should be present.'
    )
    mwas_parser.add_argument(
        '-at', '--annotation', type = str, required = True, metavar = '<str>',
        help = 'Path to the annotation file generated by MetaCAT\'s representative.\nA header line "Cluster ID<tab>Classification" should be present.'
    )
    mwas_parser.add_argument(
        '-v', '--variant', type = str, required = True, metavar = '<str>',
        help = 'Path to the variant file.\nA header line "Chromosome<tab>Position<tab>Population<tab>ID1<tab>ID2 ..." should be present.'
    )
    mwas_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    mwas_parser.add_argument(
        '-c', '--covariate', type = str, required = False, metavar = '<str>',
        help = 'Path to the covariate file.\nA header line "ID<tab>Covariate1<tab>Covariate2 ..." should be present.'
    )
    mwas_parser.add_argument(
        '-t', '--threads', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads.\nDefault: {os.cpu_count()}.'
    )
    mwas_parser.add_argument(
        '--phenotype-column', default = 2, type = int, required = False, metavar = '<int>',
        help = 'Column of phenotype.\nDefault: 2.'
    )
    mwas_parser.add_argument(
        '--mwas', type = str, required = False, metavar = '<str>',
        help = 'Path to the mwas file.\nMajor alleles will be determined in this file rather than in the variant file.\nA header line "Classification<tab>Chromosome<tab>Position<tab>Allele<tab>Beta<tab>SE<tab>P<tab>Significance" should be present.'
    )
    mwas_parser.add_argument(
        '--size', default = 0.50, type = float, required = False, metavar = '<float>',
        help = 'Minimum proportion of samples for mwas.\nDefault: 0.50.'
    )
    mwas_parser.add_argument(
        '--pca-components', default = 5, type = int, required = False, metavar = '<int>',
        help = 'Number of components of GDM to keep.\nDefault: 5.'
    )
    mwas_parser.add_argument(
        '--alpha', default = 0.05, type = float, required = False, metavar = '<float>',
        help = 'The probability of making the wrong decision when the null hypothesis is true.\nDefault: 0.05.'
    )

    # plotMwas parser #
    plotMWAS_parser = subparsers.add_parser(
        'plotMWAS', formatter_class = RawTextHelpFormatter,
        help = 'QQ and Manhattan plots for MWAS results.',
        usage = '%(prog)s [options] -m <MWAS> -o <MWAS>'
    )
    plotMWAS_parser.add_argument(
        '-m', '--mwas', type = str, required = True, metavar = '<str>',
        help = 'Path to the mwas file.\nA header line "Classification<tab>Chromosome<tab>Position<tab>Allele<tab>Beta<tab>SE<tab>P<tab>Significance" should be present.'
    )
    plotMWAS_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Prefix of the output files.\n\"*.qq.pdf\": the qq plot for the mwas.\n\"*.manhattan\": the manhattan plot for the mwas.'
    )
    plotMWAS_parser.add_argument(
        '--variant-annotation', type = str, required = False, metavar = '<str>',
        help = 'Path to the variant annotation file.\nA header line "Chromosome<tab>Position<tab>Annotation ..." should be present.'
    )
    plotMWAS_parser.add_argument(
        '--alpha', default = 0.05, type = float, required = False, metavar = '<float>',
        help = 'The probability of making the wrong decision when the null hypothesis is true.\nDefault: 0.05.'
    )
    plotMWAS_parser.add_argument(
        '--max-negative-log-p', default = 20, type = float, required = False, metavar = '<float>',
        help = 'Maximum value of -log10(P).\nDefault: 20.'
    )

    # checkm2 parser #
    checkm2_parser = subparsers.add_parser(
        'checkm2', formatter_class = RawTextHelpFormatter,
        help = 'Estimate the quality of clusters using CheckM2.',
        usage = '%(prog)s [options] -i <FASTAs> -o <CHECKM2>'
    )
    checkm2_parser.add_argument(
        '-i', '--input', nargs = '+', type = str, required = True, metavar = '<str>',
        help = 'Path to the fasta files.'
    )
    checkm2_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    checkm2_parser.add_argument(
        '-t', '--threads', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    checkm2_parser.add_argument(
        '--checkm2', default = 'checkm2', type = str, required = False, metavar = '<str>',
        help = 'Path to the "checkm2".\nDefault: environmental checkm2.'
    )

    # gtdbtk parser #
    gtdbtk_parser = subparsers.add_parser(
        'gtdbtk', formatter_class = RawTextHelpFormatter,
        help = 'Classify the clusters using GTDB-Tk.',
        usage = '%(prog)s [options] -i <FASTAs> -o <GTDBTk>'
    )
    gtdbtk_parser.add_argument(
        '-i', '--input', nargs = '+', type = str, required = True, metavar = '<str>',
        help = 'Path to the fasta files.'
    )
    gtdbtk_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = '<str>',
        help = 'Path to the output file.'
    )
    gtdbtk_parser.add_argument(
        '--mash-db', default = 'GTDB-Tk.msh', type = str, required = False, metavar = '<str>',
        help = 'Path to the Mash reference sketch database used by GTDB-Tk.\nIf not specified, it will be automatically created as "GTDB-Tk.msh".\nDefault: GTDB-Tk.msh.'
    )
    gtdbtk_parser.add_argument(
        '-c', '--checkm2', type = str, required = False, metavar = '<str>',
        help = 'Path to the file generated by MetaCAT\'s checkm2.\nA header line "Name<tab>Completeness<tab>Contamination ..." should be present.'
    )
    gtdbtk_parser.add_argument(
        '--contamination', default = 0.10, type = float, required = False, metavar = '<float>',
        help = 'Contamination threshold, only available if "--checkm2" is enabled.\nThe value should be from 0 to 1.\nDefault: 0.10.'
    )
    gtdbtk_parser.add_argument(
        '--completeness', default = 0.70, type = float, required = False, metavar = '<float>',
        help = 'Completeness threshold, only available if "--checkm2" is enabled.\nThe value should be from 0 to 1.\nDefault: 0.70.'
    )
    gtdbtk_parser.add_argument(
        '-t', '--threads', default = os.cpu_count(), type = int, required = False, metavar = '<int>',
        help = f'Number of threads.\nThe value should be a positive integer.\nDefault: {os.cpu_count()}.'
    )
    gtdbtk_parser.add_argument(
        '--gtdbtk', default = 'gtdbtk', type = str, required = False, metavar = '<str>',
        help = 'Path to the "gtdbtk".\nDefault: environmental gtdbtk.'
    )
    return parser.parse_args()


def parse_parameters(parameters):
    random = uuid4().hex
    if parameters.command == 'coverage':
        for i, file in enumerate(parameters.bam):
            parameters.bam[i] = os.path.abspath(file)
            isReadable(parameters.bam[i])
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)
        assert parameters.threads_index > 0, '"-ti|--threads-index" should be a positive integer.'
        assert parameters.threads_count > 0, '"-tc|--threads-count" should be a positive integer.'
        assert parameters.min_sequence_length > parameters.trim * 2 + 1, '"-m|--min-sequence-length" should be greater than 2 * "trim" + 1.'
        assert parameters.trim >= 0, '"--trim" should be a non-negative integer.'
        assert parameters.mapq >= 0, '"--mapq" should be a non-negative integer.'
        assert parameters.identity >= 0 and parameters.identity <= 1, '"--identity" should be from 0 to 1.'
        if parameters.long:
            parameters.identity = 0

    elif parameters.command == 'seed':
        parameters.fasta = os.path.abspath(parameters.fasta)
        isReadable(parameters.fasta)
        if parameters.fraggenescan is None:
            parameters.fraggenescan = findBinary('FragGeneScan')
        parameters.fraggenescan = findExecutablePath(parameters.fraggenescan, '--fraggenescan')
        if parameters.hmmpress is None:
            parameters.hmmpress = findBinary('hmmpress')
        parameters.hmmpress = findExecutablePath(parameters.hmmpress, '--hmmpress')
        if parameters.hmmsearch is None:
            parameters.hmmsearch = findBinary('hmmsearch')
        parameters.hmmsearch = findExecutablePath(parameters.hmmsearch, '--hmmsearch')
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)
        assert parameters.threads > 0, '"-t|--threads" should be a positive integer.'

    elif parameters.command == 'cluster':
        parameters.fasta = os.path.abspath(parameters.fasta)
        isReadable(parameters.fasta)
        parameters.seed = os.path.abspath(parameters.seed)
        isReadable(parameters.seed)
        parameters.coverage = os.path.abspath(parameters.coverage)
        isReadable(parameters.coverage)
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(f'{parameters.output}.{random}')
        assert parameters.threads > 0, '"-t|--threads" should be a positive integer.'
        if parameters.min_sequence_length is not None:
            assert parameters.min_sequence_length > 0, '"-m|--min-sequence-length" should be a positive integer.'
        assert parameters.min_cluster > 0, '"--min-cluster" should be a positive integer.'
        assert parameters.min_swdpgmm_clusters >= 2, '"--min-swdpgmm-clusters" should be an integer greater than or equal to 2.'
        assert parameters.min_swdpgmm_sequences >= 2, '"--min-swdpgmm-sequences" should be an integer greater than or equal to 2.'
        assert parameters.max_swdpgmm_iterations > 0, '"--max-swdpgmm-iterations" should be a positive integer.'
        assert parameters.min_pca_variance_ratio > 0 and parameters.min_pca_variance_ratio <= 1, '"--min-pca-variance-ratio" should be from 0 to 1.'
        assert parameters.min_pca_components > 0, '"--min-pca-components" should be a positive integer.'
        assert len(parameters.kmer_frequence_weight) == len(parameters.coverage_weight), 'The numbers of values for "--kmer-frequence-weight" and "--coverage-weight" should be equal.'
        for i, j in zip(parameters.kmer_frequence_weight, parameters.coverage_weight):
            assert i >= 0 and i <= 1, '"--kmer-frequence-weight" should be from 0 to 1.'
            assert j >= 0 and j <= 1, '"--coverage-weight" should be from 0 to 1.'
            assert i + j == 1, 'The sum of each pair of values in "--kmer-frequence-weight" and "--coverage-weight" should be equal to 1.'
        assert parameters.min_kmer_frequence_probability >= 0 and parameters.min_kmer_frequence_probability <= 1, '"--min-kmer-frequence-probability" should be from 0 to 1.'
        assert parameters.min_coverage_probability >= 0 and parameters.min_coverage_probability <= 1, '"--min-coverage-probability" should be from 0 to 1.'
        assert parameters.neighbors > 0, '"--neighbors" should be positive integers.'
        assert parameters.max_seeds > 0, '"--max-seeds" should be positive integers.'

    elif parameters.command == 'cluster2mapping':
        for i, file in enumerate(parameters.fasta):
            parameters.fasta[i] = os.path.abspath(file)
            isReadable(parameters.fasta[i])
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)

    elif parameters.command == 'mapping2cluster':
        parameters.fasta = os.path.abspath(parameters.fasta)
        isReadable(parameters.fasta)
        parameters.mapping = os.path.abspath(parameters.mapping)
        isReadable(parameters.mapping)
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output, pathType = 'directory')
        if parameters.checkm2 is not None:
            parameters.checkm2 = os.path.abspath(parameters.checkm2)
            isReadable(parameters.checkm2)
            assert parameters.contamination >= 0 and parameters.contamination <= 1, '"--contamination" should be from 0 to 1.'
            assert parameters.completeness >= 0 and parameters.completeness <= 1, '"--completeness" should be from 0 to 1.'
        assert parameters.min_cluster > 0, '"--min-cluster" should be a positive integer.'

    elif parameters.command == 'benchmarkGT':
        parameters.ground_truth = os.path.abspath(parameters.ground_truth)
        isReadable(parameters.ground_truth)
        for i, file in enumerate(parameters.mapping):
            parameters.mapping[i] = os.path.abspath(file)
            isReadable(parameters.mapping[i])
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(f'{parameters.output}.{random}')
        if parameters.label is not None:
            assert len(parameters.label) == len(parameters.mapping), 'The numbers of values for "--label" and "--mapping" should be equal.'
        for i in parameters.precision:
            assert i >= 0 and i <= 1, '"--precision" should be from 0 to 1.'
        for i in parameters.recall:
            assert i >= 0 and i <= 1, '"--recall" should be from 0 to 1.'
        isValidColor(parameters.color)
        assert parameters.color_min >= 0 and parameters.color_min <= 1, '"--color-min" should be from 0 to 1.'
        assert parameters.color_max >= 0 and parameters.color_max <= 1, '"--color-max" should be from 0 to 1.'

    elif parameters.command == 'benchmarkRW':
        parameters.input = os.path.abspath(parameters.input)
        isReadable(parameters.input)
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(f'{parameters.output}.{random}')
        for i in parameters.contamination:
            assert i >= 0 and i <= 1, '"--contamination" should be from 0 to 1.'
        for i in parameters.completeness:
            assert i >= 0 and i <= 1, '"--completeness" should be from 0 to 1.'
        if parameters.rows is not None:
            parameters.rows > 0, '"--rows" should be a positive integer.'
        if parameters.columns is not None:
            parameters.columns > 0, '"--columns" should be a positive integer.'
        if (parameters.program is not None) and (parameters.label is not None):
            assert len(parameters.label) == len(parameters.program), 'The numbers of values for "--label" and "--program" should be equal.'
        isValidColor(parameters.color)
        assert parameters.color_min >= 0 and parameters.color_min <= 1, '"--color-min" should be from 0 to 1.'
        assert parameters.color_max >= 0 and parameters.color_max <= 1, '"--color-max" should be from 0 to 1.'

    elif parameters.command == 'indexBam':
        for i, file in enumerate(parameters.input):
            parameters.input[i] = os.path.abspath(file)
            isReadable(parameters.input[i])
        assert parameters.threads > 0, '"-t|--threads" should be a positive integer.'

    elif parameters.command == 'representative':
        for i, file in enumerate(parameters.fasta):
            parameters.fasta[i] = os.path.abspath(file)
            isReadable(parameters.fasta[i])
        isReadable(parameters.checkm2)
        parameters.gtdbtk = os.path.abspath(parameters.gtdbtk)
        isReadable(parameters.gtdbtk)
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(f'{parameters.output}.{random}')
        assert parameters.threads > 0, '"-t|--threads" should be a positive integer.'
        if parameters.mash is None:
            parameters.mash = findBinary('mash')
        parameters.mash = findExecutablePath(parameters.mash, '--mash')
        assert parameters.distance >= 0 and parameters.distance <= 1, '"--distance" should be from 0 to 1.'
        assert parameters.kmer_size > 0, '"--kmer-size" should be a positive integer.'
        assert parameters.sketch_size > 0, '"--sketch-size" should be a positive integer.'
        assert parameters.contamination >= 0 and parameters.contamination <= 1, '"--contamination" should be from 0 to 1.'
        assert parameters.completeness >= 0 and parameters.completeness <= 1, '"--completeness" should be from 0 to 1.'

    elif parameters.command == 'abundance':
        parameters.annotation = os.path.abspath(parameters.annotation)
        isReadable(parameters.annotation)
        for i, file in enumerate(parameters.bam):
            parameters.bam[i] = os.path.abspath(file)
            isReadable(parameters.bam[i])
        assert len(parameters.mapping) == 1 or len(parameters.mapping) == len(parameters.bam), 'The numbers of values for "-m|--mapping" and "-b|--bam" should be equal if the number of values for "-m|--mapping" is not one.'
        for i, file in enumerate(parameters.mapping):
            parameters.mapping[i] = os.path.abspath(file)
            isReadable(parameters.mapping[i])
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)
        assert parameters.threads_index > 0, '"-ti|--threads-index" should be a positive integer.'
        assert parameters.threads_count > 0, '"-tc|--threads-count" should be a positive integer.'
        assert parameters.trim >= 0, '"--trim" should be a non-negative integer.'
        assert parameters.mapq >= 0, '"--mapq" should be a non-negative integer.'
        assert parameters.identity >= 0 and parameters.identity <= 1, '"--identity" should be from 0 to 1.'
        if parameters.long:
            parameters.identity = 0
        assert parameters.min_abundance >= 0 and parameters.min_abundance <= 1, '"--min-abundance" should be from 0 to 1.'

    elif parameters.command == 'abundanceTest':
        parameters.abundance = os.path.abspath(parameters.abundance)
        isReadable(parameters.abundance)
        parameters.group = os.path.abspath(parameters.group)
        isReadable(parameters.group)
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(f'{parameters.output}.{random}')
        assert parameters.coverage >= 0 and parameters.coverage <= 1, '"--coverage" should be from 0 to 1.'
        if parameters.rows is not None:
            parameters.rows > 0, '"--rows" should be a positive integer.'
        if parameters.columns is not None:
            parameters.columns > 0, '"--columns" should be a positive integer.'
        assert parameters.alpha >= 0 and parameters.alpha <= 1, '"--alpha" should be from 0 to 1.'

    elif parameters.command == 'variant':
        for i, file in enumerate(parameters.input):
            parameters.input[i] = os.path.abspath(file)
            isReadable(parameters.input[i])
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)
        assert parameters.threads_call > 0, '"-tc|--threads-call" should be a positive integer.'
        assert parameters.threads_index > 0, '"-ti|--threads-index" should be a positive integer.'
        assert parameters.coverage > 0 and parameters.coverage <= 1, '"--coverage" should be from 0 to 1.'
        assert parameters.mapq >= 0, '"--mapq" should be a non-negative integer.'
        assert parameters.identity >= 0 and parameters.identity <= 1, '"--identity" should be from 0 to 1.'
        assert parameters.depth > 0, '"--depth" should be a positive integer.'

    elif parameters.command == 'mwas':
        parameters.phenotype = os.path.abspath(parameters.phenotype)
        isReadable(parameters.phenotype)
        parameters.mapping = os.path.abspath(parameters.mapping)
        isReadable(parameters.mapping)
        parameters.abundance = os.path.abspath(parameters.abundance)
        isReadable(parameters.abundance)
        parameters.annotation = os.path.abspath(parameters.annotation)
        isReadable(parameters.annotation)
        if parameters.covariate is not None:
            parameters.covariate = os.path.abspath(parameters.covariate)
            isReadable(parameters.covariate)
        parameters.variant = os.path.abspath(parameters.variant)
        isReadable(parameters.variant)
        if parameters.mwas is not None:
            parameters.mwas = os.path.abspath(parameters.mwas)
            isReadable(parameters.mwas)
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)
        assert parameters.threads > 0, '"-t|--threads" should be a positive integer.'
        assert parameters.phenotype_column > 0, '"--phenotype-column" should be a positive integer.'
        assert parameters.size > 0 and parameters.size <= 1, '"--size" should be from 0 to 1.'
        assert parameters.pca_components > 0, '"--pca-components" should be a positive integer.'
        assert parameters.alpha >= 0 and parameters.alpha <= 1, '"--alpha" should be from 0 to 1.'

    elif parameters.command == 'plotMWAS':
        parameters.mwas = os.path.abspath(parameters.mwas)
        isReadable(parameters.mwas)
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(f'{parameters.output}.{random}')
        if parameters.variant_annotation is not None:
            parameters.variant_annotation = os.path.abspath(parameters.variant_annotation)
            isReadable(parameters.variant_annotation)
        assert parameters.alpha >= 0 and parameters.alpha <= 1, '"--alpha" should be from 0 to 1.'
        assert parameters.max_negative_log_p >= 0, '"--max-negative-log-p" should be a non-negative float.'

    elif parameters.command == 'checkm2':
        for i, file in enumerate(parameters.input):
            parameters.input[i] = os.path.abspath(file)
            isReadable(parameters.input[i])
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)
        assert parameters.threads > 0, '"-t|--threads" should be a positive integer.'
        parameters.checkm2 = findExecutablePath(parameters.checkm2, '--checkm2')

    elif parameters.command == 'gtdbtk':
        for i, file in enumerate(parameters.input):
            parameters.input[i] = os.path.abspath(file)
            isReadable(parameters.input[i])
        parameters.output = os.path.abspath(parameters.output)
        isWriteable(parameters.output)
        if parameters.checkm2 is not None:
            parameters.checkm2 = os.path.abspath(parameters.checkm2)
            isReadable(parameters.checkm2)
            assert parameters.contamination >= 0 and parameters.contamination <= 1, '"--contamination" should be from 0 to 1.'
            assert parameters.completeness >= 0 and parameters.completeness <= 1, '"--completeness" should be from 0 to 1.'
        assert parameters.threads > 0, '"-t|--threads" should be a positive integer.'
        parameters.gtdbtk = findExecutablePath(parameters.gtdbtk, '--gtdbtk')
    return None


def main():
    parameters = __init__()
    parse_parameters(parameters)
    if parameters.command == 'coverage':
        coverage.main(parameters)
    elif parameters.command == 'seed':
        seed.main(parameters)
    elif parameters.command == 'cluster':
        cluster.main(parameters)
    elif parameters.command == 'cluster2mapping':
        cluster2mapping.main(parameters)
    elif parameters.command == 'mapping2cluster':
        mapping2cluster.main(parameters)
    elif parameters.command == 'benchmarkGT':
        benchmark_gt.main(parameters)
    elif parameters.command == 'benchmarkRW':
        benchmark_rw.main(parameters)
    elif parameters.command == 'indexBam':
        index_bam.main(parameters)
    elif parameters.command == 'representative':
        representative.main(parameters)
    elif parameters.command == 'abundance':
        abundance.main(parameters)
    elif parameters.command == 'abundanceTest':
        abundance_test.main(parameters)
    elif parameters.command == 'variant':
        variant.main(parameters)
    elif parameters.command == 'mwas':
        mwas.main(parameters)
    elif parameters.command == 'plotMWAS':
        plot_mwas.main(parameters)
    elif parameters.command == 'checkm2':
        checkm2.main(parameters)
    elif parameters.command == 'gtdbtk':
        gtdbtk.main(parameters)
    return None
