import gzip
import os
from datetime import datetime
from math import ceil

from .fraggenescan import runFraggenescan
from .hmm import check_database, create_database, get_valid_hits, read_score_file, runHmmsearch
from .marker import nProfiles


def parse_sequence_id(sequence_id):
    '''
    sequence_id:
        FragGeneScan: sequence_start_end_strand
        Prodigal: sequence_number
    '''
    return sequence_id.rsplit('_', maxsplit = 3)[0]


def get_seeds(hit_generator, output_file):
    sequence2profile_score = dict()
    for sequence, profile, score in hit_generator:
        if score > sequence2profile_score.get(sequence, ('', -1e10))[1]:
            sequence2profile_score[sequence] = (profile, score)
    profile2sequences = dict()
    for sequence, (profile, _) in sequence2profile_score.items():
        profile2sequences.setdefault(profile, list()).append(parse_sequence_id(sequence))
    open_file = gzip.open(output_file, mode = 'wt', compresslevel = 9)
    for profile in sorted(profile2sequences.keys()):
        open_file.write('\t'.join([profile] + sorted(set(profile2sequences[profile]))) + '\n')
    open_file.close()
    return None


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Checking the database of markers.', flush = True)
    marker_file = os.path.join(os.path.dirname(__file__), 'markers.gz')
    n = ceil(nProfiles / ceil(nProfiles / os.cpu_count()))
    check_marker, database_files = check_database(marker_file, n)
    if not check_marker:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Creating the database of markers.', flush = True)
        create_database(parameters.hmmpress, marker_file, n)
        database_files = check_database(marker_file, n)[1]
    profile2score = read_score_file(database_files[0])

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Identifying protein sequences.', flush = True)
    protein = runFraggenescan(parameters.fraggenescan, parameters.fasta, parameters.threads)
    if os.path.getsize(protein):
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Mapping markers to protein sequences.', flush = True)
        hit_file = runHmmsearch(parameters.hmmsearch, database_files[1 : ], protein, parameters.threads)
        hit_generator = get_valid_hits(profile2score, hit_file)
        get_seeds(hit_generator, parameters.output)
        os.remove(hit_file)
    else:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Cannot predict any protein sequences.', flush = True)
    os.remove(protein)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
