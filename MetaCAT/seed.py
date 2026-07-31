import gzip
import os
from datetime import datetime
from math import ceil
from pathlib import Path
from uuid import uuid4

from .fraggenescan import Fraggenescan
from .hmm import dump_database, get_valid_hits, runHmmsearch
from .marker import markerHash, nProfiles


def parse_sequence_id(sequence_id):
    '''
    sequence_id:
        FragGeneScan: sequence_start_end_strand
        Prodigal: sequence_number
    '''
    return sequence_id.rsplit('_', maxsplit = 3)[0]


def get_seeds(hitGenerator, file):
    sequence2profileScore = dict()
    for sequence, profile, score in hitGenerator:
        if score > sequence2profileScore.get(sequence, ('', -1e10))[1]:
            sequence2profileScore[sequence] = (profile, score)
    sequence2profiles = dict()
    for sequence, profileScore in sequence2profileScore.items():
        sequence2profiles.setdefault(parse_sequence_id(sequence), set()).add(markerHash[profileScore[0]])
    openFile = gzip.open(file, mode = 'wt', compresslevel = 9)
    openFile.write('Sequence ID\tMarker IDs\n')
    for sequence in sorted(sequence2profiles.keys()):
        profiles = ' '.join(str(i) for i in sorted(sequence2profiles[sequence]))
        openFile.write(f'{sequence}\t{profiles}\n')
    openFile.close()
    return None


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Identifying protein sequences.', flush = True)
    protein = Path(parameters.temp).joinpath(uuid4().hex)
    fraggenescan = Fraggenescan(parameters.fraggenescan, parameters.threads, 'protein', parameters.temp)
    fraggenescan.main(parameters.fasta, str(protein))
    if protein.stat().st_size:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Initializing the database of markers.', flush = True)
        profile2score, database_files = dump_database(
            str(Path(__file__).with_name('markers.gz')), ceil(nProfiles / parameters.threads), parameters.temp
        )
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Mapping markers to protein sequences.', flush = True)
        hit_file = runHmmsearch(parameters.hmmsearch, database_files, str(protein), parameters.threads, parameters.temp)
        for i in database_files:
            os.remove(i)
        hit_generator = get_valid_hits(profile2score, hit_file)
        get_seeds(hit_generator, parameters.output)
        os.remove(hit_file)
    else:
        print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Cannot predict any protein sequences.', flush = True)
    protein.unlink()
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} Finished.', flush = True)
    return None
