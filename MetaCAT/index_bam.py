from datetime import datetime

from .bam import indexBam


def main(parameters):
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Indexing all bam files.', flush = True)
    indexBam(parameters.bam, parameters.threads)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
    return None
