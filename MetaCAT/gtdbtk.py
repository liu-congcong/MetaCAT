import os
from shutil import rmtree
from subprocess import DEVNULL, PIPE, STDOUT, run
from uuid import uuid4

from .checkm2 import readCheckm2File


def createBatchFile(inputFiles, clusters):
    outputFile = uuid4().hex
    openFile = open(outputFile, 'w')
    for inputFile in inputFiles:
        cluster = os.path.splitext(os.path.basename(inputFile))[0]
        if clusters is None:
            openFile.write(f'{inputFile}\t{cluster}\n')
        elif cluster in clusters:
            openFile.write(f'{inputFile}\t{cluster}\n')
    openFile.close()
    return outputFile


def runGtdbtk(gtdbtk, threads, inputFile, outputFile):
    temp = uuid4().hex
    command = [gtdbtk, 'classify_wf', '--batchfile', inputFile, '--out_dir', temp, '--cpus', str(threads), '--force']
    # --mash_db #
    completedProcess = run([gtdbtk, 'classify_wf', '-h'], stdout = PIPE, stderr = STDOUT)
    assert not completedProcess.returncode, 'An error has occurred while running GTDB-Tk.'
    for line in completedProcess.stdout.decode('utf-8', errors = 'ignore').splitlines():
        if line.lstrip().startswith('--mash_db'):
            command.extend(['--mash_db', f'{temp}.msh'])
            break
    completedProcess = run(command, stdout = None, stderr = None)
    assert not completedProcess.returncode, 'An error has occurred while running GTDB-Tk.'
    if os.access(f'{temp}.msh', os.R_OK):
        os.remove(f'{temp}.msh')
    open4w = open(outputFile, 'w')
    open4w.write('user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings\n')
    for i in ('bac120', 'ar53'):
        file = os.path.join(temp, f'gtdbtk.{i}.summary.tsv')
        if os.access(file, os.R_OK):
            open4r = open(file, 'r')
            open4r.readline()
            for line in open4r:
                open4w.write(line)
            open4r.close()
    open4w.close()
    rmtree(temp)
    return None


def main(parameters):
    if parameters.checkm2:
        clusters = readCheckm2File(parameters.checkm2, parameters.contamination, parameters.completeness)
    else:
        clusters = None
    batchFile = createBatchFile(parameters.fasta, clusters)
    runGtdbtk(parameters.gtdbtk, parameters.threads, batchFile, parameters.output)
    os.remove(batchFile)
    return None
