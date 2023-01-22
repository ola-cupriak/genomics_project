from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import Applications
from joblib import Parallel, delayed
import multiprocessing
import os
import sys



def create_dirs(path: str):
    elements = path.split('/')
    for i in range(len(elements)):
        tocheck = '/'.join(elements[0:i+1])
        if not os.path.exists(tocheck):
            try:
                os.mkdir(tocheck)
            except:
                sys.stderr = print('ERROR: An invalid path to the directory for output files is given.')
                sys.exit(0)
    

def run_parallel(func, input_dict, output_dir):
    create_dirs(output_dir)
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(func)(object, name, output_dir) for name, object in input_dict.items())
    # make dict
    results_dict = dict()
    for result in results:
        results_dict[result[0]] = result[1]
    return results_dict


def get_multialignment(input_file: str, cluster_name: str, output_dir: str):
    '''Funkcja przy pomocy programu clustalo tworzy multiuliniowienie 
    na podstawie sekwencji białkowych przekazanych w pliku wejściowym.'''
    # Tworzenie pliku wejściowego
    clustalw_cline = Applications.ClustalwCommandline('clustalo', infile=input_file,
                                                    outfile=f'{output_dir}/MSA_{cluster_name}.aln')
    clustalw_cline()
    proteins_msa = AlignIO.read(f'{output_dir}/MSA_{cluster_name}.aln', format='clustal')
    
    return (cluster_name, proteins_msa)


def tree_construction(MSA, name, output_dir):
    '''Funkcja na podstawie multiuliniowienia generuje drzewa filogenetyczne 
    przy pmetody neighbor-joining dla podanej macierzy blosum. 
    Rysuje je oraz zapisuje w plikach newick.'''
    cal = DistanceCalculator()
    constr = DistanceTreeConstructor()
    proteins = cal.get_distance(MSA)
    
    njtree = constr.nj(proteins)
    Phylo.write(njtree, f'{output_dir}/njtree_{name}.nwk', 'newick')

    return (name, f'{output_dir}/njtree_{name}.nwk')


def run(input_dir: str, msa_dir: str, trees_dir: str):
    paths = os.listdir(input_dir)
    paths = [path for path in paths if path[0] != '.']
    cluster_names = [input_file.split('.')[:-1] for input_file in paths]
    cluster_names = ['.'.join(name) for name in cluster_names]
    paths = {cluster_name: input_dir + '/' + path for cluster_name, path in zip(cluster_names, paths)}
    msa_dict = run_parallel(get_multialignment, paths, msa_dir)
    run_parallel(tree_construction, msa_dict, trees_dir)


if __name__ == '__main__':
    run(snakemake.input[0], snakemake.output[0], snakemake.output[1])