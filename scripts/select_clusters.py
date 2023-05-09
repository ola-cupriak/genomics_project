import os
import random
import sys


def clusters2dict(mmseq_file: str) -> dict:
    """
    Saves clusters from mmseqs2 to dictionary.
    """
    with open(mmseq_file  , 'r') as f:
        contet = f.read()
    clusters = {}
    key = ''
    for part in contet.split('>'):
        if part == '': continue
        if len(part.strip()) == 6 or len(part.strip()) == 10:
            key = part.strip()
            clusters[part.strip()] = []
        else:
            clusters[key].append(part.strip())
    return clusters


def get_species_name(description: str) -> str:
    """
    Returns species name from description.
    """
    return description.split('OS=')[1].split('OX=')[0].strip()


def split_desc_seq(protein: str) -> str:
    """
    Returns description and sequence from protein.
    """
    return protein.split('\n')[0], protein.split('\n')[1]


def ref_spec_names(species_file: str) -> set:
    """
    Returns set of species names from species_file.
    """
    with open(species_file, 'r') as f:
        species_names = f.readlines()
        species_names = [name.split('\t')[1].strip() for name in species_names]
    return set(species_names)


def check_names(ref_name: str, name: str) -> bool:
    """
    Checks if name is the same as ref_name.
    """
    ref_name = ref_name.split(' ')
    name = name.split(' ')
    name = [el.replace('(', '') for el in name]
    name = [el.replace(')', '') for el in name]
    # check
    if ref_name[0] != name[0]: return False
    if ref_name[1] != name[1]: return False
    if len(ref_name) > 2 and len(name) > 2:
        for ref_el in ref_name[2:]:
            if ref_el not in name: return False
    return True


def change_names(clusters: dict, species_names: set) -> None:
    """
    Changes species names in clusters to names from species_names.
    """
    change_spec_name = dict()
    for clust in clusters.values():
        to_remove = []
        for i, protein in enumerate(clust):
            try:
                description, seq = split_desc_seq(protein)
                spec_name = get_species_name(description)
                if spec_name not in change_spec_name.keys():
                    for ref_name in species_names:
                        if check_names(ref_name, spec_name):
                            change_spec_name[spec_name] = ref_name
                            break
                clust[i] = (change_spec_name[spec_name], seq)
            except:
                to_remove.append(protein)
        for protein in to_remove:
            clust.remove(protein)
    
def select_clusters(clusters: dict, min_cluster_size: int = 30, max_cluster_size: int = 100000, 
            max_repeat: float = 1/30, orthological: bool = False, choose_random: bool = False,
            excluded_species=set()) -> dict:
    """
    Selects clusters from clusters dictionary.
    Args:
        clusters (dict): dictionary with clusters
        min_cluster_size (int): minimal size of cluster
        max_cluster_size (int): maximal size of cluster
        max_repeat (float): maximal repeat of one species in cluster
        orthological (bool): if True, only orthological clusters are selected
        choose_random (bool): if True, random seq from species is selected
        excluded_species (set): set of species to exclude
    Returns:
        selected_clusters (dict): dictionary with selected clusters
    """
    selected_clusters = {}
    for name, clust in clusters.items():
        if len(clust) < min_cluster_size: continue
        if max_cluster_size and len(clust) > max_cluster_size: continue
        species_in_cluster = [element[0] for element in clust]
        if excluded_species:
            if len([i for i in species_in_cluster if i in excluded_species]) > 0: continue
        
        if not choose_random:
            most_common = max(set(species_in_cluster), key=species_in_cluster.count)
            if species_in_cluster.count(most_common) / len(clust) > max_repeat:
                continue

        if orthological and not choose_random:
            used_species = set()
            if_next = False
            for spec, _ in clust:
                if spec in used_species:
                    if_next = True
                    break
                used_species.add(spec)
            if if_next: continue
        
        if orthological and choose_random:
            species_dict = {}
            new_clust = []
            for spec, seq in clust:
                if spec in species_dict.keys():
                    species_dict[spec].append(seq)
                else:
                    species_dict[spec] = [seq]
            for spec, seqs in species_dict.items():
                # random choose of one seq from seqs
                random_seq = random.choice(seqs)
                new_clust.append((spec, random_seq))
            if len(new_clust) < min_cluster_size: continue
            clust = new_clust

        selected_clusters[name] = clust
    print(f'select_clusters: {len(selected_clusters)} clusters have been selected.')
    return selected_clusters


def check_unique_names(cluster: list) -> list:
    '''
    Checks name uniqueness.
    If the name is not unique, adds number to the end of name.
    '''
    names = {}
    for i, protein in enumerate(cluster):
        species, seq = protein
        if species in names.keys():
            names[species] += 1
            species += f'_paralog{names[species]+1}'
        else:
            names[species] = 0
        protein = (species, seq)
        cluster[i] = protein
    return cluster


def cluster2fasta(cluster: list, file: str) -> None:
    """
    Saves cluster to fasta file.
    """
    with open(file, 'w') as f:
        for protein in cluster:
            species, seq = protein
            species = species.replace(' ', '_')
            f.write(f'>{species}\n{seq}\n')


def create_dirs(path: str) -> None:
    """
    Creates directories for output files.
    """
    elements = path.split('/')
    for i in range(len(elements)):
        tocheck = '/'.join(elements[0:i+1])
        if not os.path.exists(tocheck):
            try:
                os.mkdir(tocheck)
            except:
                sys.stderr = print('ERROR: An invalid path to the directory for output files is given.')
                sys.exit(0)



def run(mmseq_file: str, species_file: str, output_dir: str,  
        min_cluster_size: int, max_cluster_size: int, max_repeat: float, 
        orthological: bool, choose_random: bool, excluded_species: set) -> None:
    """
    Runs selecting clusters.
    Args:
        mmseq_file (str): path to mmseq file
        species_file (str): path to file with species names
        output_dir (str): path to directory for output files
        min_cluster_size (int): minimal size of cluster
        max_cluster_size (int): maximal size of cluster
        max_repeat (float): maximal repeat of one species in cluster
        orthological (bool): if True, only orthological clusters are selected
        choose_random (bool): if True, random seq from species is selected
        excluded_species (set): set of species to exclude
    """
    create_dirs(output_dir)
    clusters = clusters2dict(mmseq_file)
    species_names = ref_spec_names(species_file)
    change_names(clusters, species_names)
    selected_clusters = select_clusters(clusters, min_cluster_size, max_cluster_size, 
                                        max_repeat, orthological, choose_random, excluded_species)
    for name, cluster in selected_clusters.items():
        cluster = check_unique_names(cluster)
        cluster2fasta(cluster, f'{output_dir}/{name}.fasta')


if __name__ == '__main__':
    min_cluster_size = snakemake.config['selecting_clusters']['min_cluster_size']
    max_cluster_size = snakemake.config['selecting_clusters']['max_cluster_size']
    max_repeat = snakemake.config['selecting_clusters']['max_repeat']
    orthological = snakemake.config['selecting_clusters']['orthological']
    choose_random = snakemake.config['selecting_clusters']['choose_random']
    excluded_species = snakemake.config['selecting_clusters']['excluded_species']
    excluded_species = set(excluded_species)
    run(snakemake.input[0], snakemake.config['species_file'], snakemake.output[0],
        min_cluster_size, max_cluster_size, max_repeat, orthological, choose_random, excluded_species)