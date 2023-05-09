import os


def proteomes2onefile(input_file: str, proteome_files: list, 
                      output_file: str, excluded_species: list=[]):
    '''
    Saves all proteomes to one fasta file.
    '''
    with open(input_file, 'r') as f:
        entries = f.readlines()
        ids = [e.split('\t')[0].strip() for e in entries]
        species = [e.split('\t')[1].strip() for e in entries]   
    # Saves all proteomes to one fasta file
    if excluded_species:
        directory = proteome_files[0].split('/')[:-1]
        directory = ('/').join(directory)
        to_rm = [directory+'/'+id+'.fasta' for i, id in enumerate(ids) if species[i] in excluded_species]
        for file in to_rm:
            proteome_files.remove(file)
    with open(output_file, 'w') as f:
        for proteome in proteome_files:
            try:
                with open(proteome, 'r') as file:
                    f.write(file.read())
            except FileNotFoundError:
                print(f'ERROR: File {proteome} not found!')


if __name__ == '__main__':
    proteomes2onefile(snakemake.input[0], snakemake.params[0], snakemake.output[0], snakemake.config['save_proteomes']['to_remove'])