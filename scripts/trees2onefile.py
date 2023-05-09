import os
import re



def replace(string: str) -> str:
    '''
    Remove floats and : from tree string. Shorten names to 4 letters. Remove _paralog and number
    '''
    string = re.sub(r'\d+\.\d+', '', string)
    string = re.sub(r':', '', string)
    string = re.sub(r'\w{4}\w+', lambda x: x.group(), string)
    string = string.replace('Inner', '')
    string = string.replace(';', '')
    string = re.sub(r'\)\d+', ')', string)
    string = re.sub(r'\w{4}_paralog\d+', lambda x: x.group()[:4], string)

    return string


def save_trees2onefile(input_dir: dict, output: str) -> None:
    '''
    Writes all sequences from the input file to a single fasta file
    '''
    paths_list = os.listdir(input_dir)
    paths_list = [input_dir+'/'+path for path in paths_list]
    trees_nwk = []
    for t in paths_list:
        with open(t, 'r') as f:
            tree = f.read()
        tree = replace(tree)
        tree = tree.strip()
        trees_nwk.append(tree)
    # sort trees from the longest to the shortest
    trees_nwk.sort(key=lambda x: len(x), reverse=True)
    with open(output, 'w') as one_file:
        for i, t in enumerate(trees_nwk):
            t = t.replace('-', '')
            if i < len(trees_nwk)-1:
                one_file.write(t+';\n')
            else:
                one_file.write(t)

if __name__ == '__main__':
    save_trees2onefile(snakemake.input[0], snakemake.output[0])