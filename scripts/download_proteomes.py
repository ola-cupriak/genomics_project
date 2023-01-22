import os
import requests


def download_proteome(entery: str, output_path: str):
    url = f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3A{entery}%29%29'
    with open(output_path, 'w') as f:
        f.write(requests.get(url).text)


def download_all(species_file: str, output_dir: str):
    output_dir = output_dir.split('/')[:-1]
    output_dir = '/'.join(output_dir)
    with open(species_file, 'r') as file:
        lines = file.readlines()
        proteomes_list = [f'{output_dir}/'+line.split('\t')[0]+'.fasta' for line in lines]
    # Download
    entries = [path.split('/')[-1].split('.') for path in proteomes_list]
    entries = ['.'.join(entry[:-1]) for entry in entries]
    for i, entry in enumerate(entries):
        download_proteome(entry.strip(), proteomes_list[i])



if __name__ == '__main__':
    download_all(snakemake.input[0], snakemake.output[0])