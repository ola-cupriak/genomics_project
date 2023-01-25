configfile: "configuration.yaml"


with open(config['species_file'], 'r') as file:
    lines = file.readlines()
    PROTEOMES = [line.split('\t')[0] for line in lines]


def get_input(wildcards):
    input_list = []
    if config['results']['supertree']:
        if config['bootstrap']:
            if config['selecting_clusters']['choose_random']:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'random'+'_b'+str(config['bootstrap_thresh'])+"_supertree.nwk")
            else:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'_b'+str(config['bootstrap_thresh'])+"_supertree.nwk")
        else:
            if config['selecting_clusters']['choose_random']:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'random'+'_nb'+"_supertree.nwk")
            else:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'_nb'+"_supertree.nwk")
    if config['results']['consensus']:
        if config['bootstrap']:
            if config['selecting_clusters']['choose_random']:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'random'+'_b'+str(config['bootstrap_thresh'])+"_consensus.nwk")
            else:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'_b'+str(config['bootstrap_thresh'])+"_consensus.nwk")
        else:
            if config['selecting_clusters']['choose_random']:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'random'+'_nb'+"_consensus.nwk")
            else:
                input_list.append("results/final/"+config['names']['subset']+"_"+config['names']['type']+'_nb'+"_consensus.nwk")
    return input_list


rule all:
    input:
        get_input


rule download_proteomes:
    input:
        config['species_file']
    output:
        expand("results/proteomes_raw/{proteome}.fasta", proteome=PROTEOMES)
    script:
        "scripts/download_proteomes.py"


rule save_to_onefile:
    input:
        config['species_file'],
        expand("results/proteomes_raw/{proteome}.fasta", proteome=PROTEOMES)
    output:
        "results/proteomes/proteomes_{subset}.fasta"
    params:
        expand("results/proteomes_raw/{proteome}.fasta", proteome=PROTEOMES)
    script:
        "scripts/save_proteomes.py"


rule run_mmseq2:
    input:
        "results/proteomes/proteomes_{subset}.fasta"
    output:
        "results/mmseq_res/mmseq_result_{subset}_all_seqs.fasta"
    params:
        c = config['clustering']['c']
    shell:
        "mmseqs easy-cluster {input} results/mmseq_res/mmseq_result_{wildcards.subset} tmp -c {params.c}"


rule select_clusters:
    input:
        "results/mmseq_res/mmseq_result_{subset}_all_seqs.fasta"
    output:
        directory("results/selected_clusters/{subset}_{type}")
    script:
        "scripts/select_clusters.py"


rule build_trees:
    input:
        "results/selected_clusters/{subset}_{type}"
    output:
        directory("results/msa_results/{subset}_{type}"),
        directory("results/nj_trees/{subset}_{type}")
    script:
        "scripts/build_trees.py"


rule bootstrap:
    input:
        "results/msa_results/{subset}_{type}",
        "results/nj_trees/{subset}_{type}"
    output:
        directory("results/nj_trees_bootstrap/{subset}_{type}")
    script:
        "scripts/bootstrap.R"


def check_bootstrap(wildcards):
    if config['bootstrap']:
        if config['selecting_clusters']['choose_random']:
            return ["results/nj_trees_bootstrap/"+config['names']['subset']+"_"+config['names']['type']+'random']
        else:
            return ["results/nj_trees_bootstrap/"+config['names']['subset']+"_"+config['names']['type']]
    else:
        if config['selecting_clusters']['choose_random']:
            return ["results/nj_trees/"+config['names']['subset']+"_"+config['names']['type']+'random']
        else:
            return ["results/nj_trees/"+config['names']['subset']+"_"+config['names']['type']]


rule trees2onefile:
    input:
        check_bootstrap
    output:
        "results/final/{subset}_{type}_{bootstrap}_trees1file.nwk"
    script:
        "scripts/trees2onefile.py"


rule build_supertree:
    input:
        "results/final/{subset}_{type}_{bootstrap}_trees1file.nwk"
    output:
        "results/final/{subset}_{type}_{bootstrap}_supertree.nwk"
    shell:
        """
        sed 's/;//' {input} > {input}_temp
        var=$(fasturec -G {input}_temp -Z | tail -1 | cut -d"-" -f2 | cut -d" " -f1)
        head -n 1 $var | cut -d' ' -f2 > {output}
        sed -i 's/$/;/' {output}
        rm $var
        rm {input}_temp
        """

rule build_consensus:
    input:
        "results/final/{subset}_{type}_{bootstrap}_trees1file.nwk"
    output:
        "results/final/{subset}_{type}_{bootstrap}_consensus.nwk"
    script:
        "scripts/get_consensus.R"
