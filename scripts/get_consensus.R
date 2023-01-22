library(ape)

get_consensus <- function(input, output, p) {
  trees <- ape::read.tree(input)
  consensus <- ape::consensus(trees, p=p, rooted = FALSE)
  ape::write.tree(consensus, output)
}

get_consensus(snakemake@input[[1]], snakemake@output[[1]], snakemake@config[['get_consensus']]['p'])