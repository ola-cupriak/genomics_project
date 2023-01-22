library(ape)
library(phangorn)


bootstrap <- function(msa_dir, njt_dir, msa_path, njt_path){
  njt_path <- paste(njt_dir, njt_path, sep="/")
  msa_path <- paste(msa_dir, msa_path, sep="/")
  tree <- ape::read.tree(njt_path)
  msa <- phangorn::read.phyDat(msa_path, format='clustal')
  fit <- phangorn::pml(tree, msa)
  tryCatch({
    bs <- bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE, mc.cores=8, control = pml.control(trace=0))
    tb <- plotBS(fit$tree, bs, type="p")
    tb.node.label <- tb$node.label[!is.na(tb$node.label)]
    mean <- mean(tb.node.label)
    if(mean>=50){
      return(njt_path)
    }
  }, error=function(e){})
}

create_dir <- function(directory){
  dirs <- strsplit(directory, '/')
  for(i in 1:length(dirs[[1]])){
    to_check <- paste(dirs[[1]][1:i], collapse="/")
    if(!dir.exists(to_check)){
      dir.create(to_check)
    }
  }
}

run <- function(msa_dir, njt_dir){
  nj_paths <- list.files(njt_dir)
  msa_paths <- list.files(msa_dir)
  bs_trees <- mapply(bootstrap, msa_dir=msa_dir, njt_dir=njt_dir, msa_path=msa_paths, njt_path=nj_paths)
  bs_trees <- Filter(Negate(is.null), bs_trees)
  return(bs_trees)
}

outdir <- snakemake@output[[1]]
res <- run(snakemake@input[[1]], snakemake@input[[2]])
create_dir(outdir)
for(file in res){
  print(file)
  file.copy(file, outdir)
}