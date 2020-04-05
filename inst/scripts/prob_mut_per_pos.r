# setwd("~/micropopgen/exp/2020/today3")
# setwd("~/micropopgen/src/phylomut/")
# devtools::document()

gene_phylomut_or <- function(aln_file, tree_file, Map){
  # Read data
  aln <- phangorn::read.phyDat(aln_file, format = 'fasta',
                               type = 'USER',
                               levels = c('a', 'c', 'g', 't', '-'))
  tree <- ape::read.tree(tree_file)
  Map <- Map %>%
    dplyr::filter(ID %in% tree$tip.label)
  
  # Ancestral reconstruction
  ml.fit <- phangorn::pml(tree, aln)
  ml.fit <- phangorn::optim.pml(ml.fit, model = "GTR", control = phangorn::pml.control(trace=1))
  anc <- phangorn::ancestral.pml(ml.fit, 'bayes')
  
  # test chn
  ids <- Map$ID[Map$Group == "CHN"]
  mono_clades <- find_mono_clades(tre = tree, ids = ids, return = "outgroup")
  
  # Mono probs
  Probs.chn <- mono_clades %>%
    purrr::map(~mut_probs(tree = tree,
                          aln = aln,
                          anc = anc,
                          root_id = .$outgroup,
                          tips = .$tips)) %>%
    purrr::map(~.x %>% dplyr::select(prob, position)) %>%
    purrr::map(~.x[attr(aln,"index"),] %>%
                 dplyr::mutate(pattern = position,
                               position = 1:length(attr(aln, "index"))))
  
  # test chn
  ids <- Map$ID[Map$Group == "USA"]
  mono_clades <- find_mono_clades(tre = tree, ids = ids, return = "outgroup")
  
  # Mono probs
  Probs.usa <- mono_clades %>%
    purrr::map(~mut_probs(tree = tree,
                          aln = aln,
                          anc = anc,
                          root_id = .$outgroup,
                          tips = .$tips)) %>%
    purrr::map(~.x %>% dplyr::select(prob, position)) %>%
    purrr::map(~.x[attr(aln,"index"),] %>%
                 dplyr::mutate(pattern = position,
                               position = 1:length(attr(aln, "index"))))
  
  # Calculate log2(OR)
  or <- log2(apply(Probs.chn %>%
                     purrr::map(~.x$prob) %>%
                     dplyr::bind_cols() %>%
                     as.matrix, 1, median)) - log2(apply(Probs.usa %>%
                                                           purrr::map(~.x$prob) %>%
                                                           dplyr::bind_cols() %>%
                                                           as.matrix, 1, median))
  or <- tibble::tibble(position = 1:length(or),
                       log2OR = or)
  
  return(or)
}


library(magrittr)
# library(phangorn)
library(phylomut)
library(ggplot2)
# library(tidyverse)

# args <- list(tree = "example/435590.9.peg.366.baseml.tre",
#              aln = "example/435590.9.peg.366.aln.fasta",
#              map = "Bacteroides_vulgatus_57955.map.txt")

args <- list(trees = "Bacteroides_vulgatus_57955/gene_trees/",
             alns = "Bacteroides_vulgatus_57955/gene_alns/",
             map = "Bacteroides_vulgatus_57955/Bacteroides_vulgatus_57955.map.txt",
             rer = "Bacteroides_vulgatus_57955/Bacteroides_vulgatus_57955.rer.fdr.txt",
             outdir = "output")

rerperms <- readr::read_tsv(args$rer) %>%
  dplyr::filter(!is.na(FDR))
Map <- readr::read_tsv(args$map)

genes <- rerperms$gene_id[rerperms$FDR < 0.01]

Res <- NULL
for(g in genes){
  # g <- genes[1]
  cat(g, "\n")
  
  aln_file <- file.path(args$alns,paste0(g, ".aln.fasta"))
  tree_file <- file.path(args$trees,paste0(g, ".baseml.tre"))
  
  res <- gene_phylomut_or(aln_file = aln_file, tree_file = tree_file, Map = Map)
  res$gene_id <- g
  
  Res <- Res %>%
    dplyr::bind_rows(res)
}

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

readr::write_tsv(Res, "sig_genes_or.txt")

### MOVE PLOTTING TO LOOP
Res %>%
  split(.$gene_id) %>%
  purrr::map(function(d, outdir){
    g <- unique(d$gene_id)
    p1 <- ggplot(d, aes(x=position, y = log2OR)) +
      # geom_bar(stat = "identity") +
      geom_point() +
      ggtitle(label = g) +
      AMOR::theme_blackbox()
    # p1
    filename <- paste0(g,".phylomut.svg")
    filename <- file.path(outdir, filename)
    ggsave(filename, p1, width = 16, height = 9)
  }, outdir = args$outdir)


# p1 <- ggplot(res, aes(x=position, y = log2OR)) +
#   # geom_bar(stat = "identity") +
#   geom_point() +
#   AMOR::theme_blackbox()
# p1

