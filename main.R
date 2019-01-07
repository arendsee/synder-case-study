# main
#  * pass filenamess to the functions that do the real work
#  * determine what is written into the final output

library(synder)
library(ape)
library(dplyr)
library(XLConnect)
library(knitr)
library(ggplot2)
library(reshape2)

source('f_collate.R')
source('f_collect.R')
source('f_excel.R')
source('f_io.R')
source('f_plots.R')
source('f_tables.R')

dir.create("output", showWarnings=FALSE)

results <- collect_results(
  files = list(
    Arabidopsis_lyrata  = list(syn="at-vs-al.tab", gff="al.gff", tblastn="al.blast.tab", blastp="yc-vs-al.blast.tab"),
    Capsella_rubella    = list(syn="at-vs-cr.tab", gff="cr.gff", tblastn="cr.blast.tab"),
    Brassica_rapa       = list(syn="at-vs-br.tab", gff="br.gff", tblastn="br.blast.tab"),
    Eutrema_salsugineum = list(syn="at-vs-es.tab", gff="es.gff", tblastn="es.blast.tab")
  ),
  fgff_file="nfyc.gff",
  gene_map_file="nfyc-map.tab",
  trans  = "d",
  k      = 0L,
  r      = 1e-4
)

result <- collate_results(results) 

tree <- read.tree(get_in('tree.newick'))

pdf(get_out('tree.pdf'))
  plot(tree)
dev.off()

make_spreadsheet(
  ds        = result,
  out_file  = get_out("nfyc.xlsx"),
  tree      = tree
)

make_tables(result)

pdf(get_out("wonky.pdf"))
  wonky_plot()
dev.off()

pdf(get_out("count_fig.pdf"))
  nameless_plot(result$count_summaries)
dev.off()
