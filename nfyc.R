library(synder)
library(ape)
library(dplyr)
library(XLConnect)
library(knitr)
library(ggplot2)
library(reshape2)

dir.create("output", showWarnings=FALSE)

get_in <- function(f){
  fin <- file.path("arendsee", "synder-nfyc", "archive", f) 
  if(!file.exists(fin)){
    stop("Cannot find file ", fin)
  }
  fin
}
get_out <- function(f){
  if(!dir.exists("output")){
    dir.create("output")
  }
  file.path("output", f)
}

tree <- read.tree(get_in('tree.newick'))
pdf(get_out('tree.pdf'))
plot(tree)
dev.off()

collect_NFYC_data <- function(
    syn_file,
    fgff_file,
    tgff_file,
    gene_map_file,
    tblastn_file,
    blastp_file = NULL,
    trans  = "d",
    k      = 0L,
    r      = 0
){
  extractNameFromAttr <- function(x){
    sub(".*Name=([^;]+).*", "\\1", x)
  }

  # load synteny map
  syn  <- synder::as_synmap(syn_file)

  # load focal GFF
  fgff <- synder:::.tidy_gff(fgff_file, prefix="f", fromAttr=extractNameFromAttr)

  # load the transcribed intervals target GFF
  tgff <- synder:::.tidy_gff(tgff_file, prefix="t", fromAttr=extractNameFromAttr)
  # FIXME: do I have to restrict this to mRNA?
  tgff <- tgff[tgff$ttype == "mRNA", ]

  # load map with name (NF-YC*), gene (locus ID), to transcript (model ID)
  nfyc <- read.table(
    gene_map_file,
    col.names=c("fgene_name", "fseqid", "fprotein"),
    stringsAsFactors=FALSE
  )
  # a map from fprot_id to fseqid
  fmap <- nfyc[, c(3,2)]

  # map focal genes to target search intervals 
  srcres <- synder::search(syn, synder::as_gff(fgff_file), k=k,r=r,trans=trans)
  srcres <- synder:::.tidy_searchResult(srcres, fromAttr=extractNameFromAttr)

  tblastn_si_map <- synder::make_tblastn_si_map(tblastn_file, srcres, fmap=fmap)
  tblastn_gene_map <- synder::make_tblastn_gene_map(tblastn_file, tgff)

  if(is.null(blastp_file)){
    blastp_raw <- NULL 
    blastp_map <- NULL
  } else {
    x <- make_blastp_map(
      blastp_file=blastp_file,
      tgff=tgff,
      fgff=fgff,
      srcres=srcres,
      fmap=fmap
    )
    blastp_raw <- x$blastraw
    blastp_map <- x$blastmap
  }

  feats <- featureMap(srcres=srcres, fgff=fgff, tgff=tgff)

  tblastn_raw <- tblastn_si_map$blastraw
  tblastn_si_map <- tblastn_si_map$blastmap
  tblastn_gene_map <- tblastn_gene_map$blastmap

  synder:::.namedlist(
    syn, tgff, fgff, nfyc, srcres,
    tblastn_raw, tblastn_si_map, tblastn_gene_map,
    blastp_raw, blastp_map, feats)
}

files <- list(
  Arabidopsis_lyrata=list(syn="at-vs-al.tab", gff="al.gff", tblastn="al.blast.tab", blastp="yc-vs-al.blast.tab"),
  Capsella_rubella=list(syn="at-vs-cr.tab", gff="cr.gff", tblastn="cr.blast.tab"),
  Brassica_rapa=list(syn="at-vs-br.tab", gff="br.gff", tblastn="br.blast.tab"),
  Eutrema_salsugineum=list(syn="at-vs-es.tab", gff="es.gff", tblastn="es.blast.tab")
)

if(file.exists(get_out("results.Rds"))){
  results <- readRDS(get_out("results.Rds"))
} else {
  results <- lapply(files, function(x){
    collect_NFYC_data(
      fgff_file     = get_in("nfyc.gff"),
      gene_map_file = get_in("nfyc-map.tab"),
      syn_file      = get_in(x$syn),
      tgff_file     = get_in(x$gff),
      tblastn_file  = get_in(x$tblastn),
      blastp_file   = if(!is.null(x$blastp)) {get_in(x$blastp)} else {NULL},
      k=0L, r=0, trans="d"
    )
  })
  saveRDS(results, get_out("results.Rds"))
}

### ===========================================================================

tabulate_from <- function(xs, key){
  ys <- lapply(xs, function(x) x[[key]])
  for(grp in names(ys)){
    ys[[grp]]$group <- grp
  }
  do.call(rbind, ys) %>% synder:::.remove_rownames()
}

nuniq <- function(x){ length(unique(x)) }

show_vec <- function(xs, nmax=Inf) {
  if(nuniq(xs) <= nmax) {
    paste(unique(xs), collapse=",")
  } else {
    "lots"
  }
}

### ==========================================================================
# joined GFFs
tgffs <- tabulate_from(results, "tgff")
tgffs$group=NULL

srcress <- tabulate_from(results, "srcres")
srcress$si_fwidth <- srcress$si_fstop - srcress$si_fstart + 1
srcress$si_twidth <- srcress$si_tstop - srcress$si_tstart + 1

fstrand_map <- tabulate_from(results, "fgff") %>%
  dplyr::select(fseqid, fstrand) %>%
  dplyr::distinct()

### ===========================================================================
# results based solely on synder
feature_maps <- tabulate_from(results, "feats")
per_feature <- feature_maps %>%
  dplyr::filter(strands_agree) %>%
  dplyr::group_by(group, fseqid) %>%
  dplyr::summarize(n_synder_hits = nuniq(tseqid), synder_hits=show_vec(tseqid))


### ===========================================================================
# results based solely on tblastn
tblastns <- tabulate_from(results, "tblastn_raw") %>%
  dplyr::filter(tn_evalue < 0.001)
tblastns$id <- 1:nrow(tblastns)

tblastn_against_proteins <- synder:::.overlaps(
  x=tblastns, xid="tchr", xa="tn_sstart", xb="tn_send",
  y=tgffs,    yid="tchr", ya="tstart",    yb="tstop"
)

# get the total number of tBLASTn hits against target proteins 
total_tblastn_against_proteins <-
  tblastn_against_proteins %>% 
  dplyr::group_by(group, fseqid) %>%
  dplyr::summarize(n_tblastn_hits = nuniq(tseqid), tblastn_hits=show_vec(tseqid))

# get the number of signifcant hits that do not overlap a target protein
nongenic_tblastn <- tblastns[-tblastn_against_proteins$id, ]
genic_tblastn <- tblastn_against_proteins %>% dplyr::select(-xid, -yid, -id)

### ===========================================================================
tblastn_si_maps <-  tabulate_from(results, "tblastn_si_map")
tblastn_si_maps$cset <- as.character(tblastn_si_maps$cset)
tblastn_si_maps <- dplyr::filter(tblastn_si_maps, tn_evalue < 0.001)

# Find overlaps with tgff protein features, this resoves the above exon issue
tblastn_synder_hits <-
  # get between tBLASTn hits and target genes
  synder:::.overlaps(
    x=tblastn_si_maps, xid="tchr", xa="tn_tstart", xb="tn_tstop",
    y=tgffs, yid="tchr", ya="tstart", yb="tstop"
  ) %>%
  # extract rows where the protein coded for by query gene i has a tBLASTn hit
  # within a search interval inferred for the same gene. 
  dplyr::filter(fseqid == si_fseqid) %>%
  # get the strand of the BLAST hit (-1,-2,-3 are '-' strand, 1,2,3 are '+')
  dplyr::mutate(hit_strand = ifelse(tn_tframe > 0, '+', '-')) %>%
  # get strand info for the focal species
  merge(fstrand_map) %>%
  # determine whether the strands of the synder intervals and the target genes agree
  dplyr::mutate(strand_agreement = 
    (orientation == '+' & (fstrand == tstrand)) |
    (orientation == '-' & (fstrand != tstrand))
  ) %>%
  dplyr::select(
    group, fseqid, fchr, fstrand, tchr, fprotein, tn_tstart, tn_tstop, prot_fstart,
    prot_fstop, hit_strand, tn_tframe, tn_evalue, tn_ppos, tstart, tstop, tseqid,
    tsource, ttype, tstrand, cset, si_fstart, si_fstop, orientation, si_score,
    l_flag, r_flag, inbetween, strand_agreement
  ) %>%
  list(
    blast_offstrand = subset(., hit_strand != tstrand),
    synder_offstrand = subset(., !strand_agreement),
    # require the tBLASTn hit and target gene are on the same strand
    # require strand agreement
    norm = subset(., strand_agreement & (hit_strand == tstrand))
  )

tblastn_synder_summary <- tblastn_synder_hits$norm %>%
  dplyr::group_by(group, fseqid) %>%
  dplyr::summarize(n_tblastn_synder_targets = nuniq(tseqid),
                   tblastn_synder_targets=show_vec(tseqid))

# summarize each group
group_summary <- dplyr::group_by(tblastn_si_maps, group) %>%
  dplyr::summarize(
    fseqids = length(fseqid),
    sets = nuniq(cset),
    agreements = mean(fseqid == si_fseqid),
    tchrs = nuniq(tchr)
  )

# summarize each cset
cset_summary <- dplyr::group_by(tblastn_si_maps, cset) %>%
  dplyr::summarize(
    fseqids = nuniq(fseqid),
    tseqids = nuniq(si_fseqid),
    min_eval = min(tn_evalue),
    fseqids_ref = show_vec(fseqid),
    si_fseqids = show_vec(si_fseqid)
  )
cset_summary_b <- dplyr::group_by(tblastn_si_maps, group, si_fseqid) %>%
  dplyr::summarize(csets = nuniq(cset))


# ### ===========================================================================

count_summaries <-
  merge(
    total_tblastn_against_proteins,
    per_feature,
    by=c("group","fseqid"),
    all=TRUE
  ) %>%
  merge(
    tblastn_synder_summary,
    by=c("group","fseqid"),
    all=TRUE
  ) %>% dplyr::select(
    group, fseqid,
    n_tblastn_synder_targets, n_synder_hits, n_tblastn_hits,
    tblastn_synder_targets, synder_hits, tblastn_hits 
  )
count_summaries <- merge(dplyr::distinct(results[[1]]$nfyc[, 1:2]), count_summaries)
count_summaries[is.na(count_summaries)] <- 0
# assert that there is one row for each focal gene / target species pair
stopifnot(nuniq(results[[1]]$nfyc$fseqid) * nuniq(names(results)) == nrow(count_summaries))

### ===========================================================================
# Make an excel spreadsheet

key_sheet <- as.data.frame(matrix(ncol=2, byrow=TRUE, c(
  "fseqid",        "focal feature ID",
  "ftype",         "focal feature type",
  "fsource",       "focal feature source",
  "fchr",          "focal chromosome or scaffold",
  "fstart",        "focal genomic start position",
  "fstop",         "focal genomic stop position",
  "fstrand",       "focal strand for this feature",
  "tseqid",        "target feature ID",
  "ttype",         "target feature type",
  "tsource",       "target feature source",
  "tchr",          "target chromosome or scaffold",
  "tstart",        "target genomic start position",
  "tstop",         "target genomic stop position",
  "tstrand",       "target strand for this feature",
  "si_fstart",     "start position on the focal (query) side of a search interval link",
  "si_fstop",      "stop position on the focal (query) side of a search interval link",
  "si_tstart",     "start position on the target side of a search interval link",
  "si_tstop",      "stop position on the target side of a search interval link",
  "orientation",   "orientation of the target syntenic interval relative to the query genome ('+' means they are on strand)",
  "si_score",      "synder score for the search interval",
  "cset",          "an ID for the contiguous set this search interval was based on",
  "l_flag",        "left-side synder flag",
  "r_flag",        "right-side synder flag",
  "inbetween",     "is the query contained entirely inbetween syntenic links?",
  "strands_agree", "does the strand of the target feature match the expected strand given the search interval?",
  "group",         "target species"
)))
colnames(key_sheet) <- c("column", "description")

tree <- read.tree(get_in('tree.newick'))
png(get_out('tree.png'))
plot(tree)
dev.off()

# increase memory to handle large files
options(java.parameters = "-Xmx1024m")
file.remove(get_out("nfyc.xlsx"))
wb <- XLConnect::loadWorkbook(get_out("nfyc.xlsx"), create=TRUE)
XLConnect::createSheet(wb, "Key")
XLConnect::writeWorksheet(wb, data=key_sheet, sheet="Key")
XLConnect::createSheet(wb, "Tree")
XLConnect::createName(
  wb,
  name    = 'Tree',
  formula = paste('Tree', XLConnect::idx2cref(c(1, 1)), sep="!")
)
XLConnect::addImage(
  wb,
  filename     = get_out("tree.png"),
  name         = 'Tree',
  originalSize = TRUE
)

addTable <- function(w, d, sheetName){
  XLConnect::createSheet(w, sheetName)
  XLConnect::writeWorksheet(w, data=d, sheet=sheetName)
}

addTable(wb, count_summaries, "Hits")
addTable(wb, tblastn_synder_hits$norm, "tblastn_synder")
addTable(wb, tblastn_synder_hits$blast_offstrand, "tblastn_synder_blast_offstrand")
addTable(wb, tblastn_synder_hits$synder_offstrand, "tblastn_synder_synder_offstrand")
addTable(wb, srcress, "search_intervals")
addTable(wb, feature_maps, "synder_featureMap")
addTable(wb, genic_tblastn, "tblastn_genic")
addTable(wb, nongenic_tblastn, "tblastn_non-genic")
XLConnect::saveWorkbook(wb)

write(knitr::kable(count_summaries[, 1:6], format='latex'), get_out("count_summaries.tex"))

m <- reshape2::melt(
  count_summaries[, 1:6],
  id.vars=c("fseqid", "fgene_name", "group")
) %>%
  dplyr::rename(target_species = group, count=value, group=variable) %>%
  dplyr::mutate(
    gene_name = factor(fgene_name, levels=paste0("NF-YC", 1:13)),
    target_species = factor(
      target_species,
      levels=c(
        "Arabidopsis_lyrata",
        "Capsella_rubella",
        "Brassica_rapa",
        "Eutrema_salsugineum"
      )
    )
  )
m$group <- factor(
  as.character(m$group),
  levels=c(
    "n_synder_hits",
    "n_tblastn_hits",
    "n_tblastn_synder_targets"
  )
)

pdf(get_out("count_fig.pdf"))

legend_labels <- c("synder", "tBLASTn", "tBLASTn+synder")
ggplot(m) +
  geom_point(aes(x=gene_name, y=count, shape=group, color=group), size=2) +
  geom_hline(aes(yintercept=0), alpha=0.2) +
  geom_hline(aes(yintercept=13), alpha=0.2) +
  facet_grid(target_species ~ .) +
  scale_shape_manual(values = c(4, 20, 0), labels=legend_labels) +
  scale_color_discrete(labels=legend_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=270, hjust=0, vjust=1),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_y_continuous(minor_breaks = seq(0 , 20, 1), breaks = seq(0, 20, 5))

dev.off()
```
