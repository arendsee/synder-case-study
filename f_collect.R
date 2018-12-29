# collect
#  * link together syntenic intervals, synder seach intervals, genome features,
#    and BLAST results
#  * sanitize the results into consistently named, tidy data.frames 

collect_results <- function(files, fgff_file, gene_map_file){
  if(file.exists(get_out("results.Rds"))){
    results <- readRDS(get_out("results.Rds"))
  } else {
    results <- lapply(files, function(x){
      .collect_result(
        fgff_file     = get_in(fgff_file),
        gene_map_file = get_in(gene_map_file),
        syn_file      = get_in(x$syn),
        tgff_file     = get_in(x$gff),
        tblastn_file  = get_in(x$tblastn),
        blastp_file   = if(!is.null(x$blastp)) {get_in(x$blastp)} else {NULL},
        k=0L, r=0, trans="d"
      )
    })
    saveRDS(results, get_out("results.Rds"))
  }
  results
}

.collect_result <- function(
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
