# collate
#  * Make cummulative tables from the species-wise data in "collect"

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

collate_results <- function(results){
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

  synder:::.namedlist(count_summaries, tblastn_synder_hits, srcress, feature_maps, genic_tblastn, nongenic_tblastn)
}
