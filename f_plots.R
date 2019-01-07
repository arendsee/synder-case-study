# plots
#  * plots !!

nameless_plot <- function(count_summaries){

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
      ),
      group = factor(
        as.character(group),
        levels=c(
          "n_synder_hits",
          "n_tblastn_hits",
          "n_tblastn_synder_targets"
        )
      )
    )

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
    scale_y_continuous(minor_breaks = seq(0, 20, 1), breaks = seq(0, 20, 5))
}

wonky_plot <- function(){
  # map from gene name (e.g. NF-YC4) to focal protein ID
  prot2name <- results[[1]]$nfyc %>%
    dplyr::select(gene=fgene_name, fprotein)
  # map from gene name (e.g. NF-YC4) to focal gene ID
  gene2name <- results[[1]]$nfyc %>%
    dplyr::select(gene=fgene_name, fseqid)
  # transform evalue
  as_levalue <- function(x) (-1) * log10(x)

  # make list of synder scores
  m_synder <-
    result$srcress %>%
    merge(gene2name, by="fseqid") %>%
    dplyr::select(fid=gene, group, score=si_score, cset=cset) %>%
    dplyr::distinct() %>%
    dplyr::mutate(xid = paste0(fid, "_synder"))

  # make table of edges between tBLASTn and synder hits
  m_links <-
    result$tblastn_si_maps %>%
    dplyr::mutate(
      levalue=as_levalue(tn_evalue),
      linkid=1:n()
    ) %>%
    merge(prot2name, by="fprotein") %>%
    dplyr::select(fid=gene, group, synder=si_score, blast=levalue, cset, linkid) %>%
    melt(id.vars=c("fid", "group", "cset", "linkid")) %>%
    dplyr::rename(score=value) %>%
    dplyr::mutate(xid = paste(fid, variable, sep="_"))

  # make list of blast scores
  .load <- function(x, isgen){
    merge(result[[x]], prot2name) %>%
    dplyr::mutate(genic=isgen, score=as_levalue(tn_evalue)) %>%
    dplyr::select(fid=gene, group, score, genic) %>%
    dplyr::mutate(xid = paste0(fid, "_blast"))
  }
  m_blast <- rbind(.load("genic_tblastn", "Genic"), .load("nongenic_tblastn", "Non-Genic"))

  normalize_blast_score <- function(x){
    with(m_blast, (x - mean(score)) / sd(score))
  }
  normalize_synder_score <- function(x){
    with(result$srcress, (log(x) - mean(log(si_score))) / sd(log(si_score)))
  }

  m_links$score[m_links$variable == "synder"] <- normalize_synder_score(m_links$score[m_links$variable == "synder"])
  m_links$score[m_links$variable == "blast"] <- normalize_blast_score(m_links$score[m_links$variable == "blast"])
  m_synder$score <- normalize_synder_score(m_synder$score)
  m_blast$score <- normalize_blast_score(m_blast$score)

  # normalize scores 

  ggplot() +
    geom_line(data=m_links, aes(x=xid, y=score, group=linkid), alpha=0.3) +
    geom_point(data=m_links,aes(x=xid, y=score), size=0.2) +
    facet_grid(group ~ .) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=270, hjust=0, vjust=1),
      legend.title = element_blank(),
      legend.position = "bottom"
    ) +
    scale_shape_manual(values = c(4, 0)) +
    geom_point(data=m_blast, aes(x=xid, y=score, color=genic, shape=genic)) +
    geom_point(data=m_synder, aes(x=xid, y=score), color="blue")
}
