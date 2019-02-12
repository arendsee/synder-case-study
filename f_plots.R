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
    dplyr::mutate(xid = "synder")

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
    dplyr::mutate(xid = variable)

  # make list of blast scores
  .load <- function(x, isgen){
    merge(result[[x]], prot2name) %>%
    dplyr::mutate(genic=isgen, score=as_levalue(tn_evalue)) %>%
    dplyr::select(fid=gene, group, score, genic) %>%
    dplyr::mutate(xid = "blast")
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

  gene_levels <- paste0("NF-YC", 1:13)
  m_links$fid <- factor(m_links$fid, levels=gene_levels)
  m_synder$fid <- factor(m_synder$fid, levels=gene_levels)
  m_blast$fid <- factor(m_blast$fid, levels=gene_levels)

  # abbreviate names (e.g Zea mays -> Z. mays) and order them relative
  # to alphabetic ordering
  clean_names <- function(x, ord){
    x <- sub("([A-Z])[a-z]+_(.*)", "\\1. \\2", x)
    x <- factor(x, levels=sort(unique(x))[ord]) 
    x
  }
  m_links$group <- clean_names(m_links$group, c(1,3,2,4))
  m_synder$group <- clean_names(m_synder$group, c(1,3,2,4))
  m_blast$group <- clean_names(m_blast$group, c(1,3,2,4))

  ggplot() +
    geom_line(data=m_links, aes(x=xid, y=score, group=linkid), alpha=0.2) +
    geom_point(data=m_links,aes(x=xid, y=score), size=0.2) +
    facet_grid(group ~ fid) +
    ylab("Normalized Score") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=270, hjust=0, vjust=1),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(0, "lines"),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(size = 6),
      strip.text.y = element_text(size = 8, face="italic"),
      strip.background = element_blank(),
      axis.title.x = element_blank()
    ) +
    scale_shape_manual(values = c(1, 4)) + # 4 : x
    scale_color_manual(values = c("darkgray", "red")) +
    geom_point(data=m_blast, aes(x=xid, y=score, shape=genic, color=genic)) +
    geom_point(data=m_synder, aes(x=xid, y=score), color="blue")
}
