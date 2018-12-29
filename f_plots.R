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
    scale_y_continuous(minor_breaks = seq(0 , 20, 1), breaks = seq(0, 20, 5))
}
