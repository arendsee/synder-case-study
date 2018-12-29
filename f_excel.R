# excel
#  * gather all tables and some figures into a single Excel document 
#  * add annotation for each table
#  * add annotation for each column in each table
# TODO: the metadata should all be included in the Data Package and I should
# just be able to extract it from the JSON file

make_spreadsheet <- function(ds, out_file, tree){

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

  png('.tree.png')
    plot(tree)
  dev.off()

  # increase memory to handle large files
  options(java.parameters = "-Xmx1024m")
  if(file.exists(out_file)){
    file.remove(out_file)
  }
  wb <- XLConnect::loadWorkbook(out_file, create=TRUE)
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
    filename     = ".tree.png",
    name         = 'Tree',
    originalSize = TRUE
  )
  file.remove('.tree.png') # this won't be needed anymore

  addTable <- function(w, d, sheetName){
    XLConnect::createSheet(w, sheetName)
    XLConnect::writeWorksheet(w, data=d, sheet=sheetName)
  }

  addTable(wb, ds$count_summaries, "Hits")
  addTable(wb, ds$tblastn_synder_hits$norm, "tblastn_synder")
  addTable(wb, ds$tblastn_synder_hits$blast_offstrand, "tblastn_synder_blast_offstrand")
  addTable(wb, ds$tblastn_synder_hits$synder_offstrand, "tblastn_synder_synder_offstrand")
  addTable(wb, ds$srcress, "search_intervals")
  addTable(wb, ds$feature_maps, "synder_featureMap")
  addTable(wb, ds$genic_tblastn, "tblastn_genic")
  addTable(wb, ds$nongenic_tblastn, "tblastn_non-genic")
  XLConnect::saveWorkbook(wb)

}
