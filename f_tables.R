# tables
#  * tables !!

make_tables <- function(ds){
  knitr::kable(ds$count_summaries[, 1:6], format='latex') %>%
    write(get_out("count_summaries.tex"))
}
