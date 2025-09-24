# table compare function

tab_comp <- function(ca_tab, sql_tab)
{
  query_command <- paste0(" SELECT * FROM ", sql_tab)
  sql_table <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names()
  
  if("substrata" %in% names(sql_table))
  {
    sql_table <-  sql_table %>% rename(stratum = substrata)
  }
  if("strata" %in% names(sql_table))
  {
    sql_table <-  sql_table %>% rename(subarea = strata)
  }
  if("sd_length" %in% names(sql_table))
  {
    sql_table <-  sql_table %>% rename(sd_length_sk = sd_length)
  }
  
  tab_compare <- full_join(ca_tab, sql_table)
  print(paste(" ca tab rows    :", dim(ca_tab)[1]))
  if("year" %in% names(ca_tab)) {
  print(paste(" ca tab w/o 2021:", dim(ca_tab %>% dplyr::filter(year != 2021))[1]))}
  print(paste(" sk tab rows    :", dim(sql_table)[1]))
  if("subarea" %in% names(ca_tab)) {
    print(paste(" ca tab w/o 2021, sub 0:", dim(ca_tab %>% dplyr::filter(subarea != 0))[1]))}
  
  return(tab_compare)
}