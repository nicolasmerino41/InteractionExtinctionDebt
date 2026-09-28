# Export fresh inputs with the existing loader; calculations and figures are Julia/Makie.
args<-commandArgs(trailingOnly=TRUE)
if(length(args))setwd(args[1])
if(.Platform$OS.type=='windows')Sys.setlocale('LC_CTYPE','.UTF-8')
source('All/scripts/00_dataset_loaders_and_helpers_all.R')
# Local correction only: drop unnamed matrix rows/columns BEFORE data.frame
# name repair turns the Gottin_HP unnamed consumer into V1.
original_clean43<-clean_web_matrix
clean_web_matrix<-function(web){
 if(!is.null(colnames(web)))web<-web[,!is.na(colnames(web))&nzchar(colnames(web)),drop=FALSE]
 if(!is.null(rownames(web)))web<-web[!is.na(rownames(web))&nzchar(rownames(web)),,drop=FALSE]
 original_clean43(web)
}
out<-'All/outputs/43 main figures/inputs';dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(ds in all_dataset_names){
 message('Preparing ',ds)
 st<-get_dataset_site_tables(ds)
 ints<-st$empirical_site_interactions %>% distinct(site,consumer,resource)
 co<-bind_rows(st$cooc_triples,ints) %>% distinct(site,consumer,resource)
 occ<-st$occupancy %>% distinct(site,species,trophic_level)
 write.csv(ints,file.path(out,paste0(ds,'_interactions.csv')),row.names=FALSE)
 write.csv(co,file.path(out,paste0(ds,'_cooccurrences.csv')),row.names=FALSE)
 write.csv(occ,file.path(out,paste0(ds,'_occupancy.csv')),row.names=FALSE)
}
