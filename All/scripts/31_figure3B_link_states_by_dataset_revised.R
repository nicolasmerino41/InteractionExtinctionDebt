## Original entry point retained; fresh inputs from the established R loader,
## exact removal analysis and all rendering in Julia/CairoMakie.
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','.UTF-8')
source('All/scripts/00_dataset_loaders_and_helpers_all.R')
out <- 'All/outputs/31_hidden_support_main_figures'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pp<-list();ss<-list()
for(ds in all_dataset_names){
 t<-get_dataset_site_tables(ds)
 ints<-t$empirical_site_interactions %>% distinct(site,consumer,resource)
 cooc<-bind_rows(t$cooc_triples,ints) %>% distinct(site,consumer,resource)
 p<-cooc %>% count(consumer,resource,name='full_n') %>% left_join(ints %>% count(consumer,resource,name='full_K'),by=c('consumer','resource')) %>% mutate(full_K=replace_na(full_K,0L),dataset=ds)
 stopifnot(all(p$full_K<=p$full_n))
 pp[[ds]]<-p;ss[[ds]]<-data.frame(dataset=ds,n_sites=n_distinct(cooc$site))
}
write.csv(bind_rows(pp),file.path(out,'31_3B_fresh_pair_support.csv'),row.names=FALSE)
write.csv(bind_rows(ss),file.path(out,'31_3B_fresh_dataset_sizes.csv'),row.names=FALSE)
status<-system2('julia',c('--startup-file=no',shQuote('All/scripts/31_figure3B_revised.jl'),shQuote(normalizePath('.'))))
if(status!=0)stop('Julia Figure 31 rendering failed')
