# Run: Rscript All/scripts/42/run42.R
# Optional first argument overrides the repository root.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else getwd()
setwd(root)
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','.UTF-8')
source('All/scripts/00_dataset_loaders_and_helpers_all.R')
out <- 'All/outputs/42';dir.create(out,recursive=TRUE,showWarnings=FALSE)
data <- list()
for(dataset in all_dataset_names){
 message('Preparing ',dataset)
 st<-get_dataset_site_tables(dataset)
 ints<-distinct(st$empirical_site_interactions,site,consumer,resource)
 co<-distinct(st$cooc_triples,site,consumer,resource)
 ints[]<-lapply(ints,as.character);co[]<-lapply(co,as.character)
 sites<-sort(unique(c(ints$site,co$site,as.character(st$occupancy$site))))
 pairs<-distinct(ints,consumer,resource)
 observed<-allowed<-vector('list',nrow(pairs))
 for(j in seq_len(nrow(pairs))){
  observed[[j]]<-match(ints$site[ints$consumer==pairs$consumer[j]&ints$resource==pairs$resource[j]],sites)-1L
  allowed[[j]]<-match(co$site[co$consumer==pairs$consumer[j]&co$resource==pairs$resource[j]],sites)-1L
  stopifnot(all(observed[[j]] %in% allowed[[j]]))
 }
 data[[dataset]]<-list(sites=sites,consumer=pairs$consumer,resource=pairs$resource,observed=observed,allowed=allowed)
}
jsonlite::write_json(data,file.path(out,'42_inputs.json'),auto_unbox=TRUE)
python<-Sys.getenv('IED_PYTHON',unset='C:/Users/MM-1/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe')
if(!file.exists(python))python<-Sys.which('python')
status<-system2(python,c(shQuote('All/scripts/42/engine42.py'),shQuote(normalizePath(out))))
if(status!=0)stop('Analysis engine failed')
source('All/scripts/42/plots42.R')
