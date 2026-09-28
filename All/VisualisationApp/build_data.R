# Run from the repository root:
# Rscript All/VisualisationApp/build_data.R
# Optional arguments: repository root, output folder.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else getwd()
out <- if(length(args)>1) args[2] else file.path(root,'All','VisualisationApp')
out <- normalizePath(out, mustWork=FALSE)
setwd(root)
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','.UTF-8')
source('All/scripts/00_dataset_loaders_and_helpers_all.R')
stopifnot(requireNamespace('jsonlite',quietly=TRUE))
all_data <- list(); audit <- list()
for(dataset in all_dataset_names){
 message('Exporting ',dataset)
 st <- get_dataset_site_tables(dataset)
 d <- unique(st$empirical_site_interactions[c('site','consumer','resource')])
 names(d)<-c('s','c','r');d[]<-lapply(d,as.character)
 ss<-sort(unique(as.character(c(st$cooc_triples$site,st$occupancy$site,d$s))))
 cc<-sort(unique(d$c));rr<-sort(unique(d$r))
 p<-unique(d[c('c','r')]);p<-p[order(p$c,p$r),]
 links<-lapply(seq_len(nrow(p)),function(i) list(match(p$c[i],cc)-1L,match(p$r[i],rr)-1L,sort(unique(match(d$s[d$c==p$c[i]&d$r==p$r[i]],ss)-1L))))
 M<-matrix(0L,length(ss),nrow(p))
 for(i in seq_along(links)) M[links[[i]][[3]]+1L,i]<-1L
 stopifnot(sum(M)==nrow(d),all(colSums(M)>0))
 cluster_order<-function(x) if(nrow(x)<2) seq_len(nrow(x))-1L else hclust(dist(x,method='binary'),method='average')$order-1L
 all_data[[dataset]]<-list(s=ss,c=cc,r=rr,l=links,so=cluster_order(M),lo=cluster_order(t(M)))
 audit[[dataset]]<-data.frame(dataset,sites=nrow(M),consumers=length(cc),resources=length(rr),links=ncol(M),occurrences=sum(M))
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(all_data,file.path(out,'datasets.json'),auto_unbox=TRUE)
write.csv(do.call(rbind,audit),file.path(out,'data_checks.csv'),row.names=FALSE)
print(do.call(rbind,audit))
