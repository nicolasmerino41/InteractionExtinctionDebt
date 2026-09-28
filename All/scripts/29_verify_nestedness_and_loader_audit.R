if(.Platform$OS.type=='windows')Sys.setlocale('LC_CTYPE','.UTF-8')
library(vegan)
s<-read.csv('All/outputs/42.1 nestedness/global_nestedness/global_nestedness.csv')
for(ds in s$dataset){
 d<-read.csv(file.path('All/outputs/42',ds,'42_incidence.csv'),colClasses='character')
 M<-unclass(table(d$site,paste(d$consumer,d$resource,sep='___')))>0
 r<-vegan::nestednodf(M,order=TRUE)$statistic
 stopifnot(abs(r['N.rows']-s$site_NODF[s$dataset==ds])<1e-8,
           abs(r['NODF']-s$whole_matrix_NODF[s$dataset==ds])<1e-8)
}
cat('All ten Julia NODF results match vegan::nestednodf\n')
library(bipartite)
d<-read.csv('Gottin/raw-data/host_para_all_interactions.csv',sep=';')
w<-bipartite::frame2webs(d,varnames=c('Genus.Species','P1.Species.Genus','Site','P1.cells'))
for(i in seq_along(w)){
 a<-w[[i]];blank<-which(colnames(a)=='')
 if(length(blank)>0 && any(a[,blank,drop=FALSE]>0,na.rm=TRUE)){
 cat('Site',names(w)[i],': blank consumer column with nonzero entries\n')
 print(a[rowSums(a[,blank,drop=FALSE],na.rm=TRUE)>0,blank,drop=FALSE])
 cat('Column names after as.data.frame: ',paste(names(as.data.frame(a)),collapse=', '),'\n')
 }
}
