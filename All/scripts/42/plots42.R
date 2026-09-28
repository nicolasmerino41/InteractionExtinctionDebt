library(ggplot2);library(dplyr);library(tidyr)

out<-'All/outputs/42'
s<-read.csv(file.path(out,'42_dataset_summary.csv'))
g<-read.csv(file.path(out,'42_null_comparisons.csv'))
z<-read.csv(file.path(out,'42_site_importance.csv'))
tr<-read.csv(file.path(out,'42_null_traces.csv'))
rem<-read.csv(file.path(out,'42_removal_diagnostics.csv'))

pair<-bind_rows(lapply(s$dataset,function(ds){
  x<-read.csv(file.path(out,ds,'42_pairwise_site_metrics.csv'))
  x$site1<-as.character(x$site1)
  x$site2<-as.character(x$site2)
  x$dataset<-ds
  x
}))

pair_long<-bind_rows(
  transmute(filter(pair,included_containment=='True'),
            dataset,metric='Containment',value=containment),
  transmute(filter(pair,pmin(richness1,richness2)>0),
            dataset,metric='Turnover',value=turnover),
  transmute(filter(pair,pmin(richness1,richness2)>0),
            dataset,metric='Nestedness',value=nestedness_resultant)
)

pair_long$metric<-factor(
  pair_long$metric,
  levels=c('Containment','Turnover','Nestedness')
)

g$assessment<-ifelse(
  g$diagnostics_ok=='True',
  'Diagnostics passed',
  'Limited / invariant null'
)

z$flag<-ifelse(
  z$diagnostics_ok!='True',
  'Unresolved',
  ifelse(!is.na(z$q_bh)&z$q_bh<.05,
         'Outside expectation',
         'No clear departure')
)

theme_set(
  theme_bw(base_size=11)+
    theme(
      panel.grid.minor=element_blank(),
      strip.background=element_rect(fill='grey95'),
      legend.position='bottom',
      plot.title=element_text(face='bold')
    )
)

cols<-c(Margins='#0072B2',Opportunity='#D55E00')

save_fig<-function(p,path,w=14,h=7){
  for(label in c('subtitle','caption'))
    if(!is.null(p$labels[[label]]))
      p$labels[[label]]<-paste(
        unlist(lapply(
          strsplit(p$labels[[label]],'\n',fixed=TRUE)[[1]],
          function(x)strwrap(x,width=floor(w*12))
        )),
        collapse='\n'
      )
  
  ggsave(paste0(path,'.png'),p,width=w,height=h,dpi=220,bg='white')
  ggsave(paste0(path,'.pdf'),p,width=w,height=h,bg='white')
}


# FIGURE 42A
p1fun<-function(d,fac=TRUE){
  p<-ggplot(d,aes(metric,value,fill=metric))+
    geom_boxplot(
      width=.55,
      outlier.shape=NA,
      linewidth=.4
    )+
    stat_summary(
      fun=mean,
      geom='point',
      shape=23,
      fill='white',
      size=2
    )+
    scale_fill_manual(
      values=c('#009E73','#0072B2','#CC79A7'),
      guide='none'
    )+
    coord_cartesian(ylim=c(0,1))+
    labs(
      x=NULL,
      y='Pairwise value',
      title='Are sites nested?'
    )+
    theme(
      axis.text.x=element_text(angle=25,hjust=1)
    )
  
  if(fac)
    p<-p+facet_wrap(~dataset,ncol=5)
  
  p
}


p2fun<-function(d,fac=TRUE){
  d<-filter(d,metric=='containment')
  
  p<-ggplot(d,aes(x=model,y=null_mean,colour=model))+
    geom_linerange(aes(ymin=null_low,ymax=null_high),linewidth=1)+
    geom_point(aes(shape=assessment),size=3)+
    geom_point(aes(y=observed),shape=4,colour='black',size=3,stroke=1)+
    scale_colour_manual(values=cols,guide='none')+
    scale_shape_manual(
      values=c(
        'Diagnostics passed'=16,
        'Limited / invariant null'=1
      )
    )+
    labs(
      x=NULL,
      y='Mean containment',
      title='42B - Is site containment unusual?',
      subtitle='Coloured point and line: null mean and central 95% range; black x: observed. Panel y-scales differ.',
      caption='Margins: site richness and interaction support fixed. Opportunity: also restricted to valid co-occurrence sites.\nOpen circles flag failed diagnostics or invariant statistics; an apparently narrow interval is not proof of a well-explored null.'
    )
  
  if(fac)
    p<-p+facet_wrap(~dataset,ncol=5,scales='free_y')
  
  p
}


p3fun<-function(d,fac=TRUE){
  p<-ggplot(
    d,
    aes(
      richness,
      100*adjusted_importance,
      colour=model,
      shape=flag,
      size=unique_links
    )
  )+
    geom_hline(yintercept=0,colour='grey65')+
    geom_point(alpha=.7)+
    scale_colour_manual(values=cols)+
    scale_shape_manual(
      values=c(
        'No clear departure'=16,
        'Outside expectation'=17,
        'Unresolved'=1
      )
    )+
    scale_size_continuous(range=c(1.5,4.5))+
    labs(
      x='Observed interactions at site',
      y='Importance minus null mean (percentage points)',
      title='42C - Which sites matter beyond their richness?',
      subtitle='Focal-site excess link loss at approximately 50% removal, adjusted against each null model',
      colour='Null model',
      shape='Site comparison',
      size='Single-site links',
      caption='Zero: importance expected under the null for this site. Triangles: approximate two-sided BH q < 0.05 within dataset x model.\nOpen circles: insufficient site-level diagnostics. Free scales; results describe recorded interactions, not environmental or sampling causes.'
    )
  
  if(fac)
    p<-p+facet_wrap(~dataset,scales='free',ncol=5)
  
  if(!fac)
    p<-p+theme(
      legend.box='vertical',
      legend.spacing.y=grid::unit(0,'pt')
    )
  
  p
}


save_fig(
  p1fun(pair_long),
  file.path(out,'42A_containment_and_turnover'),
  h=8
)

save_fig(
  p2fun(g),
  file.path(out,'42B_observed_vs_null_containment'),
  h=7
)

save_fig(
  p3fun(z),
  file.path(out,'42C_richness_adjusted_site_importance'),
  h=8
)


for(ds in s$dataset){
  
  folder<-file.path(out,ds)
  
  # Keep the same requested title for individual 42A plots
  save_fig(
    p1fun(filter(pair_long,dataset==ds),FALSE)+
      labs(title='Are sites nested?'),
    file.path(folder,'42A_containment_and_turnover'),
    8,6
  )
  
  save_fig(
    p2fun(filter(g,dataset==ds),FALSE)+
      labs(title=paste(ds,'- observed vs null')),
    file.path(folder,'42B_observed_vs_null_containment'),
    8,6
  )
  
  save_fig(
    p3fun(filter(z,dataset==ds),FALSE)+
      labs(title=paste(ds,'- richness-adjusted site importance')),
    file.path(folder,'42C_richness_adjusted_site_importance'),
    9,6
  )
  
  
  # Pairwise similarity is descriptive; no formal modules are inferred.
  pp<-filter(pair,dataset==ds)
  sites<-unique(c(pp$site1,pp$site2))
  mat<-diag(1,length(sites))
  rownames(mat)<-colnames(mat)<-sites
  
  for(k in seq_len(nrow(pp))){
    i<-match(pp$site1[k],sites)
    j<-match(pp$site2[k],sites)
    mat[i,j]<-mat[j,i]<-1-pp$sorensen[k]
  }
  
  ord<-hclust(
    as.dist(1-mat),
    method='average'
  )$order
  
  hh<-expand.grid(
    i=seq_along(ord),
    j=seq_along(ord)
  )
  
  hh$similarity<-as.vector(mat[ord,ord])
  
  save_fig(
    ggplot(hh,aes(i,j,fill=similarity))+
      geom_raster()+
      scale_fill_viridis_c(limits=c(0,1))+
      coord_equal()+
      labs(
        title=paste(ds,'- similarity among sites'),
        x='Sites ordered by similarity',
        y='Sites ordered by similarity',
        caption='Sorensen similarity; average-linkage ordering. Blocks are descriptive, not evidence of geographic barriers.'
      ),
    file.path(folder,'42D_site_similarity'),
    8,7
  )
  
  write.csv(
    data.frame(
      order=seq_along(ord),
      site=sites[ord]
    ),
    file.path(folder,'42_similarity_order.csv'),
    row.names=FALSE
  )
  
  
  tt<-filter(tr,dataset==ds)
  
  save_fig(
    ggplot(tt,aes(draw,containment,colour=factor(chain)))+
      geom_line(linewidth=.3)+
      facet_wrap(~model,ncol=1)+
      labs(
        title=paste(ds,'- null-chain traces'),
        colour='Chain',
        y='Mean containment',
        x='Saved draw'
      ),
    file.path(folder,'42E_chain_diagnostics'),
    9,6
  )
  
  
  rr<-filter(rem,dataset==ds) %>%
    group_by(model,order,fraction) %>%
    summarise(
      retention=mean(retention_mean),
      variance=mean(retention_variance),
      .groups='drop'
    )
  
  save_fig(
    ggplot(rr,aes(fraction,retention,colour=model))+
      geom_line()+
      geom_point()+
      facet_wrap(~order)+
      labs(
        title=paste(ds,'- secondary removal comparison'),
        x='Fraction of sites removed',
        y='Regional link fraction retained',
        caption='Exploratory: 100 observed removal sequences; 40 snapshots per null, 10 sequences each. Richness ties randomised.\nInterpret only where the primary null comparison is informative; fixed support implies identical exact random-removal means.'
      ),
    file.path(folder,'42F_removal_diagnostics'),
    11,5
  )
  
  write.csv(
    rr,
    file.path(folder,'42_removal_summary.csv'),
    row.names=FALSE
  )
}

message('All 42 figures saved')