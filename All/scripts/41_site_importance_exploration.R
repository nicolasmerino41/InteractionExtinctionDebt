# Run from InteractionExtinctionDebt with Rscript All/scripts/41_site_importance_exploration.R
# Exploratory structural importance; not the Ross et al. fragility formula.
# Exact expectations for uniform removal without replacement. No rewiring,
# demographic response, detection correction, or environmental interpretation.
# Compare forced loss of s + m-1 random other sites with m random sites.
# Disconnection = no observed partners AND retained occurrence, divided by the
# ORIGINAL connected species count (not a changing conditional denominator).
if(.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", ".UTF-8")
source("All/scripts/00_dataset_loaders_and_helpers_all.R")
options(warn = 1)
stopifnot(requireNamespace("patchwork", quietly = TRUE))
out_dir <- "All/outputs/41_site_importance_exploration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Probability all k supporting sites fall in a uniform m-site removal from N.
p_loss <- function(k, N, m) {
  ans <- numeric(length(k)); ok <- k <= m
  ans[ok] <- exp(lchoose(m, k[ok]) - lchoose(N, k[ok]))
  ans
}
safe_rho <- function(x, y) {
  if(length(x) < 3 || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  cor(x, y, method = "spearman")
}
incidence <- function(df, sites, ids, id_col) {
  a <- matrix(0L, length(sites), length(ids), dimnames = list(sites, ids))
  a[cbind(match(df$site, sites), match(df[[id_col]], ids))] <- 1L
  a
}
checks <- list()
check <- function(name, pass) {
  checks[[length(checks)+1L]] <<- data.frame(check = name, pass = isTRUE(pass))
  if(!isTRUE(pass)) stop("Validation failed: ", name)
}
# Exhaustive toy landscape: link sites {1,2}, focal interaction sites {1,2},
# occurrence sites {1,2,3}. Tests the joint event and every forced focal loss.
for(m in 1:4) for(s in 1:5) {
  sets <- combn(setdiff(1:5, s), m-1, simplify = FALSE)
  exact <- mean(vapply(sets, function(z) {
    removed <- c(s,z)
    all(c(1,2) %in% removed) && !all(c(1,2,3) %in% removed)
  }, logical(1)))
  pred <- p_loss(2-as.integer(s %in% 1:2),4,m-1) -
    p_loss(3-as.integer(s %in% 1:3),4,m-1)
  check(paste("toy_joint",m,s), abs(exact-pred)<1e-12)
}

results <- list(); site_tables <- list(); audits <- list()
for(dataset in all_dataset_names) {
  message("Site importance: ", dataset)
  st <- get_dataset_site_tables(dataset)
  ints <- distinct(st$empirical_site_interactions, site, consumer, resource)
  cooc <- distinct(st$cooc_triples, site, consumer, resource)
  ints[] <- lapply(ints, as.character); cooc[] <- lapply(cooc, as.character)
  occ <- st$occupancy
  occ[] <- lapply(occ, as.character)
  # Guild-specific occurrence avoids merging identical names between guilds.
  occ <- transmute(occ, site, node = paste(tolower(trophic_level), species, sep="::"))
  ii <- bind_rows(transmute(ints,site,node=paste("consumer",consumer,sep="::")),
                  transmute(ints,site,node=paste("resource",resource,sep="::"))) %>% distinct()
  check(paste(dataset,"interaction subset occupancy"), nrow(anti_join(ii,occ,by=c("site","node")))==0)
  sites <- sort(unique(c(cooc$site, ints$site, occ$site))); N <- length(sites)
  pairs <- distinct(ints,consumer,resource) %>% mutate(link=seq_len(n()))
  ints <- left_join(ints,pairs,by=c("consumer","resource"))
  M <- incidence(ints,sites,pairs$link,"link"); K <- colSums(M)
  nodes <- sort(unique(ii$node)); I <- incidence(ii,sites,nodes,"node")
  O <- incidence(filter(occ,node %in% nodes),sites,nodes,"node")
  U <- colSums(I); V <- colSums(O)
  check(paste(dataset,"occurrence contains interaction incidence"), all(I<=O))
  # O often derives from interaction webs: no independent evidence of a
  # species occurring without any partner then exists. Report, do not hide.
  audits[[dataset]] <- data.frame(dataset,N,links=ncol(M),species=length(nodes),
    species_with_occurrence_only_sites=sum(V>U), occurrence_only_records=sum(O-I))
  descriptors <- data.frame(dataset, site=sites, sites=N,
    local_links=rowSums(M), unique_links=as.vector(M %*% as.numeric(K==1)),
    inverse_support_score=as.vector(M %*% (1/K)))
  site_tables[[dataset]] <- descriptors
  counts <- sort(unique(c(1L,pmax(1L,pmin(N-1L,round(N*c(.1,.2,.4,.5,.6,.8)))))))
  for(m in counts) {
    link_base <- mean(p_loss(K,N,m))
    for(s in seq_len(N)) {
      link_loss <- mean(p_loss(K-M[s,],N-1,m-1))
      for(guild in c("consumer","resource")) {
        g <- startsWith(nodes,paste0(guild,"::"))
        gone <- p_loss(U[g]-I[s,g],N-1,m-1)
        absent <- p_loss(V[g]-O[s,g],N-1,m-1)
        disconnected <- mean(gone-absent)
        base_dis <- mean(p_loss(U[g],N,m)-p_loss(V[g],N,m))
        results[[length(results)+1L]] <- data.frame(dataset,site=sites[s],guild,
          N,removed=m,removal_fraction=m/N,link_retention=1-link_loss,
          baseline_link_retention=1-link_base,excess_link_loss=link_loss-link_base,
          disconnection=disconnected,baseline_disconnection=base_dis,
          excess_disconnection=disconnected-base_dis,
          all_partner_loss=mean(gone),baseline_all_partner_loss=mean(p_loss(U[g],N,m)),
          excess_all_partner_loss=mean(gone)-mean(p_loss(U[g],N,m)))
      }
    }
  }
}
res <- bind_rows(results); desc <- bind_rows(site_tables); audit <- bind_rows(audits)
check("probability bounds", all(res$link_retention>=-1e-12 & res$link_retention<=1+1e-12 &
  res$disconnection>=-1e-12 & res$disconnection<=res$all_partner_loss+1e-12))
avg <- res %>% group_by(dataset,guild,removed) %>% summarise(across(starts_with("excess_"),mean),.groups="drop")
check("mean forced-site contrast is zero", max(abs(as.matrix(select(avg,starts_with("excess_")))))<1e-10)
immediate <- res %>% filter(removed==1) %>% select(dataset,site,guild,
  immediate_link_retention=link_retention, immediate_disconnection=disconnection)
res <- left_join(res,immediate,by=c("dataset","site","guild")) %>%
  mutate(additional_link_loss=immediate_link_retention-link_retention,
         disconnection_change=disconnection-immediate_disconnection)
# Disconnection can fall as species themselves disappear; all_partner_loss is
# also provided to distinguish this observation constraint from robustness.
mid <- res %>% group_by(dataset,guild) %>% filter(removed==removed[which.min(abs(removal_fraction-.5))]) %>%
  ungroup() %>% left_join(desc,by=c("dataset","site")) %>%
  group_by(dataset,guild) %>% mutate(link_rank=rank(-excess_link_loss,ties.method="min"),
    disconnection_rank=rank(-excess_disconnection,ties.method="min")) %>% ungroup()
summary <- mid %>% group_by(dataset,guild) %>% summarise(
  actual_removal_fraction=first(removal_fraction),
  richness_vs_link_importance=safe_rho(local_links,excess_link_loss),
  rarity_vs_link_importance=safe_rho(inverse_support_score,excess_link_loss),
  richness_vs_disconnection=safe_rho(local_links,excess_disconnection),
  link_vs_disconnection=safe_rho(excess_link_loss,excess_disconnection),
  range_link_importance_pp=100*diff(range(excess_link_loss)),
  range_disconnection_importance_pp=100*diff(range(excess_disconnection)),.groups="drop")
write.csv2(res,file.path(out_dir,"41_exact_site_removal_outcomes.csv"),row.names=FALSE)
write.csv2(mid,file.path(out_dir,"41_site_importance_at_half_removal.csv"),row.names=FALSE)
write.csv2(summary,file.path(out_dir,"41_dataset_comparisons.csv"),row.names=FALSE)
write.csv2(audit,file.path(out_dir,"41_occurrence_information_audit.csv"),row.names=FALSE)
write.csv2(bind_rows(checks),file.path(out_dir,"41_validation_checks.csv"),row.names=FALSE)
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="bottom"))
save_plot <- function(p,name,w=12,h=7) {
  ggsave(file.path(out_dir,paste0(name,".png")),p,width=w,height=h,dpi=220,bg="white")
  ggsave(file.path(out_dir,paste0(name,".pdf")),p,width=w,height=h,bg="white")
}
p1 <- ggplot(filter(mid,guild=="consumer"),aes(local_links,100*excess_link_loss,colour=inverse_support_score))+
  geom_hline(yintercept=0,colour="grey70")+geom_point(size=2,alpha=.8)+
  facet_wrap(~dataset,scales="free",ncol=5)+scale_colour_viridis_c()+
  labs(title="Which sites matter beyond their interaction richness?",
    subtitle="About 50% of sites removed; exact removal fractions are saved in the tables",
    x="Observed interactions at focal site",y="Excess regional link loss (percentage points)",
    colour="Sum of inverse link support",
    caption="Forced focal-site loss compared with uniform removal of the same total number of sites.")
save_plot(p1,"41A_richness_and_link_importance")
unresolved <- mid %>% group_by(dataset) %>% summarise(unresolved=all(abs(disconnection)<1e-12),.groups="drop") %>% filter(unresolved) %>%
  mutate(label="Occurrence derived\nfrom interactions:\nendpoint uninformative")
p2 <- ggplot(anti_join(mid,unresolved,by="dataset"),aes(100*excess_link_loss,100*excess_disconnection,colour=guild))+
  geom_hline(yintercept=0,colour="grey70")+geom_vline(xintercept=0,colour="grey70")+
  geom_point(size=2,alpha=.8)+
  geom_text(data=unresolved,aes(x=0,y=0,label=label),inherit.aes=FALSE,size=3)+
  facet_wrap(~dataset,scales="free",ncol=5)+
  labs(title="Do the same sites maintain links and prevent species disconnection?",
    subtitle="Disconnection requires the focal species to remain recorded after removal",
    x="Excess regional link loss (percentage points)",y="Excess species disconnection (percentage points)",
    caption="Zero disconnection can reflect occurrence records derived from interaction webs; see occurrence audit.")
save_plot(p2,"41B_two_site_importance_outcomes")
# A supplementary endpoint makes structural all-partner loss visible when
# available occurrence data cannot distinguish disconnection from disappearance.
p3 <- ggplot(mid,aes(local_links,100*excess_all_partner_loss,colour=guild))+
  geom_hline(yintercept=0,colour="grey70")+geom_point(alpha=.8)+facet_wrap(~dataset,scales="free",ncol=5)+
  labs(title="Loss of every observed partner, including species disappearance",
    x="Observed interactions at focal site",y="Excess all-partner loss (percentage points)")
save_plot(p3,"41C_all_partner_loss_sensitivity")
chosen <- mid %>% group_by(dataset,guild) %>% slice_max(excess_link_loss,n=1,with_ties=FALSE) %>%
  select(dataset,guild,site)
curves <- inner_join(res,chosen,by=c("dataset","guild","site")) %>%
  select(dataset,guild,site,removal_fraction,excess_link_loss,excess_disconnection) %>%
  tidyr::pivot_longer(starts_with("excess_"),names_to="outcome",values_to="excess")
p4 <- ggplot(curves,aes(removal_fraction,100*excess,colour=guild,linetype=outcome))+
  geom_hline(yintercept=0,colour="grey70")+geom_line()+geom_point()+
  facet_wrap(~dataset,scales="free_y",ncol=5)+
  labs(title="Importance changes with the extent of site loss",
    subtitle="Each dataset uses its site with greatest link-loss importance at about 50% removal",
    x="Total fraction of sites removed",y="Excess loss (percentage points)",
    caption="Selection and evaluation use the same empirical structure; these are descriptive contrasts, not predictive validation.")
save_plot(p4,"41D_importance_across_removal")
writeLines(c("Exploratory site importance: exact uniform site removal",
  "Forced-site scenario removes the focal site plus m-1 random other sites; baseline removes m random sites.",
  "Positive excess loss means greater damage than the size-m random baseline. Averaging over focal sites returns zero.",
  "Expected link retention uses each link's support. Species all-partner loss uses the union of its interaction sites.",
  "Disconnection while present subtracts probability of losing all occurrence sites from probability of losing all interaction sites.",
  "Species denominators are the original connected species counts within each guild.",
  "Occurrence may be derived from interactions. Zero disconnection is not evidence that interactions cannot be lost before species.",
  "Immediate damage and later total damage are both saved. Disconnection change may be negative as species disappear.",
  "Inverse-support scores are descriptive rarity-weighted richness, not Ross et al. network fragility.",
  "No environmental metadata, abundance, sampling-effort correction or causal habitat-loss inference is assumed.",
  "Correlations are within-dataset descriptive comparisons; constant outcomes return NA, not zero.",
  "Plots have free panel scales. Exact removal fractions vary with dataset size. No pooled significance tests are performed."),
  file.path(out_dir,"41_README.txt"))
print(audit); print(summary,width=Inf)
message("Finished: ",out_dir)
