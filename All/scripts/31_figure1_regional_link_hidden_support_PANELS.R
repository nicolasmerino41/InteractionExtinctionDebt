## ------------------------------------------------------------
## Script: All/scripts/31_figure1_regional_link_hidden_support.R
## Figure 1 only: conceptual illustration.
## No empirical data used.
##
## This version saves:
##   - Panel A separately
##   - Panel B separately
##   - Panel C separately
##   - merged Figure 1 as before
## ------------------------------------------------------------

packages <- c("ggplot2", "dplyr", "tidyr", "patchwork")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/31_hidden_support_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(p, name, width = 8, height = 6){
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width, height = height)
}

## ---------------------------
## Colours
## ---------------------------

col_site <- "#D9D9D9"
col_cooc <- "#6B7280"
col_interaction <- "#0072B2"
col_removed <- "#EFEFEF"
col_link <- "#E69F00"
col_consumer <- "#BFD7EA"
col_resource <- "#C7E9C0"

## ---------------------------
## Toy sites
## ---------------------------

sites <- data.frame(
  site = paste0("Site ", 1:12),
  x = c(1,2,3,4,1,2,3,4,1,2,3,4),
  y = c(3,3,3,3,2,2,2,2,1,1,1,1),
  removed = paste0("Site ", 1:12) %in% c("Site 2", "Site 4", "Site 7", "Site 9", "Site 11"),
  stringsAsFactors = FALSE
)

pair_sites <- expand.grid(
  pair = c("Pair A: thin local support", "Pair B: repeated local support"),
  site = sites$site,
  stringsAsFactors = FALSE
) %>%
  left_join(sites, by = "site") %>%
  mutate(
    recorded_together = case_when(
      pair == "Pair A: thin local support" ~ site %in% paste0("Site ", c(1,2,3,4,5,6)),
      pair == "Pair B: repeated local support" ~ site %in% paste0("Site ", c(4,5,6,7,8,10)),
      TRUE ~ FALSE
    ),
    interaction_observed = case_when(
      pair == "Pair A: thin local support" ~ site == "Site 2",
      pair == "Pair B: repeated local support" ~ site %in% paste0("Site ", c(4,5,6,8,10)),
      TRUE ~ FALSE
    )
  )

## ---------------------------
## Panel A: local sites
## ---------------------------

panel_a <- ggplot() +
  geom_point(
    data = pair_sites,
    aes(x = x, y = y),
    size = 7,
    colour = col_site,
    fill = col_site,
    shape = 21
  ) +
  geom_point(
    data = subset(pair_sites, recorded_together),
    aes(x = x, y = y),
    size = 7,
    colour = col_cooc,
    fill = "white",
    shape = 21,
    stroke = 1.2
  ) +
  geom_point(
    data = subset(pair_sites, interaction_observed),
    aes(x = x, y = y),
    size = 3.4,
    colour = col_interaction,
    fill = col_interaction,
    shape = 21
  ) +
  facet_wrap(~ pair, nrow = 1) +
  coord_equal() +
  theme_void(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  ggtitle(
    "A. Local sites",
    subtitle = "Open circles: sites where species are recorded together. Blue dots: sites where interaction is observed."
  )

## ---------------------------
## Panel B: binary regional network
## ---------------------------

network <- data.frame(
  x = c(1, 2, 1, 2),
  y = c(2, 2, 1, 1),
  node = c("Consumer A", "Resource x", "Consumer B", "Resource y"),
  type = c("Consumer", "Resource", "Consumer", "Resource"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  x = c(1, 1),
  y = c(2, 1),
  xend = c(2, 2),
  yend = c(2, 1)
)

panel_b <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 1.2,
    colour = col_link
  ) +
  geom_point(
    data = network,
    aes(x = x, y = y, fill = type),
    size = 7,
    shape = 21,
    colour = "grey25"
  ) +
  geom_text(
    data = network,
    aes(x = x, y = y - 0.18, label = node),
    size = 3.2
  ) +
  scale_fill_manual(values = c("Consumer" = col_consumer, "Resource" = col_resource)) +
  coord_equal(xlim = c(0.65, 2.35), ylim = c(0.65, 2.35)) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  ggtitle(
    "B. Regional network",
    subtitle = "Both pairs appear as one regional interaction link."
  )

## ---------------------------
## Panel C: after deterministic site loss
## ---------------------------

after <- pair_sites %>%
  mutate(
    still_retained = !removed,
    interaction_after = interaction_observed & still_retained,
    support_status = case_when(
      removed ~ "Removed site",
      interaction_after ~ "Sites still supporting interaction",
      recorded_together & still_retained ~ "Species still recorded together",
      TRUE ~ "Sampled site"
    ),
    support_status = factor(
      support_status,
      levels = c(
        "Sites still supporting interaction",
        "Species still recorded together",
        "Removed site",
        "Sampled site"
      )
    )
  )

panel_c <- ggplot() +
  geom_point(
    data = after,
    aes(x = x, y = y, fill = support_status),
    size = 7,
    shape = 21,
    colour = "grey50"
  ) +
  scale_fill_manual(
    values = c(
      "Removed site" = col_removed,
      "Sites still supporting interaction" = col_interaction,
      "Species still recorded together" = "white",
      "Sampled site" = col_site
    ),
    drop = FALSE
  ) +
  facet_wrap(~ pair, nrow = 1) +
  coord_equal() +
  theme_void(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 4))) +
  ggtitle(
    "C. After site loss",
    subtitle = "Pair A loses its only supporting site; Pair B still has support left."
  )

## ---------------------------
## Save individual panels
## ---------------------------

save_both(panel_a, "31_Figure1A_local_sites", width = 8.5, height = 3.4)
save_both(panel_b, "31_Figure1B_regional_network", width = 5.2, height = 4.2)
save_both(panel_c, "31_Figure1C_after_site_loss", width = 8.5, height = 3.8)

## ---------------------------
## Save merged figure too
## ---------------------------

fig <- (panel_a / panel_b / panel_c) +
  patchwork::plot_layout(heights = c(1, 0.95, 1.1)) +
  patchwork::plot_annotation(
    title = "A regional link does not reveal its local support",
    subtitle = "Conceptual illustration: binary aggregation treats interactions as present or absent, even when local support differs strongly."
  )

save_both(fig, "31_Figure1_regional_link_hidden_support", width = 11, height = 10)

message("Saved Figure 1 separate panels and merged output in: ", out_dir)
