## ------------------------------------------------------------
## Script: All/scripts/36_Figure1_hidden_local_support_concept.R
## Figure 1: Regional persistence can hide local functional loss
## Outputs PNG only: complete figure and panels separately
## ------------------------------------------------------------

packages <- c("ggplot2", "dplyr", "tibble", "patchwork", "ggforce")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/36_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(p, name, width = 8, height = 6){
  ggsave(file.path(out_dir, name), p, width = width, height = height, dpi = 320, bg = "white")
}

cols <- list(
  site = "#E8E8E8",
  cooc = "#9E9E9E",
  interact = "#0072B2",
  removed = "#D9D9D9",
  lost = "#D55E00",
  present = "#009E73",
  text = "#222222"
)

theme_fig <- function(base_size = 12){
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, colour = cols$text, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0, colour = cols$text, size = base_size - 1),
      plot.margin = margin(8, 8, 8, 8)
    )
}

## Synthetic site layout
sites <- tibble(
  site = paste0("Site ", 1:12),
  x = c(1,2.1,3.2,4.2,1.4,2.7,3.8,4.7,1.2,2.4,3.5,4.5),
  y = c(3.2,3.8,3.1,3.7,2.1,2.6,2.0,2.7,1.0,1.4,0.9,1.5),
  removed = site %in% c("Site 2", "Site 4", "Site 7", "Site 9", "Site 11")
)

pair_sites <- tibble(
  pair = rep(c("Pair A", "Pair B"), each = 6),
  site = c("Site 1", "Site 2", "Site 4", "Site 6", "Site 7", "Site 11",
           "Site 1", "Site 3", "Site 4", "Site 8", "Site 10", "Site 12"),
  interaction = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
) %>% left_join(sites, by = "site")

## Panel A: local landscape before removal
panel_a <- ggplot() +
  geom_point(data = sites, aes(x, y), size = 9, colour = cols$site) +
  geom_point(data = pair_sites %>% filter(pair == "Pair A"), aes(x, y), size = 8, shape = 21,
             fill = "white", colour = cols$cooc, stroke = 1.4) +
  geom_point(data = pair_sites %>% filter(pair == "Pair B"), aes(x, y), size = 5.7, shape = 21,
             fill = "white", colour = cols$cooc, stroke = 1.4) +
  geom_point(data = pair_sites %>% filter(interaction), aes(x, y), size = 3.4,
             colour = cols$interact) +
  annotate("text", x = 0.65, y = 4.15, label = "Pair A: interaction observed in one site", hjust = 0, size = 3.6) +
  annotate("text", x = 0.65, y = 3.85, label = "Pair B: interaction observed repeatedly", hjust = 0, size = 3.6) +
  annotate("point", x = 0.84, y = 3.50, shape = 21, size = 4, fill = "white", colour = cols$cooc, stroke = 1.1) +
  annotate("text", x = 1.03, y = 3.50, label = "Sites where species are recorded together", hjust = 0, size = 3.3) +
  annotate("point", x = 0.84, y = 3.24, size = 2.8, colour = cols$interact) +
  annotate("text", x = 1.03, y = 3.24, label = "Sites where interaction is observed", hjust = 0, size = 3.3) +
  coord_equal(xlim = c(0.45, 5.2), ylim = c(0.55, 4.35)) +
  labs(title = "A. Local sites") + theme_fig()

## Panel B: binary regional network
nodes <- tibble(
  node = c("Plant A", "Pollinator x", "Plant B", "Pollinator y"),
  type = c("Plant", "Pollinator", "Plant", "Pollinator"),
  x = c(1, 3, 1, 3), y = c(2.7, 2.7, 1.35, 1.35)
)
edges <- tibble(x = c(1,1), y = c(2.7,1.35), xend = c(3,3), yend = c(2.7,1.35))
panel_b <- ggplot() +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1.3, colour = cols$interact) +
  geom_point(data = nodes, aes(x, y, fill = type), shape = 21, size = 10, colour = "grey25", stroke = 0.5) +
  geom_text(data = nodes, aes(x, y, label = node), size = 3.2) +
  scale_fill_manual(values = c("Plant" = "#F0E442", "Pollinator" = "#56B4E9"), guide = "none") +
  annotate("text", x = 2, y = 3.35, label = "Both pairs are shown as one regional interaction link", size = 4.1, fontface = "bold") +
  annotate("text", x = 2, y = 0.65, label = "Binary aggregation treats thin and repeated support as equivalent", size = 3.4) +
  coord_equal(xlim = c(0.3,3.7), ylim = c(0.45,3.6)) + labs(title = "B. Regional binary network") + theme_fig()

## Panel C: after site removal
post <- pair_sites %>% mutate(
  status = case_when(
    removed ~ "Removed site",
    interaction ~ "Interaction still observed",
    TRUE ~ "Recorded together only"
  )
)
panel_c <- ggplot() +
  geom_point(data = sites, aes(x, y), size = 9, colour = cols$site) +
  geom_point(data = sites %>% filter(removed), aes(x, y), size = 9, colour = cols$removed) +
  geom_text(data = sites %>% filter(removed), aes(x, y, label = "×"), size = 6, colour = "white", fontface = "bold") +
  geom_point(data = post %>% filter(pair == "Pair A" & !removed), aes(x, y), size = 8, shape = 21,
             fill = "white", colour = cols$cooc, stroke = 1.2) +
  geom_point(data = post %>% filter(pair == "Pair B" & !removed), aes(x, y), size = 5.7, shape = 21,
             fill = "white", colour = cols$cooc, stroke = 1.2) +
  geom_point(data = post %>% filter(interaction & !removed), aes(x, y), size = 3.4, colour = cols$present) +
  annotate("text", x = 0.65, y = 4.18, label = "Pair A: regional interaction disappears", hjust = 0, size = 3.7, colour = cols$lost) +
  annotate("text", x = 0.65, y = 3.88, label = "Pair B: regional interaction remains, but with less support left", hjust = 0, size = 3.7, colour = cols$present) +
  annotate("text", x = 0.65, y = 0.55, label = "Regional interaction presence does not reveal how broadly the interaction is sustained across sites.", hjust = 0, size = 3.5) +
  coord_equal(xlim = c(0.45, 5.2), ylim = c(0.35, 4.35)) +
  labs(title = "C. After site removal") + theme_fig()

## Save panels separately and complete figure
save_png(panel_a, "Figure1A_local_landscape.png", 7, 5)
save_png(panel_b, "Figure1B_regional_binary_network.png", 6, 5)
save_png(panel_c, "Figure1C_after_site_removal.png", 7, 5)

fig1 <- (panel_a | panel_b | panel_c) +
  plot_annotation(
    title = "Regional persistence can hide local functional loss",
    subtitle = "Binary regional links can remain visible even as the local sites supporting them are eroded",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 11))
  )
save_png(fig1, "Figure1_hidden_local_support_concept.png", 16, 5.8)

message("Saved Figure 1 PNG outputs in: ", out_dir)
