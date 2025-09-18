library(tidyverse)
library(janitor)
library(ggpubr)
library(extrafont)
library(patchwork)

font_style <- "Arial"
if(!font_style %in% fonts()) {
  font_import(prompt = FALSE)
}
loadfonts()

# plot_colors <- palette.colors(4, palette = "Dark2")
plot_colors <- scale_fill_viridis_d(begin = 0.2, end = 0.8, alpha = 2/3)$palette(4)
names(plot_colors) <- c("Control", "*C. stoebe*", "*A. millefolium*", "*V. villosa*")

# Soil plot - filter to read just the first 9 samples, and plot with Wilcox
# pairwise comparisons, specifically comparing to the control
df_soil <- readxl::read_xlsx("data/combined_soil_CN.xlsx") %>% 
  clean_names() %>% 
  filter(timepoint != "Baseline", !str_detect(sample_code, pattern = "10")) %>% 
  mutate(
    treatment = case_when(
      treatment == "Knapweed" ~ "*C. stoebe*",
      treatment == "Yarrow" ~ "*A. millefolium*",
      treatment == "Vetch" ~ "*V. villosa*",
      .default = treatment
    ),
    treatment = factor(
      treatment, 
      levels = c("Control", "*C. stoebe*", "*A. millefolium*", "*V. villosa*")
    ),
    title = "Soil"
  )

soil_cn_plot <- ggplot(df_soil, aes(x = treatment, y = nitrogen, fill = treatment)) +
  geom_violin(alpha = 1/3,  show.legend = F, linewidth = 0.25) +
  geom_boxplot(width = 0.1, show.legend = F, linewidth = 0.25, outlier.size = 0.5) + 
  geom_pwc(ref.group = "Control",
           method = "wilcox_test",
           label = "p.adj.signif",
           family = font_style,
           hide.ns = TRUE,
           y.position = c(0.445, 0.455, 0.5)) +
  facet_grid(. ~ title) +
  theme_bw(base_size = 14, base_family = font_style) +
  scale_fill_manual(values = plot_colors) +
  labs(y = "Total Nitrogen content (%)") +
  theme(axis.text.x = ggtext::element_markdown(color = "black", angle = 22.5, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        panel.grid = element_blank())

soil_cn_plot
ggsave("figures/soil_cn.png", soil_cn_plot, width = 4, height = 4, units = "in", dpi = 800)

# Plant data - Same thing as soil data, except making the comparison against
# knapweed since there is no control data for aboveground material
df_plant_cn <- readxl::read_xlsx("data/combined_plant_biom.xlsx") %>% 
  clean_names() %>% 
  rename(carbon = nitrogen, nitrogen = carbon) %>% # column headings incorrect, so change those here
  filter(sample_type == "UNK",
         !str_detect(sample_id, pattern = "Q|BLANK"),
         carbon > 0) %>% 
  mutate(
    plant_type = str_extract(sample_id, pattern = "^.."),
    plant_type = case_when(
      plant_type == "HK" ~ "*C. stoebe*",
      plant_type == "HY" ~ "*A. millefolium*",
      plant_type == "HV" ~ "*V. villosa*",
      .default = plant_type
    ),
    plant_type = factor(plant_type, levels = c("*C. stoebe*", "*A. millefolium*", "*V. villosa*")),
    rep = str_extract(sample_id, patter = "A|B|C"),
    rep_n = as.numeric(str_extract(sample_id, patter = "[0-9]"))
  ) %>% 
  group_by(plant_type,rep_n) %>% 
  drop_na() %>% 
  summarise(mean_c = mean(carbon),
            mean_n = mean(nitrogen), .groups = "drop") %>% 
  mutate(title = "Plant")

plant_cn_plot <- ggplot(df_plant_cn, aes(x = plant_type, y = mean_n, fill = plant_type)) +
  geom_violin(alpha = 1/3, show.legend = F, linewidth = 0.25, width = 0.9) +
  geom_boxplot(width = 0.05, show.legend = F, linewidth = 0.25, outlier.size = 0.5) +
  #geom_point()+
  geom_pwc(ref.group = "*C. stoebe*",
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_grid(. ~ title) +
  theme_bw(base_size = 14, base_family = font_style) +
  # labs(y = "Nitrogen content (%)",
  #      x = "Treatment") +
  scale_fill_manual(values = plot_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.75, 2.505)) +
  theme(axis.text.x = ggtext::element_markdown(color = "black", angle = 22.5, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(color = "black", face = "bold"),
        strip.text = element_text(color = "black", face = "bold"),
        panel.grid = element_blank())

plant_cn_plot

ggsave("figures/plant_cn.png", plant_cn_plot, width = 4, height = 4, units = "in", dpi = 800)

# Combine plots
soil_plant_cn_plot <- soil_cn_plot + plant_cn_plot
soil_plant_cn_plot
ggsave("figures/soil_plant_cn.png", soil_plant_cn_plot, width = 7, height = 4, units = "in", dpi = 800)
