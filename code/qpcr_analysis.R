library(tidyverse)
library(janitor)
library(vegan)
library(readxl)
library(ggpubr)
library(ggh4x)
library(extrafont)

font_style <- "Arial"
if(!font_style %in% fonts()) {
  font_import(prompt = FALSE)
}
loadfonts()

path <- "data/qPCR_data.xlsx"

# plot_colors <- palette.colors(4, palette = "Dark2")
plot_colors <- scale_fill_viridis_d(begin = 0.2, end = 0.8, alpha = 2/3)$palette(4)
names(plot_colors) <- c("Control", "*C. stoebe*", "*A. millefolium*", "*V. villosa*")

list_df <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(readxl::read_xlsx, path = path) %>% 
  list_rbind(names_to = "id") %>% 
  select(Sample,id, total_genecopy) %>% 
  mutate(
    id = case_when(
      id == "AOB" ~ "AOB-amoA",
      id == "AOA" ~ "AOA-amoA",
      .default = as.character(id)
    )
  )  %>% 
  mutate(
    treatment = str_extract(Sample, pattern = "[A-z]"),
    treatment = case_when(
      treatment == "C" ~ "Control",
      treatment == "K" ~ "*C. stoebe*",
      treatment == "Y" ~ "*A. millefolium*",
      treatment == "V" ~ "*V. villosa*",
      .default = NA
    ),
    treatment = factor(
      treatment, 
      levels = c("Control", "*C. stoebe*", "*A. millefolium*", "*V. villosa*")
    ),
    log_gene = log(total_genecopy),
    id = factor(
      id, 
      levels = c("nifH", "AOB-amoA", "AOA-amoA", "narG", "nirK", "nirS", "nosZ")
    )
  ) 


gene_plot <- ggplot(list_df, aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = T, linewidth = 0.25) +
  geom_boxplot(width = 0.1, show.legend = F, linewidth = 0.25, outlier.size = 0.5) +
  geom_pwc(ref.group = "Control",
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style,
           hide.ns = TRUE,
           bracket.nudge.y = 0.15,
           step.increase = 0.15
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  facetted_pos_scales(
    y = list(
      id == "nifH" ~ scale_y_continuous(limits = c(NA,20)),
      id == "AOB-amoA" ~ scale_y_continuous(limits = c(NA,16.1)),
      id == "AOA-amoA" ~ scale_y_continuous(limits = c(NA,18)),
      id == "narG"  ~ scale_y_continuous(limits = c(NA,19.5)),
      id == "nirK"  ~ scale_y_continuous(limits = c(NA,20.1)),
      id == "nirS"  ~ scale_y_continuous(limits = c(NA,20.5)),
      id == "nosZ"  ~ scale_y_continuous(limits = c(NA,20))
    )) +
  theme(legend.direction = "vertical",
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.125),
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank()) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  labs(y = "Log (gene)")

gene_plot

ggsave(gene_plot, file = "figures/gene.png", width = 8, height = 10, dpi = 800)

# Save each plot separately
nifH_plot <- list_df %>% 
  filter(id == "nifH") %>% 
  ggplot(aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F) +
  geom_boxplot(width = 0.1, show.legend = T) +
  geom_pwc(ref.group = 1,
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  scale_y_continuous(limits = c(NA, 20)) + 
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent")) +
  labs(y = "Log (gene)")

nifH_plot
ggsave(nifH_plot, file = "figures/nifH.png", width = 6, height = 4, dpi = 800)

# AOA
aoa_plot <- list_df %>% 
  filter(id == "AOA") %>% 
  ggplot(aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F) +
  geom_boxplot(width = 0.1, show.legend = T) +
  geom_pwc(ref.group = 1,
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  scale_y_continuous(limits = c(NA, 18)) + 
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent")) +
  labs(y = "Log (gene)")

aoa_plot
ggsave(aoa_plot, file = "figures/AOA.png", width = 6, height = 4, dpi = 800)

# AOB
aob_plot <- list_df %>% 
  filter(id == "AOB") %>% 
  ggplot(aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F) +
  geom_boxplot(width = 0.1, show.legend = T) +
  geom_pwc(ref.group = 1,
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  scale_y_continuous(limits = c(NA, 16.1)) + 
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent")) +
  labs(y = "Log (gene)")

aob_plot
ggsave(aob_plot, file = "figures/AOB.png", width = 6, height = 4, dpi = 800)

# narG
narG_plot <- list_df %>% 
  filter(id == "narG") %>% 
  ggplot(aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F) +
  geom_boxplot(width = 0.1, show.legend = T) +
  geom_pwc(ref.group = 1,
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  scale_y_continuous(limits = c(NA, 19.5)) + 
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent")) +
  labs(y = "Log (gene)")

narG_plot
ggsave(narG_plot, file = "figures/narG.png", width = 6, height = 4, dpi = 800)

# nirK
nirK_plot <- list_df %>% 
  filter(id == "nirK") %>% 
  ggplot(aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F) +
  geom_boxplot(width = 0.1, show.legend = T) +
  geom_pwc(ref.group = 1,
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  scale_y_continuous(limits = c(NA, 20.1)) + 
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent")) +
  labs(y = "Log (gene)")

nirK_plot
ggsave(nirK_plot, file = "figures/nirK.png", width = 6, height = 4, dpi = 800)

# nirS
nirS_plot <- list_df %>% 
  filter(id == "nirS") %>% 
  ggplot(aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F) +
  geom_boxplot(width = 0.1, show.legend = T) +
  geom_pwc(ref.group = 1,
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  scale_y_continuous(limits = c(NA, 20.5)) + 
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent")) +
  labs(y = "Log (gene)")

nirS_plot
ggsave(nirS_plot, file = "figures/nirS.png", width = 6, height = 4, dpi = 800)

# nosZ
nosZ_plot <- list_df %>% 
  filter(id == "nosZ") %>% 
  ggplot(aes(x = treatment, y = log_gene, fill = treatment)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F) +
  geom_boxplot(width = 0.1, show.legend = T) +
  geom_pwc(ref.group = 1,
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style
  ) +
  facet_wrap(~ id, scales = "free", nrow = 4) +
  theme_bw(base_size = 18, base_family = font_style) +
  scale_fill_manual(values = plot_colors, name = "Treatment") +
  scale_y_continuous(limits = c(NA, 20)) + 
  theme(legend.direction = "vertical",
        legend.position = "right",
        legend.title = element_text(hjust = 0.5, color = "black", face = "bold"),
        legend.text =  ggtext::element_markdown(color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =  element_text(color = "black", face = "bold.italic"),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent")) +
  labs(y = "Log (gene)")

nosZ_plot
ggsave(nosZ_plot, file = "figures/nosZ.png", width = 6, height = 4, dpi = 800)
