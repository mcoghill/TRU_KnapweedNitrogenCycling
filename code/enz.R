library(tidyverse)
library(ggpubr)
library(extrafont)


font_style <- "Arial"
if(!font_style %in% fonts()) {
  font_import(prompt = FALSE)
}
loadfonts()


plot_colors <- scale_fill_viridis_d(begin = 0.2, end = 0.8, alpha = 2/3)$palette(4)
names(plot_colors) <- c("Control", "*C. stoebe*", "*A. millefolium*", "*V. villosa*")


df <- readxl::read_xlsx("data/combined_enzyme.xlsx") %>% 
  filter(rep != 10) %>% 
  mutate(
    plant = case_when(
      plant == "knapweed" ~ "*C. stoebe*",
      plant == "yarrow" ~ "*A. millefolium*",
      plant == "vetch" ~ "*V. villosa*",
      .default = str_to_sentence(plant)
    ),
    plant = factor(plant, levels = c("Control", "*C. stoebe*", "*A. millefolium*", "*V. villosa*")),
    cn_enz_ratio = log(cb + ag)/log(lap + nag),
    cp_enz_ratio = log(cb + ag)/log(phos),
    np_enz_ratio = log(lap + nag)/log(phos),
    vec_length = sqrt((cn_enz_ratio)^2 + (cp_enz_ratio)^2),
    vec_ang = (atan2( cn_enz_ratio, cp_enz_ratio)*(180/pi))
  )


df_long <- df %>% 
  pivot_longer(-c(plant, rep)) %>% 
  filter(name %in% c("lap", "nag")) %>% 
  mutate(value = value/1000) ## values changed from nanomole/g/hr to micromole/g/hr



enz <- c("Leucine-aminopeptidase",
                "N-acetyl-β-Glucosaminidase")
names(enz) <- c("lap", "nag")

plt <- 
ggplot(df_long , aes(x = plant, y = value, fill = plant)) +
  geom_violin(alpha = 1/3, width = 0.9, show.legend = F, linewidth = 0.25) +
  geom_boxplot(show.legend = F, width = 0.1, linewidth = 0.25, outlier.size = 0.5) +
  geom_pwc(ref.group = "Control",
           label = "p.adj.signif",
           p.adjust.method = "BH",
           p.adjust.by = "group",
           hide.ns = TRUE,
           family = font_style,
           y.position = 0.12) +
  facet_wrap(~ name , scales = "fixed",
             labeller = labeller(name = enz)) +
  theme_bw(base_size = 14, base_family = font_style) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = ggtext::element_markdown(colour = "black", angle = 22.5, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  scale_fill_manual(values = plot_colors) +
 labs(y = "Enzyme activity (µmol/g/hr) ",
      x = NULL) 
plt
ggsave(plt, file = "figures/enz_activity.png", width = 7.5, height = 5, dpi = 800, units = "in")


