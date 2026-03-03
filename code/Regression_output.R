# Import packages ---------------------------------------------------------


library(ggpubr)
library(ggtext)
library(tidyverse)
library(glue)
library(broom)
library(sjPlot)




# Bringing together all data ----------------------------------------------


source("code/enz.R")
source("code/qpcr_analysis.R")
source("code/soil_plant_nitrogen.R")

enz_df <- df %>% 
  dplyr::rename(treatment = "plant") %>% 
  mutate(rep = as.factor(rep))

soil <- df_soil %>% 
  mutate(rep = str_extract(sample_code, pattern = ".$"))

gene_df <- list_df  %>% 
  select(-total_genecopy) %>% 
  pivot_wider(id_cols = c(treatment, Sample),
              names_from = "id",
              values_from = "log_gene")%>% 
  mutate(rep = str_extract(Sample, pattern = ".$")) %>% 
  inner_join(soil, by = c("treatment", "rep")) %>% 
  inner_join(enz_df, by = c( "treatment", "rep"))

plant_bio <- df_plant_cn %>% 
  rename(treatment = "plant_type",
         rep = "rep_n") %>% 
  mutate(rep = as.factor(rep)) %>% 
  full_join(gene_df,by = c( "treatment", "rep")) %>% 
  mutate(cn_ratio = carbon/nitrogen)


biomass <- readxl::read_xlsx("data/Moisture_Analysis_plant_biomass - Nitrogen Dynamics_20thJuly2022.xlsx",sheet = 2) %>% 
  janitor::clean_names() %>% 
  select(sample_code:root_biomass)

all_data <- inner_join(biomass, plant_bio, by = "sample_code") %>% 
  mutate(nit_mov = (mean_n*shoot_biomass)/(nitrogen*1500)*100,
         ) 


# Impact of genes and enz on plant tissues --------------------------------

md3.1 <- lm(mean_n ~  cn_ratio + nifH + `AOA-amoA` + `AOB-amoA` + narG + 
              nirK + nirS + nosZ + lap + nag, all_data )

summary(md3.1)

var_exp <- performance::r2(md3.1)$R2_adjusted %>% 
  t() %>% 
  as.data.frame() %>% 
  slice(1) %>% 
  round(.,3)


term_name <-  c("*nirK*",
                "C/N",
                "*nirS*",
                "NAG",
                "LAP",
                "*nifH*",
                "*AOA-amoA*",
                "*narG*",
                "*AOB-amoA*",
                "*nosZ*")

plt_mod <- plot_model(md3.1,
                      show.values = T,
                      value.offset = 0.4,
                      sort.est = T,
                      value.size = 5,
                      line.size = 1.25,
                      dot.size = 5,
                      title = "",
                      axis.lim = c(1, -2)) +
  theme_bw(base_size = 18) +
  annotate(geom = "text", x = 10, y = -1.5,
           label = glue::glue("adj.R sq = {var_exp$`adjusted R2`}"),
           size = 5, fontface= 3)+
    scale_x_discrete(label = term_name) +
  theme(
    panel.grid =  element_blank(),
    axis.text = element_markdown(color = "black")
  )+
  labs(y = "Regression coefficents")


plt_mod

ggsave(plt_mod, file = "figures/reg_output.jpeg", height = 6, width = 6, units = "in")



# Nitrogen update efficiency ----------------------------------------------



nupe <- all_data %>% 
  filter(treatment != "Control") %>% 
  ggplot(., aes(x = treatment, y =  log1p(nit_mov), fill = treatment)) +
  geom_violin(alpha = 1/3, show.legend = F, linewidth = 0.25, width = 0.9) +
  geom_boxplot(width = 0.05, show.legend = F, linewidth = 0.25, outlier.size = 0.5) +
  geom_pwc(ref.group = "*A. millefolium*",
           label = "p.adj.signif",
           method = "wilcox_test",
           p.adjust.method = "BH",
           family = font_style,
           hide.ns = TRUE) +
  theme_bw(base_size = 14, base_family = font_style) +
  scale_fill_manual(values = plot_colors) +
  theme(axis.text.x = ggtext::element_markdown(color = "black", angle = 22.5, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        strip.text = element_text(color = "black", face = "bold"),
        panel.grid = element_blank()) +
  labs(y = "Nitrogen update efficiency ",
       x = "Treatment") 

nupe

ggsave(nupe, file = "figures/nupe.png", width = 4, height = 5, units = "in", dpi = 800)

