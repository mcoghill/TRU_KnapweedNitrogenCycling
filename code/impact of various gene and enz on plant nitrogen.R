library(ggpubr)
library(ggtext)
library(tidyverse)
library(glue)
library(janitor)
library(readxl)
library(extrafont)

font_style <- "Arial"
if(!font_style %in% fonts()) {
  font_import(prompt = FALSE)
}
loadfonts()

################################################################################
################################################################################
##### Run all prior codes from knapweed analysis ###############################
################################################################################
################################################################################


df <- readxl::read_xlsx("data/combined_enzyme.xlsx") %>% 
  filter(rep != 10) %>% 
  mutate(
    plant = case_when(
      plant == "knapweed" ~ "*C. stoebe*",
      plant == "yarrow" ~ "*A. millefolium*",
      plant == "vetch" ~ "*V. villosa*",
      .default = str_to_sentence(plant)
    ),
    plant = factor(plant, levels = c("Control", "*C. stoebe*", "*A. millefolium*", "*V. villosa*"))
  )

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

list_df <- excel_sheets("data/qPCR_data.xlsx") %>% 
  set_names() %>% 
  map(readxl::read_xlsx, path = "data/qPCR_data.xlsx") %>% 
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
  mutate(shoot_root_ratio = shoot_biomass/root_biomass,
         cn_enz_ratio = log(cb + ag)/log(lap + nag),
         cp_enz_ratio = log(cb + ag)/log(phos),
         np_enz_ratio = log(lap + nag)/log(phos),
         plant_cn_ratio = mean_c/mean_n,
         root_shoot = root_biomass/shoot_biomass,
         shoot_n_ratio = (mean_n)/shoot_biomass,
         tot_nit = (mean_n*10)*(shoot_biomass),
         nit_mov = mean_n/nitrogen)



md3 <- lm(mean_n ~ cn_enz_ratio + cn_ratio + nifH + `AOA-amoA` + `AOB-amoA` + narG + nirK + nirS + nosZ + lap + nag + cb + phos + ag, all_data )
summary.aov(md3)


md3.1 <- lm(mean_n ~  cn_ratio + nifH + `AOA-amoA` + `AOB-amoA` + narG + nirK + nirS + nosZ + lap + nag+ cb + phos + ag, all_data )
mod <- anova(md3.1)
summ_md3.1 <- broom::tidy(mod)

AIC(md3, md3.1)

shapiro.test(resid(md3.1))
performance::check_heteroscedasticity(md3.1)
caret::varImp(md3.1)

##Selected model - normality and heterosckasticity okay
mod


#figures for model md3.1
plt_n1 <- ggplot(all_data, aes(x = cn_ratio, y = mean_n)) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(pch = 21, fill= "#414487AA", size = 5, color = "black") +
  theme_bw(base_size = 18, base_family = font_style) +
  # coord_cartesian(xlim = c(11.4, 12.1))+
  annotate(geom = "text",
           x = -Inf,
           y = -Inf,
           label = glue("{'  '}F = {round(summ_md3.1$statistic[1],2)}\n{'  '}p = {round(summ_md3.1$p.value[1], 4)}"),
           hjust = 0, vjust = -0.3,
           size = 6,
           fontface = "italic") +
  annotate(geom = "text",
           x = -Inf,
           y = Inf,
           label = "A",
           hjust = 4, vjust = 1,
           size = 6,
           fontface = "bold") +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank()) +
  labs(x = "Soil C/N ratio",
       y = "Plant nitrogen content (%)") +
  coord_cartesian(expand = TRUE, xlim = c(11.377, NA), clip = "off")

plt_n1
# ggsave(plt_n1, file = "figures/plant_nitrogen_vs_cn_ratio.png", width = , height = 6, dpi = 800, units = "in")


plt_n2 <- ggplot(all_data, aes(x = `AOA-amoA`, y = mean_n )) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(pch = 21, fill= "#414487AA", size = 5, color = "black") +
  theme_bw(base_size = 18, base_family = font_style) +
  # coord_cartesian(xlim = c(16.6, 17.75))+
  annotate(geom = "text",
           x = Inf,
           y = -Inf,
           label = glue("F = {round(summ_md3.1$statistic[3],2)}  \np = {round(summ_md3.1$p.value[3], 4)}  "),
           hjust = 1, vjust = -0.3,
           size = 6,
           fontface = "italic") +
  annotate(geom = "text",
           x = -Inf,
           y = Inf,
           label = "B",
           hjust = 4, vjust = 1,
           size = 6,
           fontface = "bold") +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_text(face = "italic"),
        panel.grid = element_blank()) +
  labs(x = "AOA-amoA",
       y = "Plant nitrogen content (%)") +
  coord_cartesian(expand = TRUE, clip = "off")

plt_n2
# ggsave(plt_n2, file = "figures/plant_nitrogen_vs_aoa.png", width = , height = 6, dpi = 800, units = "in")


plt_n3 <- ggplot(all_data, aes(x = `AOB-amoA`, y = mean_n )) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(pch = 21, fill= "#414487AA", size = 5, color = "black") +
  theme_bw(base_size = 18, base_family = font_style) +
  # coord_cartesian(xlim = c(12.75, 15.25))+
  annotate(geom = "text",
           x = Inf,
           y = -Inf,
           label = glue("F = {round(summ_md3.1$statistic[4],2)}  \np = {round(summ_md3.1$p.value[4], 3)}  "),
           hjust = 1, vjust = -0.3,
           size = 6,
           fontface = "italic") +
  annotate(geom = "text",
           x = -Inf,
           y = Inf,
           label = "C",
           hjust = 4, vjust = 1,
           size = 6,
           fontface = "bold") +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_text(face = "italic"),
        panel.grid = element_blank()) +
  labs(x = "AOB-amoA",
       y = "Plant nitrogen content (%)") +
  coord_cartesian(expand = TRUE, clip = "off")

plt_n3
# ggsave(plt_n3, file = "figures/plant_nitrogen_vs_aob.png", width = , height = 6, dpi = 800, units = "in")


plt_n4 <- ggplot(all_data, aes(x = nirK, y = mean_n )) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(pch = 21, fill= "#414487AA", size = 5, color = "black") +
  theme_bw(base_size = 18, base_family = font_style) +
  # coord_cartesian(xlim = c(18.25, 19.55))+
  annotate(geom = "text",
           x = -Inf,
           y = -Inf,
           label = glue("{'  '}F = {round(summ_md3.1$statistic[6],2)}\n{'  '}p = {round(summ_md3.1$p.value[6], 3)}"),
           hjust = 0, vjust = -0.3,
           size = 6,
           fontface = "italic") +
  annotate(geom = "text",
           x = -Inf,
           y = Inf,
           label = "D",
           hjust = 4, vjust = 1,
           size = 6,
           fontface = "bold") +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_text(face = "italic"),
        panel.grid = element_blank()) +
  labs(x = "nirK",
       y = "Plant nitrogen content (%)") +
  coord_cartesian(expand = TRUE, clip = "off")

plt_n4
# ggsave(plt_n4, file = "figures/plant_nitrogen_vs_nirK.png", width = , height = 6, dpi = 800, units = "in")


plt_n5 <- ggplot(all_data, aes(x = nosZ, y = mean_n )) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(pch = 21, fill= "#414487AA", size = 5, color = "black") +
  theme_bw(base_size = 18, base_family = font_style) +
  # coord_cartesian(xlim = c(18, 19.25))+
  annotate(geom = "text",
           x = Inf,
           y = -Inf,
           label = glue("F = {round(summ_md3.1$statistic[8],2)}  \np = {round(summ_md3.1$p.value[8], 3)}  "),
           hjust = 1, vjust = -0.3,
           size = 6,
           fontface = "italic") +
  annotate(geom = "text",
           x = -Inf,
           y = Inf,
           label = "E",
           hjust = 4.5, vjust = 1,
           size = 6,
           fontface = "bold") +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_text(face = "italic"),
        panel.grid = element_blank()) +
  labs(x = "nosZ",
       y = "Plant nitrogen content (%)") +
  coord_cartesian(expand = TRUE, clip = "off")

plt_n5
# ggsave(plt_n5, file = "figures/plant_nitrogen_vs_nosZ.png", width = , height = 6, dpi = 800, units = "in")


plt_n6 <- ggplot(all_data, aes(x = lap, y = mean_n )) +
  geom_smooth(method = "lm", color = "black") +
  geom_point(pch = 21, fill= "#414487AA", size = 5, color = "black") +
  theme_bw(base_size = 18, base_family = font_style) +
  # coord_cartesian(xlim = c(90, 142))+
  annotate(geom = "text",
           x = Inf,
           y = -Inf,
           label = glue("F = {round(summ_md3.1$statistic[9],2)}  \np = {round(summ_md3.1$p.value[9], 3)}  "),
           hjust = 1, vjust = -0.3,
           size = 6,
           fontface = "italic") +
  annotate(geom = "text",
           x = -Inf,
           y = Inf,
           label = "F",
           hjust = 4.6, vjust = 1,
           size = 6,
           fontface = "bold") +
  theme(axis.text = element_text(color = "black"),
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank()) +
  labs(x = "Leucine-aminopeptidase",
       y = "Plant nitrogen content (%)") +
  coord_cartesian(expand = TRUE, xlim = c(91.193, NA), clip = "off")

plt_n6
# ggsave(plt_n6, file = "figures/plant_nitrogen_vs_lap.png", width = 12, height = 8, dpi = 800, units = "in")


##combined plot

plant_nitroplot <- gridExtra::arrangeGrob(plt_n1, plt_n2, plt_n3, plt_n4, plt_n5, plt_n6, ncol = 3, nrow = 2)
# plant_nitroplot <- ggarrange(plt_n1, plt_n2, plt_n3, plt_n4, plt_n5, plt_n6, ncol = 3, nrow = 2,
#                              labels = "AUTO")

plant_nitroplot
ggsave(plant_nitroplot, file = "figures/plant_nitroplot.jpg",
       width = 16, height = 10, dpi = 800, units = "in")


