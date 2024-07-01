#Libraries----
library(tidyverse)
library(readxl)
library(iNEXT)

library(sbm)


main_theme = theme_minimal()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        legend.text = element_text(colour = "black", size=18),
        axis.title= element_text(size=28),
        strip.text = element_text(colour = "black", size=15, face ="italic"))

#Et sinon...
#Metadonnees grilles----
read.table("data/data_clean/Metadata_grid_CAM.txt", header = T) -> metadata_grid
read.table("data/data_clean/Metadata_quadra_CAM.txt", header = T) -> metadata_quad
#remodeler le jeu de données
metadata_quad %>%  
  filter(!str_detect(Grid_code, "22_CAM_15")) %>% 
  dplyr::select(Grid_code, Host_code, Collection_date,
                Vegetation, Biomass) %>%
  separate(., Collection_date, into = c("Year", "Month", "Day"), sep = "-") %>% 
  dplyr::select(-Host_code) %>% 
  group_by(Grid_code) %>% 
  mutate(Vegetation = 2*sum(Vegetation)/1800, 
         Biomass = sum(Biomass)) %>% 
  ungroup() %>% 
  arrange(Grid_code )%>%
  distinct() -> df


#Plants----
##Lire et formater les données ----
read.table("data/data_clean/OTU_plant_CAM.txt", header = T)%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  summarise_if(is.numeric, sum)%>%
  pivot_longer(-c(Grid_code), names_to = "Plant", values_to = "occur")%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Plant")-> OUT_plant_grid_binary_df

rbind(n = rep(9,42),OUT_plant_grid_binary_df) -> OUT_plant_grid_binary_df

read.table("data/data_clean/abund_plant_grid.txt", header = T)%>%
  filter(!str_detect(Grid_code, "22_CAM_15"))%>%
  pivot_longer(-c(Grid_code), names_to = "Plant", values_to = "occur")%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Plant")-> OUT_plant_grid_abund
rbind(n = rep(1800,length(colnames(OUT_plant_grid_abund))),OUT_plant_grid_abund) -> plant_inext_formatted_data


##Estimations de diversité ----
# y <- iNEXT(plant_inext_formatted_data, q=c(0,1), datatype="incidence_freq", se=T)
# save(plant_inext_formatted_data, y, file = "outputs/plant_diversity_estimates_with_inext.Rdata")

load(file = "outputs/plant_diversity_estimates_with_inext.Rdata")

y$iNextEst$size_based %>% 
  filter(Order.q==0) %>%
  mutate(qD_obs = if_else(Method == "Observed", qD, NA),
         qD_rar = if_else(Method == "Extrapolation", NA, qD),
         qD_ext = if_else(Method == "Rarefaction", NA, qD)) %>% 
  ggplot(., aes(x = t, y = qD, group = Assemblage))+
  geom_point(aes(x = t, y = qD_obs)) +
  geom_line(aes(x = t, y = qD_rar, group = Assemblage), linewidth = 0.5) +
  geom_line(aes(x = t, y = qD_ext, group = Assemblage), linetype = 2, linewidth = 0.5)+
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), alpha = 0.05, linetype = 0)+
  ylab("Plant richness")+
  xlab("Sample size")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )

y$iNextEst$size_based %>% 
  filter(Method == "Observed" & Order.q == 0) %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         plant_richness = qD, 
         pl_rich_LCL = qD.LCL, 
         pl_rich_UCL = qD.UCL) -> plant_richness


y$iNextEst$size_based %>% 
  filter(Method == "Observed" & Order.q == 1) %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         plant_H = qD, 
         pl_H_LCL = qD.LCL, 
         pl_H_UCL = qD.UCL) -> plant_H


#Relation richesse-richesse----
load(file="Results/SBM_abund_grid.Rdata")
length(as.factor(SBM_abund_grid$memberships$Grids))
SBM_abund_grid
df%>%
  left_join( plant_richness, by = "Grid_code")%>%
  left_join(plant_H, by = "Grid_code")%>%
  
  mutate(log_plant_richness = log(plant_richness+1),
         log_plant_H = log(plant_H+1),
         log_biomass = log(Biomass+1),
         sbm_plant_grid = as.factor(SBM_abund_grid$memberships$Grids),
         log_plant_equit = log_plant_H-log_plant_richness, 
         plant_equit = plant_H/plant_richness) -> df
df$plant_H
df$plant_richness
df$Vegetation
#write.table(df, file = "outputs/dataframe_diversities.txt")

library(GGally)

df %>% 
  dplyr::select( plant_richness, plant_H, plant_equit, Biomass) %>% 
  ggpairs()

df %>% 
  dplyr::select( log_plant_richness, log_plant_H, log_plant_equit, log_biomass, Vegetation) %>% 
  ggpairs()


read.table("data/data_clean/Metadata_grid_grp_env.txt", header = T)%>%
  dplyr::select(clust_landscp, clust_land_use, sbm_soil, Grid_code)%>%
  mutate(across(where(is.numeric), as.factor))-> metadata_grid_grp_env

df%>%
  left_join(metadata_grid_grp_env, by ="Grid_code") -> df

#GLOBAL----

library(MASS)
library(MuMIn)
library(car)
library("jtools")
library(interactions)
library(cowplot)

library(emmeans)
library(viridis)
#On transforme virus_richness en log et on tente le glm gaussien 

##### Plante diversity
full.model1 <- glm(log_plant_H ~ (Vegetation +log_biomass )+
                     (clust_landscp + clust_land_use + sbm_soil), data = df,
                   na.action = "na.fail")
#####
sexpr = expression(!(
    (clust_land_use  &&  clust_landscp)|
    (sbm_soil &&  clust_land_use)|
    (sbm_soil &&  clust_landscp)|
    (Vegetation && log_biomass)
  #  with(clust_landscp, 2) | with(sbm_plant_grid, 2) | with(clust_land_use, 2)|
  #  with(sbm_soil, 2)
)
)

dd_selection = dredge(full.model1,subset = sexpr)

length(attr(dd_selection, "model.calls"))

subset(dd_selection, delta < 5)%>%
  as.data.frame()
# library(kableExtra)
# subset(dd_selection, delta < 1)%>%
#   as.data.frame()%>%
#   knit_with_parameters()
#   map(~.x) %>%
#   discard(~all(is.na(.x))) %>%
#   map_df(~.x)%>%
#   mutate(across(where(is.numeric), ~round(.,2)))%>%
#   knitr::kable("latex",booktabs = TRUE)%>%
#   kable_styling(font_size=4)
#   writeLines("writing/images/fig_vir/table_selection_rich.tex")


summary(get.models(dd_selection, 1)[[1]])

lapply(get.models(dd_selection, c(1,2,3,4,5)), summary)

best_model = get.models(dd_selection, 1)[[1]]
best_model

#glm(formula = log_virus_richness ~ log_plant_richness + Region_code, data = df)

best_model$deviance/best_model$df.residual


dev_model <- best_model$deviance
dev_nul <- best_model$null.deviance
(dev_nul-dev_model)/dev_nul 

#part de variance expliquée 0.3476266
par(mfrow = c(2,2))
plot(best_model)
par(mfrow = c(1,1))
hist(best_model$residuals, breaks = 30)

##### second best model
second_model = get.models(dd_selection, 2)[[1]]
summary(second_model)
#glm(formula = log_virus_richness ~ log_plant_richness + Region_code, data = df)

second_model$deviance/second_model$df.residual

dev_model2 <- second_model$deviance
dev_nul2 <- second_model$null.deviance
(dev_nul2-dev_model2)/dev_nul2
#part de variance expliquée 0.3476266
par(mfrow = c(2,2))
plot(second_model)
par(mfrow = c(1,1))
hist(second_model$residuals, breaks = 20)


interact_plot(best_model, pred = log_biomass  , modx = sbm_soil ,
              interval = F, plot.points = TRUE, data = df, partial.residuals = T) -> interaction1

effect_plot(best_model, pred = sbm_soil , interval = TRUE,
            plot.points = TRUE, data = df, partial.residuals = T) -> interaction1_soil

pwc <- df %>% 
  rstatix::emmeans_test(
    log_virus_richness ~ sbm_soil, 
    model = best_model,
    p.adjust.method = "bonferroni"
  )
pwc
rstatix::get_emmeans(pwc)

plot_m1_landscape = interaction1_landscape +
  # annotate("text", y = 3.9,x=  c(1,2,3,4), label = c("a", "b", "a", "b"),
  #           size =5)+
  
  labs(x = "Classification paysages", y = "Virus richness (log)")+
  scale_y_continuous(limits = c(0,4))

plot_m1_landscape


interact_plot(second_model, pred = log_plant_H, modx = log_biomass, 
              interval = F, plot.points = TRUE, data = df, partial.residuals = T) -> interaction2
plot_m2_biom_rich = interaction2+
  scale_color_viridis(discrete=F, option="plasma")+
  scale_x_continuous(transform = "exp")+
  scale_y_continuous(transform = "exp")+
  labs(x = "Plant diversity (log Shannon)", y = "Virus richness (log)",
       col = "Biomass (log)", linetype ="Biomass fixed (log)", group  = "Biomass fixed (log)")
plot_m2_biom_rich

plot_grid(plot_m1_biom_rich, plot_m1_landscape,
          labels = c('A', 'B'),ncol = 2, label_size = 14)
plot_m2_biom_rich

##### Description Classification landscape

read.table("data/data_clean/Metadata_grid_grp_env.txt", header = T)%>%
  dplyr::select(clust_landscp ,Grid_code)%>%
  left_join(metadata_grid%>%
              dplyr::select(artificial, cultivated,
                            natural_landscape, non_emitting, wetland, Grid_code),
            by = "Grid_code")-> details_landscape
details_landscape%>%
  dplyr::select(-Grid_code)%>%
  pivot_longer(-clust_landscp)%>%
  group_by(clust_landscp, name)%>%
  summarise(means = mean(value),
            ecartype = sd(value),
            mins = min(value),
            maxs = max(value),
            count = n())%>%
  ungroup()%>%
  pivot_longer(-c(clust_landscp, name), values_to = "values", names_to = "type")%>%
  pivot_wider(names_from = name,
              values_from = values )%>%
  mutate(across(where(is.numeric),~round(.,2)))%>%
  knitr::kable()
