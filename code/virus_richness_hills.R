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

#Virus----
##Lire et formater les données ----
read.table("data/data_clean/OTU_virus_CAM.txt",header = T)%>%
  str()
read.table("data/data_clean/OTU_virus_CAM.txt",header = T)%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  summarise_if(is.numeric, sum)%>%
  pivot_longer(-Grid_code, names_to = "Virus", values_to = "occur")%>%
  group_by(Grid_code)%>%
  filter(sum(occur)!=0)%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Virus") -> OUT_virus_grid_binary_df

rbind(n = rep(9,ncol(OUT_virus_grid_binary_df)),OUT_virus_grid_binary_df) -> inext_formatted_data

##Estimations de diversité ----
#z <- iNEXT(inext_formatted_data, q=c(0,1), datatype="incidence_freq", se=TRUE)
#z
# 
#save(inext_formatted_data, z, file = "outputs/inext_estimates_of_viral_diversity.Rdata")

load("outputs/inext_estimates_of_viral_diversity.Rdata")

z$iNextEst$size_based %>% 
  filter(Order.q==0) %>% 
  mutate(qD_obs = if_else(Method == "Observed", qD, NA),
         qD_rar = if_else(Method == "Extrapolation", NA, qD),
         qD_ext = if_else(Method == "Rarefaction", NA, qD)) %>% 
  ggplot(., aes(x = t, y = qD, group = Assemblage))+
  geom_point(aes(x = t, y = qD_obs)) +
  geom_line(aes(x = t, y = qD_rar, group = Assemblage), linewidth = 0.5) +
  geom_line(aes(x = t, y = qD_ext, group = Assemblage), linetype = 2, linewidth = 0.5)+
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), alpha = 0.05, linetype = 0)+
  ylab("Virus richness")+
  xlab("Sample size")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )
z$iNextEst$size_based -> temp
z$iNextEst$size_based %>% 
  filter(Order.q==0 & Method == "Observed") %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         virus_richness = qD, 
         v_rich_LCL = qD.LCL, 
         v_rich_UCL = qD.UCL)-> virus_richness

z$iNextEst$size_based %>% 
  filter(Method == "Observed" & Order.q == 1) %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         virus_H = qD, 
         v_H_LCL = qD.LCL, 
         v_H_UCL = qD.UCL) -> virus_H

z$AsyEst %>% 
  filter(Diversity == "Species richness") %>% 
  ggplot(., aes(x=Observed, y = Estimator))+
  geom_point()


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
  left_join( plant_H, by = "Grid_code")%>%
  left_join( virus_richness, by = "Grid_code")%>%
  left_join( virus_H, by = "Grid_code")%>%
 
  mutate(virus_richness = replace(virus_richness, is.na(virus_richness), 0),
         virus_H = replace(virus_H, is.na(virus_H), 0),
         
         plant_richness2 = plant_richness^2, 
         p_richness_class5 = cut(plant_richness, breaks = seq(0,30, 5)),
         p_richness_class2 = cut(plant_richness, breaks = seq(0,30, 2)),
         log_virus_richness = log(virus_richness+1),
         log_plant_richness = log(plant_richness+1),
         log_plant_richness2 = (log(plant_richness+1))^2,
         log_plant_H = log(plant_H+1),
         log_virus_H = log(virus_H+1),
         log_biomass = log(Biomass+1),
         sbm_plant_grid = as.factor(SBM_abund_grid$memberships$Grids),
         log_plant_equit = log_plant_H-log_plant_richness, 
         plant_equit = plant_H/plant_richness) -> df

df$plant_H
df$plant_richness
df$Vegetation
df%>%
  group_by(sbm_plant_grid )%>%
  summarise(mean(plant_H),
            n(),max(plant_H), min(plant_H))
df%>%
  left_join(read.table("outputs/cluster_df.txt"), by = "Grid_code")%>%
  group_by(sbm_virus)%>%
  summarise(median(virus_richness),
            n(),max(virus_richness), min(virus_richness))
summary(df$plant_H)
summary(df$Biomass)
summary(df$virus_richness)
view(df)
#write.table(df, file = "outputs/dataframe_diversities.txt")

library(GGally)

df %>% 
  dplyr::select(virus_richness, plant_richness, plant_H, plant_equit, Biomass) %>% 
  ggpairs()

df %>% 
  dplyr::select(log_virus_richness, log_plant_richness, log_plant_H, log_plant_equit, log_biomass, Vegetation) %>% 
  ggpairs()

df %>% 
  ggplot(., aes(x= log(plant_richness+1), y = log(virus_richness+1)))+
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = log(pl_rich_LCL+1), xmax = log(pl_rich_UCL+1)), width = 0.05)+
  geom_errorbar(aes(ymin = log(v_rich_LCL+1), ymax = log(v_rich_UCL+1)), width = 0.05)+
  ylab("Virus richness (log)")+
  xlab("Plant richness (log)")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )
df %>% 
  ggplot(., aes(x= log_plant_H, y = log_virus_H))+
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = log(pl_H_LCL+1), xmax = log(pl_H_UCL+1)), width = 0.05)+
  geom_errorbar(aes(ymin = log(v_H_LCL+1), ymax = log(v_H_UCL+1)), width = 0.05)+
  ylab("Virus shannon (log)")+
  xlab("Plant shannon (log)")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )

read.table("data/data_clean/Metadata_grid_grp_env.txt", header = T)%>%
  dplyr::select(clust_landscp, clust_land_use, sbm_soil, Grid_code)%>%
  mutate(across(where(is.numeric), as.factor))-> metadata_grid_grp_env

df%>%
  left_join(metadata_grid_grp_env, by ="Grid_code") -> df
df %>% 
  ggplot(., aes(x= log(plant_richness+1), y = log(virus_richness+1)))+
  facet_wrap(~clust_land_use)+
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = log(pl_rich_LCL+1), xmax = log(pl_rich_UCL+1)), width = 0.05)+
  geom_errorbar(aes(ymin = log(v_rich_LCL+1), ymax = log(v_rich_UCL+1)), width = 0.05)+
  ylab("Virus richness (log)")+
  xlab("Plant richness (log)")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )

df %>% 
  ggplot(., aes(x= log(plant_richness+1), y = log(virus_richness+1)))+
  facet_wrap(~sbm_soil)+
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = log(pl_rich_LCL+1), xmax = log(pl_rich_UCL+1)), width = 0.05)+
  geom_errorbar(aes(ymin = log(v_rich_LCL+1), ymax = log(v_rich_UCL+1)), width = 0.05)+
  ylab("Virus richness (log)")+
  xlab("Plant richness (log)")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )

df %>% 
  ggplot(., aes(x= log(plant_richness+1), y = log(virus_richness+1)))+
  facet_wrap(~clust_landscp)+
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = log(pl_rich_LCL+1), xmax = log(pl_rich_UCL+1)), width = 0.05)+
  geom_errorbar(aes(ymin = log(v_rich_LCL+1), ymax = log(v_rich_UCL+1)), width = 0.05)+
  ylab("Virus richness (log)")+
  xlab("Plant richness (log)")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )

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

##### Virus Richness
full.model1 <- glm(log_virus_richness ~ (Vegetation +log_biomass  )*(log_plant_richness +log_plant_H + plant_equit +log_plant_richness2)+
                    (clust_landscp + clust_land_use + sbm_soil + sbm_plant_grid), data = df,
                   na.action = "na.fail")
#####
sexpr = expression(!(
  (clust_landscp && sbm_plant_grid)|
  (clust_land_use  &&  sbm_plant_grid)|
  (clust_land_use  &&  clust_landscp)|
  (sbm_soil  && sbm_plant_grid)|
  (sbm_soil &&  clust_land_use)|
  (sbm_soil &&  clust_landscp)|
  (Vegetation && log_biomass)|
  (log_plant_richness && log_plant_H)|
  (plant_equit && log_plant_H)|
  (log_plant_richness2 && log_plant_H)
#  with(clust_landscp, 2) | with(sbm_plant_grid, 2) | with(clust_land_use, 2)|
#  with(sbm_soil, 2)
  )
)

dd_selection = dredge(full.model1,subset = sexpr)

length(attr(dd_selection, "model.calls"))

subset(dd_selection, delta < 2)%>%
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
summary(best_model)
#glm(formula = log_virus_richness ~ log_plant_richness + Region_code, data = df)

best_model$deviance/best_model$df.residual
vif(best_model)

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
mean(df$log_biomass)
mean(df$log_biomass) +sd(df$log_biomass)
mean(df$log_biomass) -sd(df$log_biomass)
interact_plot(best_model, pred = log_plant_H, modx = log_biomass,
              interval = F, plot.points = T, data = df, partial.residuals = T,
              point.size = 5, 
              legend.main = "Biomasse fixée (log)",
              modx.labels = c("- 1 Sd. (4.0)",
                              "Moyenne (4.6) ",
                              " + 1 Sd.(5.1)")) -> interaction1

plot_m1_biom_rich = interaction1+
  scale_color_viridis(discrete=F, option="plasma")+
  labs(x = "Hétérogénéité du couvert végétal (log)", y = "Richesse virale (log)",
       col = "Biomasse (log)")+
  main_theme

plot_m1_biom_rich
str(plot_m1_biom_rich[[1]])
str(interaction1)

effect_plot(best_model, pred = clust_landscp, interval = TRUE,x.label = "A",
            plot.points = TRUE, data = df, partial.residuals = T,
            int.type = "confidence",
            point.size = 5) -> interaction1_landscape

pwc <- df %>% 
  rstatix::emmeans_test(
    formula = log_virus_richness ~ clust_landscp, 
    model = best_model,
    p.adjust.method = "BH"
  )
pwc
rstatix::get_emmeans(pwc)
str(interaction1_landscape[[1]])
plot_m1_landscape = interaction1_landscape +
 # annotate("text", y = 3.9,x=  c(1,2,3,4), label = c("a", "b", "a", "b"),
#           size =5)+
  scale_x_discrete(labels=  c("Cultivés dominant\n + Anthropisés",
                              "Naturels non Humide\n dominant",
                              "Mixtes", 
                              "Naturels Humide\n dominant"))+
  labs(x = "Classification du paysage", y = "Richesse virale (log)")+
  scale_y_continuous(limits = c(0,4))+
  main_theme

plot_m1_landscape


interact_plot(second_model, pred = log_plant_H, modx = log_biomass, 
              interval = F, plot.points = TRUE, data = df, partial.residuals = T,
              point.size = 5, 
              legend.main = "Biomasse fixée (log)",
              modx.labels = c("- 1 Sd. (4.0)",
                              "Moyenne (4.6) ",
                              " + 1 Sd.(5.1)")) -> interaction2
plot_m2_biom_rich = interaction2+
  #scale_color_viridis(discrete=F, option="plasma")+
  # scale_x_continuous(transform = "exp")+
  # scale_y_continuous(transform = "exp")+
  scale_color_viridis(discrete=F, option="plasma")+
  labs(x = "Divérsité végétale d'ordre 1 (log)", y = "Richesse virale (log)",
       col = "Biomasse (log)")+
  main_theme

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

# crPlots(best_model,col.lines=c("red", "blue"))
# library(effects)
# plot(predictorEffects(best_model), lines=list(multiline=TRUE))
# 
# plot(predictorEffects(best_model, partial.residuals=TRUE))
# 
# 
# df%>%
#   mutate(biom_perc = cut(ecdf(log_biomass)(log_biomass)* 100, breaks = c(0,25, 50, 75, 100),
#                          labels = c("0-25%", "25-50%", "50-75%", "75-100%")) )%>%
#   ggplot(., aes(x= log_plant_H, y = log_virus_H, col = biom_perc))+
#   geom_point(size = 2) +
#   #facet_wrap(~clust_landscp)+
#   geom_errorbar(aes(xmin = log(pl_H_LCL+1), xmax = log(pl_H_UCL+1)), width = 0.05)+
#   geom_errorbar(aes(ymin = log(v_H_LCL+1), ymax = log(v_H_UCL+1)), width = 0.05)+
#   ylab("Virus shannon (log)")+
#   xlab("Plant shannon (log)")+
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 15, face = "bold"),
#     axis.title.y = element_text(size = 15, face = "bold"),
#     legend.title = element_text(size=12, face="bold"),
#     legend.text = element_text(size=12, face="italic")
#   )
# ##### second best model
# second_model = get.models(dd_selection, 2)[[1]]
# summary(second_model)
# #glm(formula = log_virus_richness ~ log_plant_richness + Region_code, data = df)
# 
# second_model$deviance/second_model$df.residual
# vif(second_model)
# 
# dev_model <- second_model$deviance
# dev_nul <- second_model$null.deviance
# (dev_nul-dev_model)/dev_nul
# #part de variance expliquée 0.3476266
# par(mfrow = c(2,2))
# plot(second_model)
# par(mfrow = c(1,1))
# hist(second_model$residuals, breaks = 20)

# plot(predictorEffects(second_model), lines=list(multiline=TRUE))
# plot(predictorEffects(second_model, partial.residuals=TRUE))

# ##### Virus Shannon
# 
# full.model1_H <- glm(log_virus_H ~ (Vegetation +log_biomass  )*(log_plant_richness +log_plant_H)+
#                      (clust_landscp + clust_land_use + sbm_soil + sbm_plant_grid), data = df,
#                    na.action = "na.fail")
# 
# dd_selection_H = dredge(full.model1_H,subset = sexpr)
# 
# subset(dd_selection_H, delta < 5)%>%
#   as.data.frame()%>%
#   map(~.x) %>%
#   discard(~all(is.na(.x))) %>%
#   map_df(~.x)%>%
#   mutate(across(where(is.numeric), ~round(.,2)))%>%
#   knitr::kable()
# 
# length(attr(dd_selection_H, "model.calls"))
# 
# summary(get.models(dd_selection_H, 1)[[1]])
# 
# lapply(get.models(dd_selection_H, c(1,2,3,4)), summary)
# 
# get.models(dd_selection_H, 1)[[1]]
# 
# best_model_H = get.models(dd_selection_H, 1)[[1]]
# dev_model_H <- best_model_H$deviance
# dev_nul_H <- best_model_H$null.deviance
# (dev_nul_H-dev_model_H)/dev_nul_H
# summary(best_model_H)
# plot(predictorEffects(best_model_H), lines=list(multiline=TRUE))
##### rerpresentation des 2 modéles finaux
# 
# 
# #####Richesse
# 
# effect_plot(best_model, pred = clust_landscp, interval = TRUE, plot.points = TRUE, data = df, partial.residuals = T)
# 
# interact_plot(best_model, pred = log_biomass, modx = log_plant_H, mod2 = clust_landscp,
#               interval = F, plot.points = TRUE, data = df, partial.residuals = T) -> interaction1
# interaction1
# interact_plot(best_model,pred = log_plant_H, modx =log_biomass,mod2 = clust_landscp,
#               interval = F, plot.points = TRUE, data = df, partial.residuals = T)-> interaction2
# richness_plot = interaction2+
#   labs(x = "Plant diversity (log Shannon)", y = "Virus richness (log)")
# 
# effect_plot(best_model, pred = clust_landscp, interval = TRUE, plot.points = TRUE, data = df, partial.residuals = T)
# ##### Shanon
# summary(best_model_H)
# effect_plot(best_model_H, pred = clust_landscp, interval = TRUE, plot.points = TRUE, data = df, partial.residuals = F)
# interact_plot(best_model_H,pred = log_biomass, modx = log_plant_H, mod2 =clust_landscp ,
#               interval = F, plot.points = TRUE, data = df, partial.residuals = T)-> interaction1_H
# interaction1_H
# interact_plot(best_model_H,pred = log_plant_H, modx =log_biomass , mod2 =clust_landscp ,
#               interval = F, plot.points = TRUE, data = df, partial.residuals = T)-> interaction2_H
# interaction2_H
# 
# shannon_plot = interaction2_H+
#   labs(x = "Plant diversity (log Shannon)", y = "Virus diversity (log Shannon)")
# 
# plot_grid(richness_plot, shannon_plot, labels = c('A', 'B'),nrow = 2, label_size = 14)
# 
# summary(best_model)%>%coef()%>%
#   round(2)%>%
#   as.data.frame()%>%
#   knitr::kable()
# summary(best_model_H)%>%coef()%>%
#   round(2)%>%
#   as.data.frame()%>%
#   knitr::kable()
# 
# 
# summary(best_model)
# (dev_nul-dev_model)/dev_nul
# summary(best_model_H)
# (dev_nul_H-dev_model_H)/dev_nul_H


##### Comparaison 

best_model
get.models(dd_selection2, 1)[[1]]
best_model_H
get.models(dd_selection_H, 1)[[1]]
summary(best_model_H)
best_model_H2 = get.models(dd_selection_H, 1)[[1]]
interact_plot(best_model_H2,pred = log_biomass, modx = clust_landscp, 
              interval = F, plot.points = TRUE, data = df, partial.residuals = T) -> interaction1_H2
interaction1_H[[1]]%>%
  cbind(y_type = "Virus_shanon")-> temp
names(temp) = c("log_plant_richness", "log_virus", "clust_landscp", "log_biomass",     
                "modx_group", "y_type")
temp = temp[,c(2,1,3,4,5,6)]

interaction1[[1]]%>%
  cbind(y_type = "Virus_rich")-> temp2
names(temp2) = c("log_virus", "log_plant_richness", "clust_landscp", "log_biomass",     
                "modx_group", "y_type")

rbind(temp, temp2) -> compar_rich_shannon

ggplot(compar_rich_shannon)+
  geom_point(aes(log_biomass, log_virus, col = clust_landscp ))+

  facet_wrap(~y_type, scales = "free")
ggplot(compar_rich_shannon)+
  geom_point(aes(log_plant_richness, log_virus, col = clust_landscp ))+
  facet_wrap(~y_type, scales = "free")

