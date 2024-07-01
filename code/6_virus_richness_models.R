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

# Relation richesse-richesse ----

library(GGally)

metadata_grid %>% 
  dplyr::select(virus_richness, plant_richness, plant_H, plant_equit, Biomass) %>% 
  ggpairs()

metadata_grid %>% 
  dplyr::select(log_virus_richness, log_plant_richness, log_plant_H, log_plant_equit, log_biomass, Vegetation) %>% 
  ggpairs()

metadata_grid %>% 
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
metadata_grid %>% 
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

metadata_grid %>% 
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

metadata_grid %>% 
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

metadata_grid %>% 
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

# Models selection ----

library(MASS)
library(MuMIn)
library(car)
library("jtools")
library(interactions)
library(cowplot)

library(emmeans)
library(viridis)
library(RColorBrewer)
#On transforme virus_richness en log et on tente le glm gaussien 

## Virus Richness ----
metadata_grid%>%
  mutate(clust_landscp = as.factor(clust_landscp), 
         clust_land_use = as.factor(clust_land_use),
         sbm_soil = as.factor(sbm_soil),
         clust_smb_grid_plant = as.factor(clust_smb_grid_plant)) -> metadata_grid
full.model1 <- glm(log_virus_richness ~ (Vegetation +log_biomass  )*(log_plant_richness +log_plant_H + plant_equit +log_plant_richness2)+
                     (clust_landscp + clust_land_use + sbm_soil + clust_smb_grid_plant), data = metadata_grid,
                   na.action = "na.fail")
# Conditions
sexpr = expression(!(
  (clust_landscp && clust_smb_grid_plant)|
    (clust_land_use  &&  clust_smb_grid_plant)|
    (clust_land_use  &&  clust_landscp)|
    (sbm_soil  && clust_smb_grid_plant)|
    (sbm_soil &&  clust_land_use)|
    (sbm_soil &&  clust_landscp)|
    (Vegetation && log_biomass)|
    (log_plant_richness && log_plant_H)|
    (plant_equit && log_plant_H)|
    (log_plant_richness2 && log_plant_H) &
    dc(log_plant_richness, log_plant_richness2)
)
)
 
#run the selection
dd_selection = dredge(full.model1,subset = sexpr)

length(attr(dd_selection, "model.calls")) # n models

# Best models
subset(dd_selection, delta < 2)%>%
  as.data.frame()

#The best one
summary(get.models(dd_selection, 1)[[1]])

lapply(get.models(dd_selection, c(1,2,3,4,5)), summary)

best_model = get.models(dd_selection, 1)[[1]]
best_model
summary(best_model)

best_model$deviance/best_model$df.residual


dev_model <- best_model$deviance
dev_nul <- best_model$null.deviance
(dev_nul-dev_model)/dev_nul #R² approx

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
mean(metadata_grid$log_biomass)
mean(metadata_grid$log_biomass) +sd(metadata_grid$log_biomass)
mean(metadata_grid$log_biomass) -sd(metadata_grid$log_biomass)

interact_plot(best_model, pred = log_plant_H, modx = log_biomass,
              interval = F, plot.points = T, data = metadata_grid, partial.residuals = T,
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


effect_plot(best_model, pred = clust_landscp, interval = TRUE,x.label = "A",
            plot.points = TRUE, data = metadata_grid, partial.residuals = T,
            int.type = "confidence",
            point.size = 5) -> interaction1_landscape

pwc <- metadata_grid %>% 
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
              interval = F, plot.points = TRUE, data = metadata_grid, partial.residuals = T,
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

## Description Classification landscape ----

metadata_grid%>%
  dplyr::select(artificial, cultivated, clust_landscp, 
                natural_landscape, non_emitting, wetland, Grid_code)-> details_landscape
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

details_landscape%>%
  pivot_longer(-c(Grid_code,clust_landscp), values_to = "landscp_val", names_to = "landscp_name") %>%
  ggplot() +
  facet_wrap(~clust_landscp) +
  geom_col(aes(x = Grid_code,y =landscp_val, fill = landscp_name))+
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(4,2,1,3,5)],
                    labels = c("Antropisés", "Cultivés", "Naturels non humide",
                               "Non emmeteur de propagules", "Naturel humide")) +
  labs(colour = "Grid cluster", fill ="Types d'occupation des sols", x =NULL, y = NULL)+
  main_theme+
  theme(axis.text.x = element_text(colour = "black", size=6,angle = 90, vjust= 1,
                                   hjust = 1))+
  labs(x = "Grilles", y = "Portion paysages",fill = "Types de paysages")

