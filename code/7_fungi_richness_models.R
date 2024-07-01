# Packages and data ----
library(MASS)
library(MuMIn)
library(car)
library("jtools")
library(interactions)
library(cowplot)

library(emmeans)
library(viridis)

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
read.table("data/data_clean/Metadata_grid_CAM.txt", header = T) -> metadata_grid
metadata_grid%>%
  mutate(clust_landscp = as.factor(clust_landscp), 
         clust_land_use = as.factor(clust_land_use),
         sbm_soil = as.factor(sbm_soil),
         clust_smb_grid_plant = as.factor(clust_smb_grid_plant)) -> metadata_grid

# Model selection ----
## full model ----


full.model1 <- glm(log(fungi_richness) ~ (Vegetation +log_biomass  )*(log_plant_richness +log_plant_H + plant_equit +log_plant_richness2)+
                     (clust_landscp + clust_land_use + sbm_soil + clust_smb_grid_plant), data = metadata_grid,
                   na.action = "na.fail")

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
    (log_plant_richness2 && log_plant_H)&
    dc(log_plant_richness, log_plant_richness2)
)
)

## selection----
dd_selection = dredge(full.model1,subset = sexpr)

length(attr(dd_selection, "model.calls")) # n mods

subset(dd_selection, delta < 5)%>%
  as.data.frame()

get.models(dd_selection, 1)[[1]]
summary(get.models(dd_selection, 1)[[1]])

# 5 best
lapply(get.models(dd_selection, c(1,2,3,4,5)), summary) 

# best one
best_model = get.models(dd_selection, 1)[[1]]
summary(best_model)

best_model$deviance/best_model$metadata_grid.residual

dev_model <- best_model$deviance
dev_nul <- best_model$null.deviance
(dev_nul-dev_model)/dev_nul # R² approx

par(mfrow = c(2,2))
plot(best_model)
par(mfrow = c(1,1))
hist(best_model$residuals, breaks = 30)

library(pals)
interact_plot(best_model, pred = log_plant_richness,modx = clust_smb_grid_plant ,
              interval = F, plot.points = TRUE, data = metadata_grid, partial.residuals = T,
              point.size = 5,
              colors = brewer.set1(7),
              legend.main ="Classification végtale"
              #legend.main = "Biomasse fixée (log)",
) -> interaction1
plot_m1_biom_rich = interaction1+
  #scale_x_continuous(transform = "exp")+
  #scale_y_continuous(transform = "exp")+
  labs(x = "Richesse végétale (log)", y = "Richesse fongique (log)")+
  main_theme
plot_m1_biom_rich

pwc <- metadata_grid %>% 
  rstatix::emmeans_test(
    formula = log(fungi_richness) ~ clust_smb_grid_plant, 
    model = best_model,
    p.adjust.method = "BH"
  )
pwc
rstatix::get_emmeans(pwc)

###
read.table(file = "data/data_clean/abund_plant_grid.txt", header = T)%>%
  filter(!str_detect(Grid_code, "22_CAM_15"))%>%
  column_to_rownames("Grid_code")%>%
  t() -> otu_plant_grid
reshape2::melt(otu_plant_grid)%>%
  dplyr::rename(plant = "Var1", Grid_code = "Var2" )%>%
  full_join(metadata_grid, by = "Grid_code")-> plant_profil

which(is.na(plant_profil$clust_smb_grid_plant))
plant_profil%>%
  dplyr::select(plant, value,clust_smb_grid_plant  )%>%
  filter(value>0)%>%
  
  group_by(clust_smb_grid_plant , plant)%>%
  summarise(n = sum(value))%>%
  group_by(clust_smb_grid_plant )%>%
  add_tally(n)%>%
  ungroup()%>%
  group_by(clust_smb_grid_plant ,plant)%>%
  arrange(desc(n))%>%
  group_by(clust_smb_grid_plant )%>%
  slice_head( n = 10)%>%
  mutate(prop = n/nn) -> plantsp_df
plantsp_df%>%tail
plantsp_df %>%
  mutate(plant =  str_replace_all(plant,"[.]", " "),
         plant = str_replace_all(plant,"sp", "sp. "),
         prop = round(prop,3)*100,
         plant =paste(plant," (", prop, "%)"))%>%
  group_by(clust_smb_grid_plant)%>%
  mutate(id_row = 1:10)%>%
  ungroup()%>%
  dplyr::select(clust_smb_grid_plant , plant, id_row)%>%
  pivot_wider(values_from = plant, names_from =id_row )%>%
  as.data.frame()%>%
  t()%>%
  as_tibble()%>%
  tail(10)%>%
  knitr::kable()
