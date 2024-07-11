# Packages ----
library(readxl)
library(tidyverse)
library(iNEXT)
library(sbm)
# Import and format data ----
read_xlsx("data/S1_Analyse_20CAM.xlsx")%>%
  dplyr::select(Host_code, Species) %>%
  distinct() %>% 
  mutate(Presence = 1) %>% 
  distinct()%>% 
  pivot_wider(names_from = Host_code, values_from = Presence, values_fill = 0)%>%
  filter(!is.na(Species))%>%
  column_to_rownames("Species")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("Host_code")-> otu_virus
#richness per replicat ---- 
otu_virus%>%
  reframe(Host_code = Host_code,
    richness = rowSums(across(where(is.numeric))),
    quadra_code =str_extract(Host_code, ".._CAM_...."),
    Grid_code =str_extract(Host_code, ".._CAM_.."))%>%
    group_by(quadra_code)%>%
  mutate(richness_med = median(richness))%>%
  ungroup()-> richness_df
richness_df%>%
  group_by(quadra_code)%>%
  summarise(richness_med = median(richness))%>%
  arrange(richness_med)%>%pull(quadra_code) -> order_quadra

richness_df%>%group_by(richness_med)%>%
  summarise(richness_var_by_med = var(richness))
richness_df%>%
  mutate(quadra_code = factor(quadra_code, levels = order_quadra))%>%
  ggplot()+
  facet_wrap(~richness_med)+
  geom_boxplot(aes(quadra_code, richness), outliers = F)+
  geom_jitter(aes(quadra_code, richness),width = 0.1, height = 0.05)
  
otu_virus%>%column_to_rownames("Host_code")%>%rowSums
otu_virus%>%
  mutate( quadra_code = str_extract(Host_code, ".._CAM_...."))%>%
  group_by(quadra_code) %>%
  summarise_if(is.numeric, sum)%>%
  column_to_rownames("quadra_code")%>%
  t()%>%as.data.frame%>%
  dplyr::select(where(~sum(.)!=0))%>%
  t()%>%as.data.frame%>%
  rownames_to_column("quadra_code")%>%
  pivot_longer(-quadra_code, names_to = "Virus", values_to = "occur")%>%
  group_by(quadra_code)%>%
  filter(sum(occur)!=0)%>%
  pivot_wider(values_from = occur, names_from = quadra_code )%>%
  column_to_rownames("Virus") -> OUT_virus_grid_binary_df
str(OUT_virus_grid_binary_df)
table(as.matrix(OUT_virus_grid_binary_df))
rbind(n = rep(5,ncol(OUT_virus_grid_binary_df)),OUT_virus_grid_binary_df) -> inext_formatted_data1
str(OUT_virus_grid_binary_df)
otu_virus$quadra_code


nrow(otu_virus)

OUT_virus_grid_binary_df%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("Host_code")%>%
  pivot_longer(-Host_code)%>%
  group_by(Host_code, value)%>%
  summarise(sum(value))

prob_true_value = OUT_virus_grid_binary_df/5
prob_true_value[prob_true_value == 0] =1
str(prob_true_value)
prob_true_value%>%
  as.matrix()%>%
plotMyMatrix(dimLabels = c("Virus", "Ech"), plotOptions = list(rowNames = T,colNames = F))+
  theme(element_text(angle = 90, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))
# INext estimation ----

if(file.exists("outputs/inext_estimates_of_viral_diversity_20CAM.Rdata")){
  load("outputs/inext_estimates_of_viral_diversity_20CAM.Rdata")
  
}else{
  z <- iNEXT(inext_formatted_data, q=c(0), datatype="incidence_freq", se=TRUE)
  z
  save(inext_formatted_data, z, file = "outputs/inext_estimates_of_viral_diversity_20CAM.Rdata")
}
str(z)
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
  main_theme

# Error probabilities ----
str(otu_virus)
table(as.matrix(OUT_virus_grid_binary_df))
OUT_virus_grid_binary_df%>%

  as.matrix()%>%
  plotMyMatrix( dimLabels = c("Virus", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

plot(0:5, dbinom(0:5,5,0.05))
points(0:5, dbinom(0:5,5,0.5))
prod(sapply(1:3,function(x) x))
apply(as.matrix(OUT_virus_grid_binary_df),1:2, function(x)ifelse(x==0,1, prod(sapply(0:x,
                                                                  function(y) dbinom(y, x, 0.5)))))-> alpha_error
alpha_error%>% # probabilité de détecter un virus sachant qu'il n'y en a pas avec une erreur d'echantillonage supposée de 10 %
  as.matrix()%>%
  plotMyMatrix( dimLabels = c("Virus", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

apply(-(as.matrix(OUT_virus_grid_binary_df)-5),1:2, function(x)ifelse(x==0,1, prod(sapply(0:x,
                                                                                      function(y) dbinom(y, x, 0.5))))) -> beta_error
beta_error%>%
  as.matrix()%>%
  plotMyMatrix( dimLabels = c("Virus", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))
full_error = alpha_error*beta_error

full_error%>%
  as.matrix()%>%
  plotMyMatrix( dimLabels = c("Virus", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

otu_virus%>%
  column_to_rownames("rep_code")%>%
  as.matrix()%>%
  t()%>%
  plotMyMatrix(dimLabels = c("Virus", "Ech"), plotOptions = list(rowNames = T,colNames = F))+
  theme(element_text(angle = 90, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))
OUT_virus_grid_binary_df%>%
  as.matrix()%>%
  t()%>%
  estimateBipartiteSBM(
    model = 'poisson', 
    dimLabels = c(row = "Quadrats", col = "Virus")) -> test_sbm
plot(test_sbm)

## prop ---
n_sp = nrow(OUT_virus_grid_binary_df)
OUT_virus_grid_binary_df%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("Host_code")%>%
  pivot_longer(-Host_code)%>%
  mutate(value = as.factor(value))%>%
  filter(value %in% c(4,5))%>%
  group_by(Host_code, value)%>%

  summarise(prop_value = n())%>%
  ggplot()+
  geom_col(aes(Host_code,prop_value,fill = value))

OUT_virus_grid_binary_df%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("Host_code")%>%
  mutate(Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  pivot_longer(-c(Grid_code,Host_code))%>%
  filter(value %in% c(0,1,2,3,4,5))%>%
  group_by(Grid_code,name)%>%
  filter(sum(value)!=0)%>%
  group_by(Grid_code)%>%
  mutate(n_sp_grid = n_distinct(name))%>%
  ungroup()%>%
  filter(value %in% c(1,2,3,4,5))%>%
  group_by(Host_code, value)%>%
  summarise(prop_value = n()/n_sp_grid,
            n_sp_grid = n_sp_grid)%>%
  distinct()%>%
  mutate(Grid_code = str_extract(Host_code, ".._CAM_.."))-> temp
  
temp%>%group_by(Grid_code)%>%summarise(n_grid = unique(n_sp_grid))
library(RColorBrewer)
ggplot(temp)+
  facet_wrap(~Grid_code)+
  geom_col(aes(Host_code,prop_value,fill = as.factor(value)))+
  scale_fill_manual(values = c("lightblue", "blue", "darkblue", "orange","pink"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"))


OUT_virus_grid_binary_df%>%
  as.data.frame()%>%
  rownames_to_column("virus")%>%
  pivot_longer(-virus)%>%
  filter(value %in% c(0,1,2,3,4,5))%>%
  group_by(virus,name)%>%
  filter(sum(value)!=0)%>%
  group_by(virus)%>%
  mutate(n_sp_grid = n_distinct(name))%>%
  ungroup()%>%
  filter(value %in% c(1,2,3,4,5))%>%
  group_by(virus, value)%>%
  summarise(prop_value = n()/n_sp_grid,
            n_sp_grid = unique(n_sp_grid),
            name = unique(name))%>%
  distinct()%>%
  mutate(Grid_code = str_extract(name, ".._CAM_.."))-> temp2


ggplot(temp2)+
  facet_wrap(~Grid_code)+
  geom_col(aes(prop_value,virus,fill = as.factor(value)))+
  scale_fill_manual(values = c("lightblue", "blue", "darkblue", "green","darkgreen"))
                    
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"))
##

  
  
## Beta div ----