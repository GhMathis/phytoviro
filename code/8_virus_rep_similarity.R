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
  rownames_to_column("rep_code")-> otu_virus


otu_virus%>%
  mutate( Host_code = str_extract(rep_code, ".._CAM_...."))%>%
  group_by(Host_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  summarise_if(is.numeric, sum)%>%
  column_to_rownames("Host_code")%>%
  t()%>%as.data.frame%>%
  dplyr::select(where(~sum(.)!=0))%>%
  t()%>%as.data.frame%>%
  rownames_to_column("Host_code")%>%
  pivot_longer(-Host_code, names_to = "Virus", values_to = "occur")%>%
  group_by(Host_code)%>%
  filter(sum(occur)!=0)%>%
  pivot_wider(values_from = occur, names_from = Host_code )%>%
  column_to_rownames("Virus") -> OUT_virus_grid_binary_df

rbind(n = rep(5,ncol(OUT_virus_grid_binary_df)),OUT_virus_grid_binary_df) -> inext_formatted_data
str(OUT_virus_grid_binary_df)
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
apply(as.matrix(OUT_virus_grid_binary_df),1:2, function(x) dbinom(as.numeric(x),3,0.8))%>%
  as.matrix()%>%
  plotMyMatrix( dimLabels = c("Quadrats", "Virus"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

apply(as.matrix(OUT_virus_grid_binary_df),1:2, function(x) dbinom(as.numeric(x),3,0.8))%>%
  as.matrix()%>%
  plotMyMatrix( dimLabels = c("Quadrats", "Virus"), plotOptions = list(rowNames = T,colNames = F))+
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
plotMyMatrix(as.matrix(otu_virus), dimLabels = c("Quadrats", "Virus"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))
