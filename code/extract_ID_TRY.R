library(readxl)
library(tidyverse)
library(rtry)
try_data_base = read.delim("data/TryAccSpecies.txt", header = T)
tax_plant =read_xlsx("data/TAX_Plant.xlsx")
AM
str(tax_plant)
try_data_base
str(try_data_base)
tax_plant$Plant_species
tax_plant %>%
  mutate(Sp_names = str_extract(Plant_species, "^\\w+\\s\\w+")) -> tax_plant


tibble(try_data_base)%>%
  filter(AccSpeciesName %in% tax_plant$Sp_names) -> identified_sp

write.table(identified_sp$AccSpeciesID, row.names = F, col.names = F, file = "data/IDsp_TRY_EDGG.txt",  eol =", " )



rtry_import("data/life_span/31810.txt")%>%
  filter(OriglName %in% c("Lifespan", "longevity")) -> life_span



##################

str_replace_all(tax_plant$Plant_species[1], "[ -]", ".")

read_xlsx("data/TAX_plant_CAM_HF.xlsx")%>%
  mutate( Plant_species= str_replace_all(species, "[ -]", "."))%>%
  select(Plant_species, CLASSE, ORDRE, FAMILLE, SOUS_FAMILLE, TRIBU, RANG)%>%
  rename(Plant_class = "CLASSE", Plant_order = "ORDRE", Plant_family ="FAMILLE",
         Plant_sub_family = "SOUS_FAMILLE", Plant_tribu = "TRIBU", Plant_rank = "RANG")-> tax_plant_Cam_temp

OTU_plant = read.table("data/data_clean/OTU_plant_CAM.txt")
tax_plant%>%
  mutate(Plant_species = str_replace_all(Plant_species, "[ -]", "."))%>%
  filter(Plant_species %in% names(OTU_plant))%>%
  select(Plant_species, Plant_genus, Sp_names)%>%
  left_join(tax_plant_Cam_temp, by = "Plant_species")-> tax_plant_CAM






tibble(try_data_base)%>%
  filter(AccSpeciesName %in% tax_plant_CAM$Sp_names) -> identified_sp_CAM

rtry_import("data/multiple_traits/31758.txt")-> traits

c(6, 7, 8, 13, 14, 15, 22, 26, 30, 33, 37, 40, 42, 47, 128, 129, 145, 196, 231, 320, 385, 1254, 3106, 3117, 3364) -> trait

traits %>% 
  filter(AccSpeciesID %in% identified_sp_CAM$AccSpeciesID, TraitID %in% trait)%>%
  dplyr::select(AccSpeciesName, OriglName,TraitID, OrigValueStr)->traits_temp
str(traits)
unique(traits_temp$OrigValueStr)

traits_temp%>%
  group_by(TraitID,OriglName)%>%
  summarise(n = n(),first = first(OrigValueStr))%>%
  filter(n == max(n))%>%
  ungroup()%>%
  mutate(clean_name = c("Rooting_depth", "Mycorrhiza_Type", "Nitrogen_fixaxtion_cap",                                         
                         "LCC", "LNC", "LPC",                                              
                         "Plant_photosynthetic_pathway", "Seed_Weight", "Drought_tolerance",                                
                         "seed_longevity", "Leaf_phenology", "Photosynthesis_rate_per_leaf_dry_mass",                                                 
                         "Plant_Growth_Form", "LDMC", "SDM",                                            
                         "Dry_mass_leaves", "	Leaf_width", "GRIME",                                            
                         "dispercal_unit_type", "Vital_attributes_of_persistence_and_establishment", "Salinity_tolerance",                                         
                         "Plant_height" , "SLA", "photosynthetic_rate_per_leaf_area	")
         )%>%
  select(TraitID, clean_name) -> ID_names_clean

traits_temp%>%
  full_join(ID_names_clean, by = join_by(TraitID))%>%
  filter(clean_name == "Salinity_tolerance")%>%
  mutate(OrigValueStr_modif = case_when(
    OrigValueStr %in% c("Higth", "halophyte", "saline", "brackish", "800", "400","25000" ,"12000",
                        "56", "52", "77", "34", "40.6", "43", "100", "118", "75", "SW") ~ "High",
    OrigValueStr %in% c("Medium", "salt-affected", "25", "225", "200", "300", "160") ~ "Medium",
    OrigValueStr %in% c("Low", "fresh") ~ "Low",
    .default = OrigValueStr))%>%
  
  select(-c(OriglName, TraitID ,OrigValueStr))%>%
  filter(!duplicated(AccSpeciesName)) -> traits_clean
head(traits_clean)

traits_temp%>%
  full_join(ID_names_clean, by = join_by(TraitID))%>%
  filter(clean_name == "GRIME")%>%
  mutate(OrigValueStr_modif = case_when(
    OrigValueStr %in% c("c","C") ~ "C",
    OrigValueStr %in% c("r","R") ~ "R",
    OrigValueStr %in% c("s","S") ~ "S",

    OrigValueStr %in% c("sr","SR", "RS") ~ "RS",
    OrigValueStr %in% c("cs","CS","SC") ~ "CS",
    OrigValueStr %in% c("cr","CR","RC") ~ "CS",
    OrigValueStr == "csr" ~ "CRS",
    .default = OrigValueStr))%>%
  
  select(-c(OriglName, TraitID ,OrigValueStr))%>%
  filter(!duplicated(AccSpeciesName)) -> GRIME_clean

traits_temp%>%
  full_join(ID_names_clean, by = join_by(TraitID))%>%
  filter(clean_name == "Plant_photosynthetic_pathway")%>%
  mutate(OrigValueStr_modif = case_when(
    OrigValueStr %in% c("C3", "c3", "3","C3?") ~ "C3",
    OrigValueStr %in% c("C4", "c4", "C4?") ~ "C4",
    AccSpeciesName =="Sonchus maritimus" ~ "C4",
    .default = OrigValueStr))%>% 
  
  select(-c(OriglName, TraitID ,OrigValueStr))%>%
  filter(!duplicated(AccSpeciesName)) ->photo_path

unique(photo_path$OrigValueStr_modif)
##################


unique(life_span$OrigValueStr)
life_span%>%
  mutate(OrigValueStr_modif = case_when(
    OrigValueStr %in% c("perennial", "Moderate", "Long", "25 - 150", "5-25", "2 - 5", "> 150", "Perennial") ~ "Long",
    OrigValueStr %in% c("annual", "Annual", "<2", "Short") ~ "Short",
    .default = OrigValueStr
    ),
    clean_name = "life_style"
  )%>% 
  filter(AccSpeciesID %in% identified_sp_CAM$AccSpeciesID)%>%
  select(AccSpeciesName, clean_name, OrigValueStr_modif)%>%
  filter(!duplicated(AccSpeciesName))-> life_span_CAM_clean

traits_clean = rbind(traits_clean,life_span_CAM_clean)
traits_clean = rbind(traits_clean,GRIME_clean)
traits_clean = rbind(traits_clean,photo_path)



#####
tax_plant_CAM%>%
  full_join(traits_clean, by = join_by("Sp_names" == "AccSpeciesName"), relationship = "many-to-many")%>%
  mutate(OrigValueStr_modif  = case_when(OrigValueStr_modif == "None" ~ NA ,
                                         .default = OrigValueStr_modif
                                         )
         )%>%
  pivot_wider(names_from = clean_name, values_from = OrigValueStr_modif)%>%
  select(-`NA`)-> tax_plant_CAM_clean


write.table(tax_plant_CAM_clean, "data/data_clean/TAX_plant_CAM.txt")


