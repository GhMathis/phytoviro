# prairiome

## How to use it

### Code :

Contain all scripts. All script can be run independently

-   `0_cleaning_data` compile and rearrange most of the data used in the 7 others scripts.

-   `1_buffer_size` compute buffer size with the `SILAND` package, when extract percentages of land cover in this buffer around each grid.

-   `2_environmental_structure` compute classification for all environmental variables (soil, landscape, agricultural practice)

-   `3_plant_community_structure` SBM for plant community at grid level. Then congruence analysis with NMI and variance partitioning

-   `4_virus_community_structure` same but for virus

-   `5_fungi_community_structure` same but for fungi

-   `6_virus_richness_models` model selection for virus richness at grid level

-   `7_fungi_richness_models` same but for fungi

-   `8_virus_replicat_exploration.Rmd` data exploration quadra replicates for virus 

-   `functions_modified.R` some functions used in other script

### data :

Contain all the raw data that I've used

#### data_clean

Contains rearranged and cleaned data (mostly from `0_cleaning_data`)

-   `abund_plant_grid.txt` : Plant community matrix with vegetale cover heterogeneity per grid

-   `Metadata_grid_CAM.txt` : Metadata at grid level

    -   `Grid_code` : grid ID
    -   `Collection_date`
    -   `Ecosystem` : ecosystem identified during the data collection
    -   `Country`
    -   `Locality`
    -   `X` : Latitude
    -   `Y` : longitude
    -   `Num` : ?
    -   `pHwater` : Soil pH
    -   `lime_tot` : Total lime in the soil (%)
    -   `MO` : Organic material in the soil (%)
    -   `Phos` : $P_2O_5$ (mg/kg)
    -   `K` : $K_2O$ (mg/kg)
    -   `Mg` : $MgO$ (mg/kg)
    -   `Ca` : $CaO$ (mg/kg)
    -   `Na` : $Na_2O$ (mg/kg)
    -   `N` : nitrogen (mg/kg)
    -   `C` : carbone (g/kg)
    -   `CN` : $C/N$
    -   `clay` : (%)
    -   `SiltF` : Silt fine (%)
    -   `SiltC` : Silt coarse (%)
    -   `SandF` : (%)
    -   `SandC` : (%)
    -   `Cl` : $Cl-$
    -   `Res` : resistivity (ohm/com)
    -   `Cond` : conductivity (mS/m)
    -   `Profondeur` : oxidation depth (character) (cm)
    -   `dep_oxy_num` : oxidation depth (numerical) (cm)
    -   `depth_oxy` : oxidation depth (factor class) (cm)
    -   `Pature` : Pasture (binary)
    -   `Fauche` : mowing (binary)
    -   `plant_richness`
    -   `pl_rich_LCL` : plant richness lower confidence interval
    -   `pl_rich_UCL` : plant richness upper confidence interval
    -   `plant_H` : vegetale cover heterogeneity (*Chao1* equivalent)
    -   `pl_H_LCL` : vegetale cover heterogeneity lower confidence interval
    -   `pl_H_UCL` : vegetale cover heterogeneity upper confidence interval
    -   `virus_richness`
    -   `v_rich_LCL` : virus richness lower confidence interval
    -   `v_rich_UCL` : virus richness upper confidence interval
    -   `virus_H` : virus diversity (*Chao1* approximation)
    -   `v_H_LCL` : virus diversity lower confidence interval
    -   `v_H_UCL` : virus diversity upper confidence interval
    -   `Biomass` : Biomass sum across 9 quadra for each grid (g)
    -   `Vegetation` : Vegetation cover proportion across 9 quadra for each grid
    -   `fungi_richness`
    -   `plant_richness2` : $`plant_richness`^2$
    -   `p_richness_class5` : 5 classes of richness
    -   `p_richness_class2` : 2 classes of richness
    -   `log_virus_richness`
    -   `log_plant_richness`
    -   `log_plant_richness2`
    -   `log_plant_H`
    -   `log_virus_H`
    -   `log_biomass`
    -   `log_plant_equit` : log of plant equity (`log_plant_H`-`log_plant_richness`)
    -   `plant_equit` : plant equity `plant_H`/`plant_richness`
    -   `buffer_size` : radius used for landscape buffer
    -   `wetland` : type of landscape (%)
    -   `natural_landscape` : // (%)
    -   `non_emitting` : // (%)
    -   `cultivated` : // (%)
    -   `artificial` : // (%)
    -   `clust_landscp` : Landscape cluster compute with PCA -\> distance in the 2 main axis -\> hierarchical classification
    -   `clust_land_use` : Agricultural usages classification (`Fauche` and `Pature`)
    -   `sbm_soil` : Soil type classification (with SBM)
    -   `clust_smb_grid_plant` : Plant community classification (with SBM)
    -   `sbm_virus` : Virus community classification (with SBM)
    -   `sbm_fungi` : Fungi community classification (with SBM)

-   `OTU_fungi.txt` : Community matrix for fungi at grid level

-   `OTU_fungi_for_sbm.txt` : Community matrix for fungi with fungi removed when abundances are less than 10 at grid level (to compute SBM faster)

-   `OTU_plant_CAM.txt` : Plant community matrix with cover per quadra

-   `OTU_virus_CAM.txt` : Virus community matrix with presence/absence per quadra

#### shapefile

Shapefile used in `1_buffer_size` (`sql_statement_d551600.shp`) coming from <https://trouver.datasud.fr/dataset/parc-de-camargue-occupation-du-sol-2016/resource/a8f4ca8e-72b8-47fa-b03a-6d2566bfae5a>. And the same cropped shapefile used to make maps in others scripts (`crop_shapefile.shp`).

#### outputs

Pre-compute objects to save time while running scripts. Like estimation of plant, viral diversity, `SILAND` or Stochastic Block Models outputs.
