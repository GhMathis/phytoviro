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

### data :
Contain all the raw data that I've used

#### data_clean 
Contains rearranged and cleaned data (mostly from `0_cleaning_data`)

#### shapefile 

Shapefile used in `1_buffer_size` (`sql_statement_d551600.shp`) coming from https://trouver.datasud.fr/dataset/parc-de-camargue-occupation-du-sol-2016/resource/a8f4ca8e-72b8-47fa-b03a-6d2566bfae5a.
And the same cropped shapefile used to make maps in others scripts (`crop_shapefile.shp`).

#### outputs

Pre-compute objects to save time while running scripts. Like estimation of plant, viral diversity, `SILAND` or Stochastic Block Models outputs. 