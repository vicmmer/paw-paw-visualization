#Creating map with our data 
# Create the map with topographic background  
library(mapmixture)
library(terra)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Create SpatRaster object
earth <- terra::rast("NE1_50M_SR_W.tif")
relief <- terra::rast("SR_50M.tif")
#Match qmat info and also coord information 
admixpp <- read.csv("admixture_and_coordinates.csv") #file with combined qmat and  coordinate data 

admixpp <- admixpp[admixpp$State != "Florida", ] #Remove Florida since precise coordinates unknown 

coord <- admixpp %>%
  select(latitude, longitude, Site) %>%
  rename(Lat = latitude, Lon = longitude) %>%
  filter(Site != "" & !is.na(Site)) %>%  # Remove entries with missing Site
  group_by(Site) %>%
  summarise(
    Lat = mean(Lat, na.rm = TRUE),  # Calculate the average Lat for each Site
    Lon = mean(Lon, na.rm = TRUE)   # Calculate the average Lon for each Site
  )
coord <- coord[coord$Site != "Florida1", ]#Remove Florida since precise coordinates unknown 

admixture <- admixpp %>%
  select(Site, corresponding.sample, V1, V2) %>% 
  rename(Ind = corresponding.sample, Cluster1 = V1, Cluster2 = V2) %>%
  filter(Site != "" & !is.na(Site)) 
 
topographic_map <-  mapmixture(admixture, 
                    coord, 
                    basemap = earth, 
                    cluster_cols =c("#0072B2", "#E69F00"),  
                    cluster_names = c("Ancestry 1", "Ancestry 2"),
                    boundary =c(xmin = -98, xmax = -71, ymin = 27, ymax = 43),  
                    pie_size = 0.8)  # Adjust pie chart size for visibility

topographic_map


##Gettign ready to plot in GIS 
#The current admixture inclides every individual/site so i want to summarize and average per site 

options(scipen = 999)  # Turn off scientific notation

admix_site <- admixture %>%
  group_by(Site) %>%
  summarise(
    Avg_Cluster1 = mean(Cluster1, na.rm = TRUE),
    Avg_Cluster2 = mean(Cluster2, na.rm = TRUE),
  )
admix_site_coords <- left_join(admix_site, coord, by = "Site")
write.csv(admix_site_coords, "admix_site_coords.csv")
#admix_site_coords.csv is the file I use as input in GIS 

