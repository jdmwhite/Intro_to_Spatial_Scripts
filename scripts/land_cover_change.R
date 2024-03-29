#### Install packages ----
# install.packages('ggalluvial') 
# install.packages('patchwork')
# install.packages('mapview')

#### Load libraries ----
library(sf)
library(terra)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(mapview)

#### Load data ----
lc1990 <- terra::rast("data/land_cover_change/SANLC_1990_COJ_extent.tif")
lc2020 <- rast("data/land_cover_change/SANLC_2020_COJ_extent.tif")
coj <- vect("data/land_cover_change/COJ_boundary.shp")

# Plot our raw data
par(mfrow = c(1,3))
plot(lc1990)
plot(lc2020)
plot(coj)

#### Check projections ----
print(paste0("Do the layers have the same CRS? ", crs(lc1990)==crs(lc2020)))

lc1990 <- project(lc1990, lc2020)
print(paste0("Do the layers have the same CRS? ", crs(lc1990)==crs(lc2020)))

#Reproject the coj vector to match the land cover rasters
coj <- project(coj, lc2020)

# Plot maps to view our data
par(mfrow = c(1,1))
plot(lc1990)
plot(coj, add = TRUE)

#### Crop & Mask ----
# First crop and then mask the land cover data to our area of interest
lc1990_aoi <- mask(crop(lc1990, coj), coj)
lc2020_aoi <- mask(crop(lc2020, coj), coj)
print(paste0("Do the layers have the same extent? ", ext(lc1990_aoi) == ext(lc2020_aoi)))

# Check to see if this worked
plot(lc2020_aoi)
plot(coj, add = TRUE)

#### Reclassify ----
# We have more than 70 classes for each of the land cover datasets. At this point we want to reclassify these classes to simpler classes
# Let's merge all of our classes into 4 broad categories:
# Key
# 1 = Water
# 2 = Agriculture
# 3 = Artificial surfaces
# 4 = Vegetation

# Create a reclassification matrix for 1990
m1990 <- rbind(c(0, 3, 1),
          c(36, 38, 1),
          c(9, 31, 2),
          c(34, 36, 3),
          c(38, 51, 3),
          c(52, 56, 3),
          c(60, 72, 3),
          c(3, 9, 4),
          c(31, 34, 4),
          c(51, 52, 4),
          c(56, 60, 4))

# Reclassify using the terra::classify function
lc1990_rcl <- classify(lc1990_aoi, m1990)
plot(lc1990_rcl)

# Create a reclassification matrix for 2020
m2020 <- rbind(c(13, 24, 1),
               c(31, 46, 2),
               c(24, 31, 3),
               c(46, 60, 3),
               c(64, 73, 3),
               c(0, 13, 4),
               c(13, 19, 4),
               c(19, 25, 4),
               c(60, 64, 4))

# Reclassify 2020
lc2020_rcl <- classify(lc2020_aoi, m2020)
plot(lc2020_rcl)

#### Crosstab ----
# We now want to calculate the pairwise changes between the 1990 and 2020 land cover data.

# stack the land cover
landcover_stack <- c(lc2020_rcl, lc1990_rcl)
# Run a change analysis using the terra::crosstab function. long = T returns a data frame instead of a table.
lc_changes <- crosstab(landcover_stack, long = TRUE)
head(lc_changes)

lc2020_rcl
# tidy up this output by changing the raster names, calculating the area of each class and assigning the full names back to the numbers
lc_changes %>% 
  rename(rcls_2020 = SANLC_2020_COJ_extent,
        rcls_1990 = SANLC_1990_COJ_extent) %>%
  mutate(area = Freq*400/1e6,
  lc1990 = case_when(
    rcls_1990 == 1 ~ 'Water',
    rcls_1990 == 2 ~ 'Agriculture',
    rcls_1990 == 3 ~ 'Artificial',
    rcls_1990 == 4 ~ 'Vegetation'),
  lc2020 = case_when(
    rcls_2020 == 1 ~ 'Water',
    rcls_2020 == 2 ~ 'Agriculture',
    rcls_2020 == 3 ~ 'Artificial',
    rcls_2020 == 4 ~ 'Vegetation')
) %>% select(lc1990, lc2020, area) -> lc_changes_labelled

# Calculate the area of COJ
paste(terra::expanse(coj)/1000000, 'km2')
# Calculate the area of all changes
paste(sum(lc_changes_labelled$area), 'km2')

#### Visualisation ----
# First let's run a sankey/alluvium plot
alluv_plot <- ggplot(lc_changes_labelled, aes(axis1 = lc1990, axis2 = lc2020, y = area)) +
  geom_alluvium(aes(fill = lc1990)) +
  scale_fill_manual(values = c('#7E6148B2','#F39B7FB2','#00A087B2','#4DBBD5B2'), guide = 'none') +
  geom_stratum(fill = c('#4DBBD5B2','#00A087B2','#F39B7FB2','#7E6148B2','#4DBBD5B2','#00A087B2','#F39B7FB2','#7E6148B2'), col = NA, alpha = 0.8) +
  geom_text(stat = 'stratum', aes(label = paste(after_stat(stratum),'\n',round(after_stat(prop)*100,1))), size = 2.5) +
  scale_x_continuous(breaks = c(1, 2), labels = c('1990','2020'), position = 'top') +
  theme_void() +
  theme(axis.text.x = element_text())

alluv_plot

#### Plot maps
# downsample the rasters using aggregate & convert to a data frame for plotting with ggplot2
lc1990_a <- aggregate(lc1990_rcl, 5, fun = 'modal', by = 'SANLC_1990_COJ_extent')
lc1990_df <- as.data.frame(lc1990_a, xy = TRUE)
names(lc1990_df)[3] <- 'land_cover'

lc2020_a <- aggregate(lc2020_rcl, 5, fun = 'modal', by = 'SANLC_2020_COJ_extent')
lc2020_df <- as.data.frame(lc2020_a, xy = TRUE)
names(lc2020_df)[3] <- 'land_cover'

# plot the reclassified land cover for 1990
lc1990_plot <- ggplot(lc1990_df) +
  geom_tile(aes(x = x, y = y, fill = as.factor(land_cover), col = as.factor(land_cover))) +
  scale_fill_manual(values = c('#4DBBD5B2','#7E6148B2','#F39B7FB2','#00A087B2'),labels = c('Water', 'Agriculture', 'Artificial', 'Vegetation'), guide = 'none') +
  scale_colour_manual(values = c('#4DBBD5B2','#7E6148B2','#F39B7FB2','#00A087B2'), guide = 'none') +
  geom_sf(data = st_as_sf(coj), fill = NA, col = 'black', lwd = 0.2) +
  labs(title = '1990', fill = 'Land Cover') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) 
  
# plot the reclassified land cover for 2020
lc2020_plot <- ggplot(lc2020_df) +
  geom_tile(aes(x = x, y = y, fill = as.factor(land_cover))) +
  geom_tile(aes(x = x, y = y, fill = as.factor(land_cover), col = as.factor(land_cover))) +
  scale_fill_manual(values = c('#4DBBD5B2','#7E6148B2','#F39B7FB2','#00A087B2'),labels = c('Water', 'Agriculture', 'Artificial', 'Vegetation'), guide = 'none') +
  scale_colour_manual(values = c('#4DBBD5B2','#7E6148B2','#F39B7FB2','#00A087B2'), guide = 'none') +
  geom_sf(data = st_as_sf(coj), fill = NA, col = 'black', lwd = 0.2) +
  labs(title = '2020', fill = 'Land Cover') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Combine the plots together using patchwork and add on plot labels
lc_plots <- lc1990_plot + lc2020_plot + alluv_plot & plot_annotation(tag_levels = 'a', tag_suffix = ')')
lc_plots

# Save the output
ggsave('output/figs/land_cover_change/land_cover_plots.png', lc_plots,
       width = 180, height = 100, units = c('mm'), dpi = 'retina')

#### Interactive map ----
cls <- c('Water', 'Agriculture', 'Artificial', 'Vegetation')
lc1990_raster <- as.factor(raster::raster(lc1990_a))
lc1990_raster[] = factor(cls[lc1990_raster[]])
lc2020_raster <- as.factor(raster::raster(lc2020_a))
lc2020_raster[] = factor(cls[lc2020_raster[]])

m <- mapview(lc1990_raster, na.color = NA, layer.name = 'Land Cover 1990', alpha = 1) +
  mapview(lc2020_raster, na.color = NA, layer.name = 'Land Cover 2020', alpha = 1, legend = FALSE)
m

mapshot(m, "output/figs/land_cover_change/interactive_map.html")
