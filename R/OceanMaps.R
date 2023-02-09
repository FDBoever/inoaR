#install.packages('terra')
#install.packages('terra', repos='https://rspatial.r-universe.dev')
#install.packages(
#  "ggOceanMapsData",
#  repos = c("https://mikkovihtakari.github.io/drat",
#            "https://cloud.r-project.org")
#)


library(ggOceanMapsData)
library(ggOceanMaps)
basemap(limits = 60) # A synonym: basemap(60)

library(ggmap)
wtr.smpls <- read.delim('~/DATA/Predator_lab/water_samples.txt', header=TRUE)
CCAPcolors=c('#59A86C')

worldmap.d = map_data('world')

worldmap <- ggplot() +
  geom_polygon(data = worldmap.d, aes(x=long,
                                      y = lat,
                                      group = group),
               size=.4,color='white',
               fill="grey",
               alpha=0.3) +
  geom_point(data=wtr.smpls, aes(x=Longitude,
                                 y=Latitude),
             color=CCAPcolors[1]) +
  theme_void() + coord_fixed()+xlim(-50,30)+ylim(0,90)

ggsave('~/DATA/Predator_lab/water_world_map.pdf',worldmap)
