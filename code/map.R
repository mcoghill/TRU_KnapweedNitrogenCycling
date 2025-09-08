library(tidyverse)
library(bcdata)
library(bcmaps)

# Download BC polygon
bc <- bc_bound() %>% 
  st_union() %>% 
  st_as_sf() %>% 
  mutate(NAME = "British Columbia")

# Download city point
kam <- bc_cities() %>% 
  filter(NAME == "Kamloops") %>% 
  select(NAME) %>% 
  mutate(NAME = "Kamloops, BC")

# Download road data for around the sampling location
pt <- data.frame(x = -120.446611, y = 50.814083, id = "Collection site") %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3153)

# pt_buf <- st_buffer(pt, 1000)
# roads <- bcdc_query_geodata("bb060417-b6e6-4548-b837-f9060d94743e", crs = 3153) %>% 
#   filter(INTERSECTS(pt_buf), ROAD_NAME_FULL == "Lac Du Bois Rd") %>% 
#   collect()

# Write spatial data
write_sf(bc, "map/bc.gpkg")
write_sf(kam, "map/kamloops.gpkg")
write_sf(pt, "map/sample_site.gpkg")
# write_sf(roads, "map/roads.gpkg")
