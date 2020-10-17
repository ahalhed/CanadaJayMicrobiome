# addition plots scripts
# not associated with a give hypothesis


# map of smapling locations
# combine location data with metadata
comb <- merge(XY_gj, meta_sub, by=0) 
# plot the locations (test run)
ggplot(comb, aes(y = Latitude, x = Longitude,
                 shape=as.character(CollectionYear))) +
  geom_count() +  #facet_grid(CollectionYear~.) +
  labs(y = "Latitude", x = "Longitude", shape = "Collection Year", size = "Number of Samples")


# these locations account for ALL sample locations (2016-2020)
map_gj <- get_map(
  location = c(left = -79.5, bottom = 45.2, right = -77.9, top = 46.2),
  source = "osm",
  force = TRUE) # adding force = TRUE to get_map to force a re-rendering of the map
# Let's look at our map
ggmap(map_gj)


# Let's just throw everything onto the map to see what we've got
ggmap(map_gj) + 
  # switch count to jitter to get around being right on top of one another
  geom_count(data = comb, 
             aes(y = Latitude, x = Longitude, shape = as.character(CollectionYear))) + 
  #facet_grid(CollectionYear~.) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  labs(shape = "Collection Year", size = "Number of Samples",
       title = "Map of Canada Jay Sampling Locations",
       subtitle = "Algonquin Park, Ontario (2016-2020)")