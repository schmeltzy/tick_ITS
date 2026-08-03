## make map of sites in VT

library(sf)
library(ggplot2)
library(cowplot)
library(tigris)
library(ggrepel)
library(ggspatial)  # Required for scale bar and north arrow

# 1. Get Vermont state and town boundaries
options(tigris_use_cache = TRUE)
vt_state <- states(cb = TRUE) |> subset(STUSPS == "VT")
vt_towns <- county_subdivisions(state = "VT", cb = TRUE)

# Transform to a common projected CRS (NAD83 / Vermont zone)
vt_state_prj <- st_transform(vt_state, 32145)
vt_towns_prj <- st_transform(vt_towns, 32145)
manchester_prj <- subset(vt_towns_prj, NAME == "Manchester")

# NEW: Download roads for Bennington County (where Manchester is located)
# The roads data is then reprojected and clipped strictly to Manchester's border
bennington_roads <- roads(state = "VT", county = "Bennington")
roads_prj <- st_transform(bennington_roads, 32145)
manchester_roads <- st_intersection(roads_prj, manchester_prj)

# NEW: Download linear water features (rivers/streams) for Bennington County & clip
bennington_water_lines <- linear_water(state = "VT", county = "Bennington")
water_lines_prj <- st_transform(bennington_water_lines, 32145)
manchester_rivers <- st_intersection(water_lines_prj, manchester_prj)

# NEW (Optional): Download area water features (lakes/ponds) for Bennington County & clip
bennington_water_areas <- area_water(state = "VT", county = "Bennington")
water_areas_prj <- st_transform(bennington_water_areas, 32145)
manchester_lakes <- st_intersection(water_areas_prj, manchester_prj)


# 2. Define specific Lat/Long sites in Manchester
sites_df <- data.frame(
  name = c("managed", "unmanaged", "burned","unburned","Village of Manchester"),
  longitude = c(-73.07550, -73.076264, -73.077647, -73.075342, -73.07111),
  latitude = c(43.14423, 43.142600, 43.131661, 43.131564, 43.16111)
)

sites_sf <- st_as_sf(sites_df, coords = c("longitude", "latitude"), crs = 4326)
sites_prj <- st_transform(sites_sf, 32145)

sites_coords <- st_coordinates(sites_prj)
sites_prj$X <- sites_coords[,1]
sites_prj$Y <- sites_coords[,2]

# 3. Create Main Vermont Map (with Scale Bar & North Arrow)
p_main_extended <- ggplot() +
  geom_sf(data = vt_state_prj, fill = "ivory", color = "black") +
  geom_sf(data = manchester_prj, fill = "lightblue", color = "red", linewidth = 0.8) +
  
  # Add text label for the town of Manchester
  geom_sf_text(
    data = manchester_prj, 
    aes(label = NAME), 
    fontface = "bold",
    size = 4, 
    color = "black",
    nudge_x = 2000,  # Shifts label to the right (units are in meters for CRS 32145)
    nudge_y = -11000   # Shifts label slightly down
  ) +
  
  # Add North Arrow to main map (bottom-left)
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.2, "in"), pad_y = unit(0.6, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Add Scale Bar to main map (bottom-left)
  annotation_scale(location = "br", width_hint = 0.2) +
  theme_minimal() +


# FORCE AXIS TICK LABELS TO BE DECIMAL DEGREES (LAT/LONG)
coord_sf(datum = st_crs(4326)) +
  
  theme_minimal() +
  
  # RENAME THE AXIS LABELS HERE
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme(plot.margin = margin(t = 5, r = 240, b = 5, l = 5, unit = "pt"))

# 4. Create Inset Map (with independent small-scale bar)
p_inset <- ggplot() +
  geom_sf(data = manchester_prj, fill = "cornsilk", color = "grey50") +
  
  # NEW: Draw the underlying road street paths inside the town boundary
  geom_sf(data = manchester_roads, color = "grey75", linewidth = 0.4) +
  
  # NEW: Draw linear river channels (Soft blue lines)
  geom_sf(data = manchester_rivers, color = "skyblue3", linewidth = 0.3) +
  
  # NEW: Draw area water features like ponds/lakes (Soft blue fill shapes)
  geom_sf(data = manchester_lakes, fill = "skyblue1", color = "skyblue3", linewidth = 0.2) +
  
  geom_sf(data = sites_prj, aes(color = name), size = 2, show.legend = FALSE) +
  geom_label_repel(
    data = sites_prj,
    aes(x = X, y = Y, label = name, color = name),
    box.padding = 0.5, point.padding = 0.3,
    segment.color = "grey30", segment.linewidth = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(values = c("managed" = "#009E73", "unmanaged" = "#CC79A7", "burned" = "#D55E00", "unburned" = "#0072B2", "Village of Manchester" = "black")) +
  # Add a highly localized scale bar inside Manchester town bounds
  annotation_scale(location = "tr", width_hint = 0.10) + 
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1.2),
    plot.margin = margin(5, 5, 5, 5)
  )

# 5. Define Inset Map Positioning
inset_x <- 0.55
inset_y <- 0.05
inset_w <- 0.43
inset_h <- 0.90


# =========================================================================
# 6. MANUAL DASHED LINE ADJUSTMENTS (Change these coordinates manually!)
# =========================================================================
# Values represent percentages of the whole image workspace canvas (0.0 to 1.0)

# Where the lines start on the main Vermont map (pointing at Manchester)
manc_x <- 0.205
manc_y <- 0.335

# Top dashed line target (Points to top-left corner of the inset frame)
top_target_x <- 0.550  
top_target_y <- 0.751  

# Bottom dashed line target (Points to bottom-left corner of the inset frame)
bottom_target_x <- 0.555  
bottom_target_y <- 0.24  
# =========================================================================


# 7. Assemble Canvas and Plot Lines
final_map_floating <- ggdraw(p_main_extended) +
  # Render Inset Map
  draw_plot(p_inset, x = inset_x, y = inset_y, width = inset_w, height = inset_h) +
  
  # Top Connection Line
  draw_line(
    x = c(manc_x, top_target_x),
    y = c(manc_y, top_target_y),
    color = "black", linewidth = 0.8, linetype = "dashed"
  ) +
  
  # Bottom Connection Line
  draw_line(
    x = c(manc_x, bottom_target_x),
    y = c(manc_y, bottom_target_y),
    color = "black", linewidth = 0.8, linetype = "dashed"
  )

# Preview map
print(final_map_floating)

# Save map with locked down dimensions so your manual coordinates don't warp
#ggsave("tick_sites_vermont.png", final_map_floating, width = 11, height = 8, dpi = 300)

