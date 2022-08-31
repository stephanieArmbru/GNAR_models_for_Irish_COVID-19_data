### MODELLING COVID-19 INCIDENCE IN IRISH COUNTIES IN GNAR MODELS WITH 
### DIFFERENTLY CONSTRUCTED NETWORKS 


# Install packages --------------------------------------------------------

# install.packages("geosphere")
# install.packages("spdep")
# install.packages("raster")
# install.packages("gganimate")
# install.packages("gifski")
# install.packages("flexclust")
# install.packages("rgeos")
# install.packages("dbscan")
# install.packages("latex2exp")
# devtools::install_github("JosephCrispell/basicPlotteR")
# install.packages("qpcR")

# set seed to guarantee reproducibility 
set.seed(1234)

# Load libraries ----------------------------------------------------------
library(readr)
library(igraph)
library(GNAR)
library(MASS) # for box cox 
library(tidyverse)
library(magrittr) # for pipping 
library(xtable) # for tables 
library(geosphere) # for long / lat distance calculation

# for MAPS: 
library(sp)
library(raster) 
library(leaflet)
library(mapview) # for saving

# for neighbourhood construction 
library(spdep) 
# to create neighbourhood data frame
library(expp) 

# for sphere of influence network 
library(rgeos)
library(dbscan)

# for visualising the network
library(sf)
library(lwgeom)

# to transform to shapefile 
library(raster)

# for pairwise distances
library(flexclust)

# for latex expressions in graph titles 
library(latex2exp)

# for non-overlapping labels in base plot 
library(basicPlotteR)


# Load functions 
source(file = "functions.R")

# Turn off warnings
options(warn = -1)

# Set ggplot theme 
theme_set(theme_bw(base_size = 16))

# Load pre-processed data -------------------------------------------------
COVID_weekly_data <- read_csv(file = "data/COVID/ireland_covid_weekly.csv", 
                       show_col_types = FALSE)

# correct format for GNAR models
covid_cases <-  COVID_weekly_data %>% 
  dplyr::select(CountyName,
                yw, 
                weeklyCases) %>% 
  spread(CountyName, 
         weeklyCases) %>% 
  column_to_rownames(var = "yw") %>% 
  as.matrix()

# transform to data frame with time column for plotting 
covid_cases_df <- covid_cases %>% 
  as.data.frame() %>% 
  mutate(time = as.Date(rownames(covid_cases))) 


# Parameters --------------------------------------------------------------
# for COVID-19 reproduction number 
# transmission rate
beta_covid <- 0.89

# recovery rate 
gamma_covid <- 1/6.0

# Economic hubs
hubs <- c("Dublin", 
          "Cork",
          "Limerick", 
          "Galway", 
          "Waterford")

# create vector after which tables for latex output are ordered 
network_ordering <- c("delaunay", 
                      "Gabriel",
                      "SOI", 
                      "Relative", 
                      "Queen", 
                      "Eco. hub", 
                      "")

# Load maps data ----------------------------------------------------------
# download Ireland data level 1 from GADM database
ireland <- raster::getData('GADM',
                   country='IRL',
                   level = 1)

# comparison of county spelling in maps data and COVID data 
setdiff(COVID_weekly_data$CountyName %>% unique(), ireland$NAME_1)
setdiff(ireland$NAME_1, COVID_weekly_data$CountyName %>% unique())
# different writing: Laois and Laoighis (the former is the newer spelling)

# correct spelling in ireland data
ireland$NAME_1[ireland$NAME_1 == "Laoighis"] <- "Laois"

# save Spatial Polygon Data Frame as shapefile 
# shapefile(x = ireland, 
#           file = "data/COVID/ireland_shapefile.shp")


# read Ireland shapefile 
ireland_shp <- st_read("data/COVID/ireland_shapefile.shp")

# construct centroid for each county 
coord <- ireland_shp %>%  
  st_geometry() %>% 
  st_centroid() %>% 
  st_coordinates()

# assign county names to their centroid coordinates 
rownames(coord) <- ireland_shp$NAME_1


# read data frame including county towns and their coordinates 
county_towns <- read_delim("data/COVID/county_towns_ireland.csv", 
                           delim = ";", 
                           escape_double = FALSE, 
                           trim_ws = TRUE, 
                           show_col_types = FALSE, 
                           locale = locale(decimal_mark = ","))

# match county towns and COVID incidence 
c1 <- county_towns %>% dplyr::pull(admin_name) %>% unique()
c2 <- COVID_weekly_data$CountyName %>% unique()

# check if all counties are included 
diff_c <- setdiff(c1, c2)
same_c <- intersect(c1, c2)

# filter out all relevant county towns 
county_towns_lim <- county_towns %>% 
  filter(admin_name %in% same_c)

# investigate counties with multiple county towns 
double_county <- county_towns_lim %>% 
  dplyr::pull(admin_name) %>% 
  duplicated()

county_towns_lim %>% 
  filter(admin_name %in% county_towns_lim[double_county, ]$admin_name)

# manually set county towns for counties with multiple county towns 
county_towns_lim_unique <- county_towns_lim %>% 
  filter(city != "Drogheda", 
         city != "Nenagh")

# check if right number of counties included 
county_towns_lim_unique %>% nrow()

# construct matrix with coordinates 
coord_central_towns <- county_towns_lim_unique %>% 
  dplyr::select(lat, lng)

# copy data frame for two different rownames 
coord_central_towns_counties <- coord_central_towns

# assign county towns as rownames 
rownames(coord_central_towns) <- county_towns_lim_unique %>% 
  pull(city)
# assign county names as rownames 
rownames(coord_central_towns_counties) <- county_towns_lim_unique %>% 
  pull(admin_name)

# transform into matrix
coords_ct <- as.matrix(coord_central_towns)[, c(2,1)]
coords_ct_counties <- as.matrix(coord_central_towns_counties)[, c(2,1)]


# Classification ----------------------------------------------------------
# classification of counties into urban and rural according to Census 2016 data
urban <- c("Dublin", "Louth", "Meath", 
           "Kildare", "Wicklow", "Waterford", 
           "Limerick", "Galway", "Cork", "Longford")
rural <- c("Wexford", "Kilkenny", "Carlow", 
           "Laois", "Kerry", "Tipperary", 
           "Offaly", "Westmeath", "Cavan", 
           "Monaghan", "Leitrim", "Mayo", 
           "Sligo", "Donegal", "Clare",
           "Roscommon")

# check that no counties are double classified, and all counties are included 
(rural %in% urban) %>% any()
c(rural, urban) %>% length()

colnames(covid_cases) %in% c(rural, urban) %>% all()

# create county name vector in the correct order  
counties <- covid_cases %>% colnames()

# create factor vector indicating urbanisation status for each county
urbanisation_factor <- c()
for (name in counties) {
  if (name %in% rural) {
    urbanisation_factor <- c(urbanisation_factor, "rural")
  } 
  if (name %in% urban) {
    urbanisation_factor <- c(urbanisation_factor, "urban")
  }
}

# reformat as factor 
urbanisation_factor <- urbanisation_factor %>% as.factor()


# create mixture coordinate matrix 
# rural: centroid coordinates as node representation
# urbal: county town coordinates as node representation 

# guarantee same ordering 
coord_ord <- coord[match(counties, rownames(coord)), ]
coords_ct_counties_ord <- coords_ct_counties[match(counties, 
                                                   rownames(coords_ct_counties)), 
                                             ]

# filter coordinates for each county
coord_rural <- coord_ord[urbanisation_factor == "rural", ]
coord_urban <- coords_ct_counties_ord[urbanisation_factor == "urban", ]

# create matrix with coordinates for every county 
coord_urbanisation_prelim <- rbind(coord_rural, 
                             coord_urban)
# guarantee same ordering 
coord_urbanisation <- coord_urbanisation_prelim[match(counties, 
                                                      rownames(coord_urbanisation_prelim)), 
                                                ]

# check if order is correct 
all(rownames(coord_urbanisation) == counties)

# Great Circle distance ---------------------------------------------------
# compute pairwise Great Circle distances between all counties 
dist_urbanisation <- circle_distance(coord_urbanisation)

# construct data frame including population density in Tsd. for each county
population_weight <- COVID_weekly_data %>% 
  dplyr::select(CountyName, 
                PopulationCensus16) %>% 
  unique() %>% 
  mutate(weight = PopulationCensus16 / 1000) %>% 
  dplyr::select(-PopulationCensus16)

# General maps ------------------------------------------------------------
# plot maps of counties
pal <- colorFactor("Reds", ireland$NAME_1)

leaflet(ireland) %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(NAME_1),
              highlightOptions = highlightOptions(color = "white", weight =2,
                                                  bringToFront = TRUE),
              label = ~NAME_1)

# visualise centres according to urbanisation state 
plot(st_geometry(ireland_shp),
     border="grey",
     col = c("#78BE21", "#516FC4")[urbanisation_factor],
     legend = c("rural", "urban"))
points(coord_urbanisation[, 1],
       coord_urbanisation[, 2])
addTextLabels(xCoords = coord_urbanisation[, 1],
              yCoords = coord_urbanisation[, 2],
              labels=rownames(coord_urbanisation),
              cex.label = 0.9,
              col.label = "black",
              col.line = "black")


# COVID maps --------------------------------------------------------------
# plot 1-lag COVID-19 ID for weeks spread across the duration of the pandemic 
ggplot(COVID_weekly_data %>% 
         filter(yw == "2020-03-01" |  
                  yw == "2020-09-27" |  
                  yw == "2021-04-11" | 
                  yw == "2021-12-05" |
                  yw == "2022-06-12"), 
       aes(x = CountyName,
           y = weeklyCases, 
           color = as.factor(yw), 
           group = yw)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1), 
        legend.title = element_blank(), 
        legend.position = "bottom") +
  ylab("COVID-19 ID") +
  xlab("County")
ggsave("plots/dataVisualisation/covid_compare_across_time.pdf", 
       width = 23, 
       height = 14, 
       unit = "cm")
# counties with highest weekly cases change during the pandemic 
# not one hub from which COVID-19 was spread identifiable 

# Network construction ----------------------------------------------------
# 1. Queen's contiguity: connection to all neighbouring counties 
nb_list_queen <- poly2nb(ireland, 
                         queen = TRUE, 
                         row.names = ireland$NAME_1)

# extract neighbourhood (nb) list and convert to igraph object 
covid_net_queen_igraph <- neighborsDataFrame(nb = nb_list_queen) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() # simplify to obtain no self-loops and no multiple edges

# create GNAR object 
covid_net_queen_gnar <- covid_net_queen_igraph %>% igraphtoGNAR()

# create ordered county index data frame
county_index_queen <- data.frame("CountyName" = covid_net_queen_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))

# visualise network 
plot(st_geometry(ireland_shp),
     border="grey")
plot(nb_list_queen,
     coord_urbanisation,
     pch = 19, cex = 0.6,
     add = TRUE)
text(coord_urbanisation[, 1],
     coord_urbanisation[, 2],
     labels = rownames(coord_urbanisation),
     cex = 0.8,
     font = 2,
     pos = 1)

# 2. Rook's contiguity
nb_list_rook <- poly2nb(ireland, 
                        queen = FALSE, 
                        row.names = ireland$NAME_1)

# analyse differences between Queen's and Rook's network
diffnb(x = nb_list_queen, 
       y = nb_list_rook)
# produces same links due to non-grid like structure 

# 3. Economic hubs
# network with edges to the direct neighbors and nearest economic hub
covid_net_eco_hubs <- create_eco_hub_igraph(dist_urbanisation, 
                                            coord_urbanisation)
# create igraph object
covid_net_eco_hubs_igraph <- covid_net_eco_hubs$igraph_net
# create GNAR object
covid_net_eco_hubs_gnar <- covid_net_eco_hubs_igraph %>% igraphtoGNAR()

# extract nb list from igraph object 
nb_eco_hubs <- igraph2nb(gr = covid_net_eco_hubs_igraph)

ord_coord_eco_hubs <- covid_net_eco_hubs$ordered_coord
coord_hubs <- covid_net_eco_hubs$ordered_hubs

# create county index data frame 
county_index_eco_hubs <- data.frame("CountyName" = ord_coord_eco_hubs %>%
                                      rownames(), 
                                   "index" = seq(1, 26))

# visualise network 
plot(st_geometry(ireland_shp),
     border = "grey")
plot(nb_eco_hubs,
     ord_coord_eco_hub,
     add = TRUE,
     pch = 19, cex = 0.6)
text(coord_hubs$X,
     coord_hubs$Y,
     labels = rownames(coord_hubs),
     cex = 0.8, font = 2, pos = 1,
     col = "#516FC4")


# 3. Railway-based network 
# Newbridge is biggest town in Kildare, Naas is county town, 
train_connections <- matrix(c("Dublin", "Dún Dealgan", 
                          "Dublin", "Trim", 
                          "Dublin", "Naas", 
                          "Dublin", "Wicklow", 
                          "Cork", "Tralee", 
                          "Cork", "Limerick", 
                          "Cork", "Clonmel", 
                          "Galway", "Ennis", 
                          "Galway", "Tullamore", 
                          "Galway", "Ros Comáin", 
                          "Limerick", "Ennis", 
                          "Limerick", "Clonmel", 
                          "Limerick", "Port Laoise", 
                          "Limerick",   "Cork", 
                          "Waterford", "Kilkenny", 
                          "Waterford", "Clonmel", 
                          "Dún Dealgan", "Dublin", # Dundalk 
                          "Tralee", "Cork", 
                          "Carlow", "Naas", 
                          "Carlow", "Kilkenny", 
                          "Carlow", "Tullamore", 
                          "Carlow", "Port Laoise",
                          "Ennis", "Galway", 
                          "Ennis",    "Limerick", 
                          "Kilkenny", "Carlow", 
                          "Kilkenny", "Waterford", 
                          "Naas", "Dublin", 
                          "Naas", "Carlow", 
                          "Naas",   "Port Laoise", 
                          "Naas",   "Tullamore", # Newbridge
                          "Sligo", c("Carrick on Shannon"), 
                          "Ros Comáin", "Castlebar", 
                          "Ros Comáin",  "Tullamore", 
                          "Ros Comáin", "Galway", # Roscommon
                          "Mullingar", "Trim", 
                          "Mullingar", "Longford", 
                          "Wicklow", "Dublin", 
                          "Wicklow", "Wexford", 
                          "Clonmel", "Waterford", 
                          "Clonmel", "Limerick", 
                          "Clonmel", "Cork", 
                          "Clonmel", "Tralee", 
                          "Clonmel", "Port Laoise",
                          "Wexford", "Wicklow", 
                          "Longford", "Carrick on Shannon", 
                          "Longford", "Mullingar", 
                          "Trim", "Mullingar", 
                          "Trim", "Dublin", # Enfield
                          "Carrick on Shannon", "Sligo", 
                          "Carrick on Shannon", "Longford", 
                          "Tullamore", "Ros Comáin", 
                          "Tullamore", "Port Laoise", 
                          "Tullamore", "Naas", 
                          "Tullamore", "Carlow", 
                          "Tullamore", "Galway", 
                          "Port Laoise", "Tullamore", 
                          "Port Laoise", "Naas",
                          "Port Laoise", "Carlow", 
                          "Port Laoise", "Limerick", 
                          "Port Laoise", "Cork", 
                          "Port Laoise", "Tralee", 
                          "Castlebar", "Ros Comáin", 
                          "Monaghan", "Dún Dealgan", 
                          "An Cabhán", "Longford", 
                          "Lifford", "Sligo"
                          ), ncol = 2, byrow = TRUE) %>% 
  as.data.frame()
colnames(train_connections) <- c("city", "city2")


no_train_counties <- c("Cavan", 
                       "Donegal", 
                       "Monaghan")
# "An Cabhán" = NA, # Cavan, not reachable by train 
# "Monaghan" = NA, # not reachable by train 
# "Lifford" = NA # Donegal, not reachable by train


# substitute county towns with county index for nb list generation
county_to_town <- county_towns_lim_unique %>% 
  dplyr::select("city", "admin_name") %>%
  mutate("name_index" = seq(1, 26))

# save order to order coordinates adequately for plot 
order_names <- county_to_town %>% 
  arrange(name_index) %>% 
  pull(admin_name)
  
town_to_index <- county_to_town %>% 
  dplyr::select(-"admin_name")

# substitute town names with corresponding indices since function edgelist() 
# requires numeric vertices
first_step <- left_join(train_connections, town_to_index, by = "city")
colnames(first_step) <- c("sub", "city", "county1")
second_step <- left_join(first_step, town_to_index, by = "city")
colnames(second_step) <- c("sub", "sub2", "county1", "county2")

# create edge list
train_connections_m <- second_step %>% 
  dplyr::select("county1",
                "county2") %>% 
  as.matrix(ncol = 2, byrow = TRUE)

# create igraph object from edge list
covid_net_train_igraph <- graph_from_edgelist(train_connections_m, 
                                              directed = FALSE)
# create GNAR object 
covid_net_train_gnar <- covid_net_train_igraph %>% igraphtoGNAR()

# create nb list object 
nb_list_train <- igraph2nb(gr = covid_net_train_igraph)

# reorder so the ordering of the nb list and the coordinates match
coord_urb_ord <- coord_urbanisation[match(order_names, 
                                         rownames(coord_urbanisation)), ] %>% 
  as.data.frame()

# create ordered county index data frame 
county_index_train <- data.frame("CountyName" = coord_urb_ord %>% rownames(), 
                                 "index" = seq(1, 26))

# visualise network 
plot(st_geometry(ireland_shp),
     border = "grey")
plot(nb_list_train,
     coord_urb_ord,
     add = TRUE,
     pch = 19, cex = 0.6)
text(coord_urb_ord$X,
     coord_urb_ord$Y,
     labels = rownames(coord_urb_ord),
     cex = 0.8, font = 2, pos = 1)


# 4. Delaunay triangulation 
nb_list_delaunay <- tri2nb(coords = coord_urbanisation, 
                           row.names = coord_urbanisation %>% 
                             rownames())

# create igraph object
covid_net_delaunay_igraph <- neighborsDataFrame(nb = nb_list_delaunay) %>% 
    graph_from_data_frame(directed = FALSE) %>%
    igraph::simplify()

# create GNAR object 
covid_net_delaunay_gnar <- covid_net_delaunay_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_delaunay <- data.frame("CountyName" = covid_net_delaunay_igraph %>%
                                      V() %>% 
                                      names(), 
                                    "index" = seq(1, 26))
# visualise network 
plot(st_geometry(ireland_shp),
     border="grey")
plot(nb_list_delaunay,
     coord_urbanisation,
     pch = 19, cex = 0.6,
     add=TRUE)
text(coord_urbanisation[, 1],
     coord_urbanisation[, 2],
     labels = rownames(coord_urbanisation),
     cex = 0.8, font = 2, pos = 1)


# 5. Gabriel neighbourhoods 
nb_list_gabriel <- gabrielneigh(coords = coord_urbanisation) %>% 
  graph2nb(row.names = coord_urbanisation %>% 
             rownames(), sym = TRUE)
# symmetry required for every region to have edges 

# create igraph object
covid_net_gabriel_igraph <- neighborsDataFrame(nb = nb_list_gabriel) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object
covid_net_gabriel_gnar <- covid_net_gabriel_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_gabriel <- data.frame("CountyName" = covid_net_gabriel_igraph %>%
                                     V() %>% 
                                     names(), 
                                   "index" = seq(1, 26))

# visualise network
plot(st_geometry(ireland_shp),
     border="grey")
plot(nb_list_gabriel,
     coord_urbanisation,
     pch = 19, cex = 0.6,
     add=TRUE)
text(coord_urbanisation[, 1],
     coord_urbanisation[, 2],
     labels = rownames(coord_urbanisation),
     cex = 0.8, font = 2, pos = 1)

# 6. Relative neighbourhood
nb_list_relative <- relativeneigh(coords = coord_urbanisation) %>% 
  graph2nb(row.names = coord_urbanisation %>% 
             rownames(), 
           sym = TRUE)

# create igraph object 
covid_net_relative_igraph <- neighborsDataFrame(nb = nb_list_relative) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
covid_net_relative_gnar <- covid_net_relative_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_relative <- data.frame("CountyName" = covid_net_relative_igraph %>%
                                      V() %>% 
                                      names(), 
                                    "index" = seq(1, 26))

# visualise network 
plot(st_geometry(ireland_shp),
     border="grey")
plot(nb_list_relative,
     coord_urbanisation,
     pch = 19, cex = 0.6,
     add=TRUE)
text(coord_urbanisation[, 1],
     coord_urbanisation[, 2],
     labels = rownames(coord_urbanisation),
     cex = 0.8, font = 2, pos = 1)

# 7. Sphere of Influence (SOI) neighbourhood 
nb_list_soi <- soi.graph(nb_list_delaunay, 
                         coord_urbanisation) %>% 
  graph2nb(row.names = coord_urbanisation %>% 
             rownames())

# create igraph object 
covid_net_soi_igraph <- neighborsDataFrame(nb = nb_list_soi) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
covid_net_soi_gnar <- covid_net_soi_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame 
county_index_soi <- data.frame("CountyName" = covid_net_soi_igraph %>%
                                 V() %>% 
                                 names(), 
                               "index" = seq(1, 26))

# visualise network 
plot(st_geometry(ireland_shp),
     border="grey")
plot(nb_list_soi,
     coord_urbanisation,
     pch = 19, cex = 0.6,
     add=TRUE)
text(coord_urbanisation[, 1],
     coord_urbanisation[, 2],
     labels = rownames(coord_urbanisation),
     cex = 0.8, font = 2, pos = 1)


# Complete network --------------------------------------------------------
# generate Complete network as KNN network with the maximum number of neighbours
complete_net <- knearneigh(x = coord_urbanisation, 
                           k = 25, 
                           longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames())

# create igraph object 
complete_net_igraph <- neighborsDataFrame(complete_net) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
complete_net_gnar <- complete_net_igraph %>% 
  igraphtoGNAR()

# create order county index data frame 
county_index_complete <- data.frame("CountyName" = complete_net_igraph %>%
                                      V() %>% 
                                      names(), 
                                    "index" = seq(1, 26))

# compute network characteristics 
graph_complete <- network_characteristics(complete_net_igraph, 
                                          "complete")

# compute reproduction number 
R_complete_lv <- reproduction_rate_lv(complete_net_igraph)

# analyse sensitivity of the reproduction number against random attacks 
# absolute deletion of edges 
sense_complete <- reproduction_sensitivity(complete_net_igraph)
complete_overview <- sense_complete %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_complete_lv, 
            sd = sd(reproduction_rate))

# proportional deletion of edges 
sense_complete_percent <- reproduction_sensitivity_percentage(complete_net_igraph)
complete_percent_overview <- sense_complete_percent %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_complete_lv, 
            sd = sd(reproduction_rate))

# fit GNAR models for INV-D weighting 
results_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                             county_index = county_index_complete, 
                                             inverse_distance = TRUE)

return_best_model(results_complete)
# GNAR-4-1100-TRUE 191.9565
# we do not require the further distant neighbors 

# for latex 
strCaption <- "\\code{GNAR} models with INV-D weighting and without vertex
classification for the \\textbf{Complete} network and the corresponding BIC 
for model selection; the best performing model is highlighted in bold red"
print(xtable(results_complete[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_complete", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_complete[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# fit GNAR models for PB weighting
results_pop_weighted_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                          county_index = county_index_complete, 
                                                          inverse_distance = FALSE)

return_best_model(results_pop_weighted_complete)
# GNAR-5-11110-TRUE 193.7542

# fit GNAR models for INV-D weighting with classification 
results_class_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                   county_index = county_index_complete, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_class_complete)
# GNAR-5-10000-TRUE 194.119

# fit GNAR models for PB weighting with classification  
results_pop_weighted_class_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                                county_index = county_index_complete, 
                                                                inverse_distance = FALSE, 
                                                                weight_factor = urbanisation_factor, 
                                                                globalalpha = TRUE)
return_best_model(results_pop_weighted_class_complete)
# GNAR-5-10000-TRUE 194.5039

# fit GNAR models for SPL weighting 
results_old_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                 county_index = county_index_complete, 
                                                 old = TRUE)
return_best_model(results_old_complete)
# GNAR-5-11110-TRUE 192.7288

# fit GNAR models for SPL weighting with classification 
results_old_class_complete <- fit_and_predict_for_many(net = complete_net_gnar, 
                                                 county_index = county_index_complete, 
                                                 old = TRUE, 
                                                 weight_factor = urbanisation_factor, 
                                                 globalalpha = TRUE)
return_best_model(results_old_class_complete)
# GNAR-4-1000-TRUE 193.6324


# visualise change in BIC fit 
plot_BIC_for_GNAR(results_complete, 
                  results_class_complete, 
                  results_pop_weighted_complete, 
                  results_pop_weighted_class_complete, 
                  results_old_complete, 
                  results_old_class_complete)
ggsave("plots/modelfit/BIC/complete_bic.pdf", 
       width = 27, height = 14, unit = "cm")


# Network characteristics -------------------------------------------------
# Queen
graph_queen <- network_characteristics(covid_net_queen_igraph, 
                                            "queen")

# for latex 
strCaption <- "Summary for the \\textbf{Queen's contiguity} network, 
av. short for average, s.d. short for standard deviation"
print(xtable(graph_queen,
             digits=2,
             caption=strCaption,
             label="tab:summary_queen", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_queen)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Eco hub 
graph_eco_hub <- network_characteristics(covid_net_eco_hubs_igraph, 
                                         "eco_hub")

# for latex 
strCaption <- "Summary for the \\textbf{economic-hub} network, 
av. short for average, s.d. short for standard deviation"
print(xtable(graph_eco_hub,
             digits=2,
             caption=strCaption,
             label="tab:summary_eco_hub", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_eco_hub)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# Railway-based 
graph_train <- network_characteristics(covid_net_train_igraph, 
                                       "train")

# for latex 
strCaption <- "Summary for the \\textbf{rail-based} network, 
av. short for average, s.d. short for standard deviation"
print(xtable(graph_train,
             digits=2,
             caption=strCaption,
             label="tab:summary_train", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_train)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Delaunay triangulation
graph_delaunay <- network_characteristics(covid_net_delaunay_igraph, 
                                          "delaunay")

# for latex 
strCaption <- "Summary for the \\textbf{Delaunay triangulation} network, 
av. short for average, s.d. short for standard deviation"
print(xtable(graph_delaunay,
             digits=2,
             caption=strCaption,
             label="tab:summary_delaunay", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_delaunay)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Gabriel 
graph_gabriel <- network_characteristics(covid_net_gabriel_igraph, 
                                         "gabriel")

# for latex 
strCaption <- "Summary for the \\textbf{Gabriel} network, 
av. short for average, s.d. short for standard deviation"
print(xtable(graph_gabriel,
             digits=2,
             caption=strCaption,
             label="tab:summary_gabriel", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_gabriel)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Relative 
graph_relative <- network_characteristics(covid_net_relative_igraph, 
                                          "relative")

# for latex 
strCaption <- "Summary for the \\textbf{Relative neighbourhood} network, 
av. short for average, s.d. short for standard deviation"
print(xtable(graph_relative,
             digits=2,
             caption=strCaption,
             label="tab:summary_relative", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_relative)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Sphere of influence 
graph_soi <- network_characteristics(covid_net_soi_igraph, 
                                     "soi")

# for latex 
strCaption <- "Summary for the \\textbf{SOI} network, 
av. short for average, s.d. short for standard deviation"
print(xtable(graph_soi,
             digits=2,
             caption=strCaption,
             label="tab:summary_soi", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_soi)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# compare all networks (but KNN, DNN and Complete)
graph_overview <- cbind("metric" = graph_queen$metric,
                        "delauany" = graph_delaunay$delaunay %>% round(2), 
                        "gabriel" = graph_gabriel$gabriel %>% round(2), 
                        "soi" = graph_soi$soi %>% round(2),                         
                        "relative" = graph_relative$relative %>% round(2), 
                        "queen" = graph_queen$queen %>% round(2), 
                        "eco_hub"= graph_eco_hub$eco_hub %>% round(2), 
                        "train" = graph_train$train %>% round(2)
                        )

# for latex 
strCaption <- "Overview of network characteristics for \\textbf{Delaunay triangulation}, 
\\textbf{Gabriel}, \\textbf{SOI}, \\textbf{Relative neighbourhood (Rel. neigh.)}, 
\\textbf{Queen's contiguity}, \\textbf{Economic (Eco.) hub}, \\textbf{Railway-based} network,  
including average (av.) degree, density, average (av.) shortest path length (SPL), 
global and average (av.) local clustering (clust.) as well as average (av.) 
betweenness (betw.) and its standard deviation (s.d.)"
print(xtable(graph_overview,
             digits=2,
             caption=strCaption,
             label="tab:network_char", 
             align = c("", "l", "|", "r", "r", "r", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_soi)),
                        command = c(paste("\\toprule \n",
                                          " Metric & Delaunay & Gabriel & 
                                           SOI & Rel. neigh. & Queen & 
                                           Eco. hub & Railway \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# Spatial autocorrelation  ------------------------------------------------
# spatial autocorrelation measured as Moran's I for each network

# Queen 
moran_queen <- moran_I(data = COVID_weekly_data, 
                       nb_list = nb_list_queen)

ggplot(moran_queen, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"), 
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_queen.pdf", 
       width = 27, height = 14, unit = "cm")


# Economic hubs 
moran_eco_hubs <- moran_I(data = COVID_weekly_data, 
                       nb_list = nb_eco_hubs)

ggplot(moran_eco_hubs, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"),
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_eco_hubs.pdf", 
       width = 27, height = 14, unit = "cm")

# Railway-based
moran_train <- moran_I(data = COVID_weekly_data, 
                       nb_list = nb_list_train)

ggplot(moran_train, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"),
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_train.pdf", 
       width = 27, height = 14, unit = "cm")

# Delaunay 
moran_delaunay <- moran_I(data = COVID_weekly_data, 
                          nb_list = nb_list_delaunay)
ggplot(moran_delaunay, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I")  +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"),
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_delaunay.pdf", 
       width = 27, height = 14, unit = "cm")

# Gabriel 
moran_gabriel <- moran_I(data = COVID_weekly_data, 
                         nb_list = nb_list_gabriel)

ggplot(moran_gabriel, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"),
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_gabriel.pdf", 
       width = 27, height = 14, unit = "cm")

# Relative 
moran_relative <- moran_I(data = COVID_weekly_data, 
                          nb_list = nb_list_relative)

ggplot(moran_relative, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"),
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_relative.pdf", 
       width = 27, height = 14, unit = "cm")


# SOI
moran_soi <- moran_I(data = COVID_weekly_data, 
                     nb_list = nb_list_soi)
ggplot(moran_soi, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"),
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_soi.pdf", 
       width = 27, height = 14, unit = "cm")

# Scale-free --------------------------------------------------------------
# analyse log-log behaviour and regression R squared
queen_scale_free <- is_scale_free(covid_net_queen_igraph, 
                                  network_name = "queen")
queen_scale_free$graph
queen_scale_free$R_squared

eco_hub_scale_free <- is_scale_free(covid_net_eco_hubs_igraph, 
                                    network_name = "eco_hubs")
eco_hub_scale_free$graph
eco_hub_scale_free$R_squared

train_scale_free <- is_scale_free(covid_net_train_igraph, 
                                  network_name = "train")
train_scale_free$graph
train_scale_free$R_squared

knn_scale_free <- is_scale_free(opt_knn_net_igraph,
                                network_name = "knn")
knn_scale_free$graph
knn_scale_free$R_squared

dnn_scale_free <- is_scale_free(opt_dnn_net_igraph, 
                                network_name = "dnn")
dnn_scale_free$graph
dnn_scale_free$R_squared

delaunay_scale_free <- is_scale_free(covid_net_delaunay_igraph, 
                                     network_name = "delaunay")
delaunay_scale_free$graph
delaunay_scale_free$R_squared

gabriel_scale_free <- is_scale_free(covid_net_gabriel_igraph, 
                                    network_name = "gabriel")
gabriel_scale_free$graph
gabriel_scale_free$R_squared

relative_scale_free <- is_scale_free(covid_net_relative_igraph, 
                                     network_name = "relative")
relative_scale_free$graph
relative_scale_free$R_squared

soi_scale_free <- is_scale_free(igraph_net = covid_net_soi_igraph, 
                                network_name = "soi")
soi_scale_free$graph
soi_scale_free$R_squared

# for latex
scale_free_overview <- data.frame("name" = c("Delaunay triangulation", 
                                             "Gabriel", 
                                             "SOI",                                           
                                             "Relative neighbourhood", 
                                             "KNN", 
                                             "DNN",
                                             "Queen's", 
                                             "Economic hub", 
                                             "Railway-based"), 
                                  "R_squared" = c(delaunay_scale_free$R_squared, 
                                                  gabriel_scale_free$R_squared,
                                                  soi_scale_free$R_squared, 
                                                  relative_scale_free$R_squared,
                                                  knn_scale_free$R_squared,
                                                  dnn_scale_free$R_squared, 
                                                  queen_scale_free$R_squared, 
                                                  eco_hub_scale_free$R_squared, 
                                                  train_scale_free$R_squared))

strCaption <- "Test for scale-free property for each constructed network, 
$R^2$ for regression of log empirical cumulative distribution on log transformed 
degree"
print(xtable(scale_free_overview,
             digits=2,
             caption=strCaption,
             label="tab:scale_free", 
             align = c("", "l", "|", "c")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(scale_free_overview)),
                        command = c(paste("\\toprule \n",
                                          "Network & $R^2$ \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Reproduction rate -------------------------------------------------------

# compute reproduction number according to Pastor-Satorras and Vespignani 2002 
R_queen_psv <- reproduction_rate_psv(covid_net_queen_igraph)
R_economic_hub_psv <- reproduction_rate_psv(covid_net_eco_hubs_igraph)
R_train_psv <- reproduction_rate_psv(covid_net_train_igraph)
R_knn_psv <- reproduction_rate_psv(opt_knn_net_igraph)
R_dnn_psv <- reproduction_rate_psv(opt_dnn_net_igraph)
R_delaunay_psv <- reproduction_rate_psv(covid_net_delaunay_igraph)
R_gabriel_psv <- reproduction_rate_psv(covid_net_gabriel_igraph)
R_relative_psv <- reproduction_rate_psv(covid_net_relative_igraph)
R_soi_psv <- reproduction_rate_psv(covid_net_soi_igraph)

# create overview data frame 
R_comparison_psv <- data.frame("name" = c("Delaunay triangulation", 
                                          "Gabriel",
                                          "SOI", 
                                          "Relative neighbourhood", 
                                          "KNN", 
                                          "DNN",
                                          "Queen", 
                                          "Economic hub",
                                          "Railway-based"), 
                               "r_value" = c(R_delaunay_psv, 
                                             R_gabriel_psv, 
                                             R_soi_psv, 
                                             R_relative_psv,
                                             R_knn_psv, 
                                             R_dnn_psv,
                                             R_queen_psv,
                                             R_economic_hub_psv, 
                                             R_train_psv))

# compute reproduction number according to Lloyd and Valeika 2007 
R_queen_lv <- reproduction_rate_lv(covid_net_queen_igraph)
R_eco_hubs_lv <- reproduction_rate_lv(igraph_object = covid_net_eco_hubs_igraph, 
                                      numeric_vertices = TRUE, 
                                      county_index = county_index_eco_hubs)
R_train_lv <- reproduction_rate_lv(igraph_object = covid_net_train_igraph, 
                                   numeric_vertices = TRUE,
                                   county_index = county_index_train)
R_knn_lv <- reproduction_rate_lv(opt_knn_net_igraph)
R_dnn_lv <- reproduction_rate_lv(opt_dnn_net_igraph)
R_delaunay_lv <- reproduction_rate_lv(covid_net_delaunay_igraph)
R_gabriel_lv <- reproduction_rate_lv(covid_net_gabriel_igraph)
R_relative_lv <- reproduction_rate_lv(covid_net_relative_igraph)
R_soi_lv <- reproduction_rate_lv(covid_net_soi_igraph)

# create overview data frame 
R_comparison_lv <- data.frame("name" = c("Delaunay triangulation", 
                                         "Gabriel",
                                         "SOI", 
                                         "Relative neighbourhood", 
                                         "Complete",
                                         "KNN", 
                                         "DNN",
                                         "Queen", 
                                         "Economic hub", 
                                         "Railway-based"), 
                              "r_value" = c(R_delaunay_lv, 
                                            R_gabriel_lv, 
                                            R_soi_lv, 
                                            R_relative_lv,
                                            R_complete_lv, 
                                            R_knn_lv, 
                                            R_dnn_lv,
                                            R_queen_lv,
                                            R_eco_hubs_lv, 
                                            R_train_lv
                                            ))

# for latex
strCaption <- "Computation of basic reproduction number according to
\\cite{lloyd2007network} for all COVID-19 networks"
print(xtable(R_comparison_lv,
             digits=2,
             caption=strCaption,
             label="tab:rep_number_adapted", 
             align = c("", "l", "|", "c")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(R_comparison_lv)),
                        command = c(paste("\\toprule \n",
                                          "Network & $R_0$ \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# compute reproduction number according to Pastor-Satorras and Vespignani 2002 
R_queen_t <- reproduction_rate_t(covid_net_queen_igraph)
R_economic_hub_t <- reproduction_rate_t(igraph_object = covid_net_eco_hubs_igraph)
R_train_t <- reproduction_rate_t(covid_net_train_igraph)
R_knn_t <- reproduction_rate_t(opt_knn_net_igraph)
R_dnn_t <- reproduction_rate_t(opt_dnn_net_igraph)
R_delaunay_t <- reproduction_rate_t(covid_net_delaunay_igraph)
R_gabriel_t <- reproduction_rate_t(covid_net_gabriel_igraph)
R_relative_t <- reproduction_rate_t(covid_net_relative_igraph)
R_soi_t <- reproduction_rate_t(covid_net_soi_igraph)

# create overview data frame 
R_comparison_t <- data.frame("name" = c("Queen", 
                                        "Economic hub", 
                                        "Railway-based", 
                                        "KNN", 
                                        "DNN", 
                                        "Delaunay triangulation", 
                                        "Gabriel",
                                        "Relative neighbourhood", 
                                        "SOI"), 
                             "r_value" = c(R_queen_t,
                                           R_economic_hub_t, 
                                           R_train_t, 
                                           R_knn_t, 
                                           R_dnn_t, 
                                           R_delaunay_t, 
                                           R_gabriel_t, 
                                           R_relative_t,
                                           R_soi_t))




# Sensitivity of reproduction number --------------------------------------
# randomly delete 1-5 edges and observe how the reproduction number changes 
# repetition for 100 times 

# Queen 
sens_queen <- reproduction_sensitivity(covid_net_queen_igraph)
# compute average difference and standard deviation
queen_overview <- sens_queen %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_queen_lv, 
            sd = sd(reproduction_rate))

# Eco hubs 
sens_eco_hubs <- reproduction_sensitivity(covid_net_eco_hubs_igraph, 
                                          numeric_vertices = TRUE, 
                                          county_index = county_index_eco_hubs)
# compute average difference and standard deviation
eco_hubs_overview <- sens_eco_hubs %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_eco_hubs_lv, 
            sd = sd(reproduction_rate))


# Train
sens_train <- reproduction_sensitivity(covid_net_train_igraph, 
                                       numeric_vertices = TRUE, 
                                       county_index = county_index_train)
# compute average difference and standard deviation
train_overview <- sens_train %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_train_lv, 
            sd = sd(reproduction_rate))

# Delaunay
sens_delaunay <- reproduction_sensitivity(covid_net_delaunay_igraph)
# compute average difference and standard deviation
delaunay_overview <- sens_delaunay %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_delaunay_lv, 
            sd = sd(reproduction_rate))

# Gabriel 
sens_gabriel <- reproduction_sensitivity(covid_net_gabriel_igraph)
# compute average difference and standard deviation
gabriel_overview <- sens_gabriel %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_gabriel_lv, 
            sd = sd(reproduction_rate))

# Relative
sens_relative <- reproduction_sensitivity(covid_net_relative_igraph)
# compute average difference and standard deviation
relative_overview <- sens_relative %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_relative_lv, 
            sd = sd(reproduction_rate))

# SOI
sens_soi <- reproduction_sensitivity(covid_net_soi_igraph)
# compute average difference and standard deviation
soi_overview <- sens_soi %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_soi_lv, 
            sd = sd(reproduction_rate))


# KNN
sens_knn <- reproduction_sensitivity(opt_knn_net_igraph)
# compute average difference and standard deviation
knn_overview <- sens_knn %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_knn_lv, 
            sd = sd(reproduction_rate))

# DNN
sens_dnn <- reproduction_sensitivity(opt_dnn_net_igraph)
# compute average difference and standard deviation
dnn_overview <- sens_dnn %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_dnn_lv, 
            sd = sd(reproduction_rate))

# order according to density of network 
deletion_overview <- rbind(knn_overview, 
                           dnn_overview,
                           queen_overview, 
                           eco_hubs_overview, 
                           train_overview,
                           delaunay_overview, 
                           gabriel_overview, 
                           relative_overview, 
                           soi_overview, 
                           complete_overview)

deletion_overview$name <- c(rep("KNN", 5),
                            rep("DNN", 5), 
                            rep("Queen", 5), 
                            rep("Eco. hub", 5), 
                            rep("Railway-based", 5),  
                            rep("Delaunay", 5), 
                            rep("Gabriel", 5),
                            rep("Relative", 5),
                            rep("SOI", 5),
                            rep("Complete", 5)
)

# for latex 
strCaption <- "Overview over change in reproduction number $R_0$ for random edge
deletions, difference between mean $R_0$ and $R_0$ for initial network, 
standard deviation (sd) of new $R_0$,
random deletion of 1-5 edges repeated 100 times"
print(xtable(deletion_overview[, c(4, 1, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:rep_number_change_overview", 
             align = c("", "l", "|", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(deletion_overview)),
                        command = c(paste("\\toprule \n",
                                          " Network & \\# of deleted edges & 
                                          mean difference & sd \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# visualise sensitivity in a plot
ggplot(deletion_overview, aes(x = deletions, 
                              y = mean, 
                              color = name, 
                              group = name)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Deletion") +
  ylab(TeX("$\\Delta R_0$")) +
  guides(color = guide_legend(title = "Network"),
         shape = guide_legend(title = "Network"))  +
  theme(legend.position = "bottom")
ggsave("plots/reproductionNumber/sensitivity_absolute.pdf", 
       width = 23, height = 14, unit = "cm")


# Sensitivity of reproduction number - percentage -------------------------
# based on percentages, delete proportional to all edges
# repetition for 100 times 

# Queen 
sens_percent_queen <- reproduction_sensitivity_percentage(covid_net_queen_igraph)
# compute average difference and standard deviation
queen_percent_overview <- sens_percent_queen %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_queen_lv, 
            sd = sd(reproduction_rate))

# Eco hub
sens_percent_eco_hubs <- reproduction_sensitivity_percentage(covid_net_eco_hubs_igraph, 
                                                             numeric_vertices = TRUE, 
                                                             county_index = county_index_eco_hubs)
# compute average difference and standard deviation
eco_hubs_percent_overview <- sens_percent_eco_hubs %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_eco_hubs_lv, 
            sd = sd(reproduction_rate))

# Train 
sens_percent_train <- reproduction_sensitivity_percentage(covid_net_train_igraph, 
                                                          numeric_vertices = TRUE, 
                                                          county_index = county_index_train)
# compute average difference and standard deviation
train_percent_overview <- sens_percent_train %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_train_lv, 
            sd = sd(reproduction_rate))

# KNN
sens_percent_knn <- reproduction_sensitivity_percentage(opt_knn_net_igraph)
# compute average difference and standard deviation
knn_percent_overview <- sens_percent_knn %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_knn_lv, 
            sd = sd(reproduction_rate))

# DNN
sens_percent_dnn <- reproduction_sensitivity_percentage(opt_dnn_net_igraph)
# compute average difference and standard deviation
dnn_percent_overview <- sens_percent_dnn %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_dnn_lv, 
            sd = sd(reproduction_rate))


# Delaunay
sens_percent_delaunay <- reproduction_sensitivity_percentage(covid_net_delaunay_igraph)
# compute average difference and standard deviation
delaunay_percent_overview <- sens_percent_delaunay %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_delaunay_lv, 
            sd = sd(reproduction_rate))

# Gabriel 
sens_percent_gabriel <- reproduction_sensitivity_percentage(covid_net_gabriel_igraph)
# compute average difference and standard deviation
gabriel_percent_overview <- sens_percent_gabriel %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_gabriel_lv, 
            sd = sd(reproduction_rate))


# Relative 
sens_percent_relative <- reproduction_sensitivity_percentage(covid_net_relative_igraph)
# compute average difference and standard deviation
relative_percent_overview <- sens_percent_relative %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_relative_lv, 
            sd = sd(reproduction_rate))

# SOI
sens_percent_soi <- reproduction_sensitivity_percentage(covid_net_soi_igraph)
# compute average difference and standard deviation
soi_percent_overview <- sens_percent_soi %>% 
  group_by(deletions) %>% 
  summarise(mean = mean(reproduction_rate) - R_soi_lv, 
            sd = sd(reproduction_rate))



# order according to density of network 
deletion_percent_overview <- rbind(knn_percent_overview, 
                                   dnn_percent_overview,
                                   queen_percent_overview, 
                                   eco_hubs_percent_overview, 
                                   train_percent_overview,
                                   delaunay_percent_overview, 
                                   gabriel_percent_overview, 
                                   relative_percent_overview, 
                                   soi_percent_overview, 
                                   complete_percent_overview)

deletion_percent_overview$name <- c(rep("KNN", 4),
                                    rep("DNN", 4), 
                                    rep("Queen", 4), 
                                    rep("Eco. hub", 4), 
                                    rep("Rail-based", 4),  
                                    rep("Delaunay", 4), 
                                    rep("Gabriel", 4),
                                    rep("Relative", 4),
                                    rep("SOI", 4),
                                    rep("Complete", 4)
)

strCaption <- "Overview over change in reproduction number $R_0$ for random 
edge deletions, difference between mean $R_0$ and $R_0$ for initial network, 
standard deviation (sd) of new $R_0$, random deletion of 1\\%, 10\\%, 20\\% 
and 25\\% of edges in the network repeated 100 times"
print(xtable(deletion_percent_overview[, c(4, 1, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:rep_number_change_percent_overview", 
             align = c("", "l", "|", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(deletion_percent_overview)),
                        command = c(paste("\\toprule \n",
                                          " Network & \\# of deleted edges & mean difference & sd \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# visualise sensitivity in a plot
ggplot(deletion_percent_overview, 
       aes(x = deletions, 
           y = mean, 
           color = name, 
           group = name)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("% deletions") +
  ylab(TeX("$\\Delta R_0$")) +
  guides(color = guide_legend(title = "Network"),
         shape = guide_legend(title = "Network")) +
  theme(legend.position = "bottom")
ggsave("plots/reproductionNumber/sensitivity_percent.pdf", 
       width = 23, height = 14, unit = "cm")


# Weighting ---------------------------------------------------------------
# initialise data frame to store constructed weights and differences in 
# 1-lag COVID-19 ID
diff_all_counties <- data.frame("time" = NA, 
                                "CountyName" = NA, 
                                "cases" = NA, 
                                "dist" = NA,  
                                "pop_product" = NA, 
                                "pop_diff" = NA)

for (county in counties) {
  # construct data frame for  current county 
  cases_county <- covid_cases[, county]
  
  # compute difference in 1-lag COVID-19 ID to all the other counties 
  diff_cases <- covid_cases[, counties != county] - cases_county
  
  # format as a data frame with time column 
  diff_cases <- diff_cases %>% as.data.frame()
  diff_cases$time <- diff_cases %>% rownames()
  
  diff_cases_df <- diff_cases %>% gather("CountyName", "cases", -time)
  
  diff_county_names <- diff_cases %>% colnames()
  
  # filter row corresponding to the county in question from the
  # Great Circle distance data frame 
  dist_counties <- dist_urbanisation[county, 
                                     diff_county_names[diff_county_names != "time"]]
  
  # extract the population size for the county in question 
  county_pop <- population_weight[population_weight$CountyName == county, ]$weight
  
  
  
  dist_df <- data.frame("dist" = as.numeric(dist_counties[1, ]) / 1000, 
                        "CountyName" = dist_counties %>% colnames(), 
                        "pop_product" = population_weight[match(dist_counties %>% colnames(),
                                                         population_weight$CountyName), ]$weight * county_pop, 
                        "pop_diff" = population_weight[match(dist_counties %>% colnames(),
                                                             population_weight$CountyName), ]$weight - county_pop
                        ) 
  
  # assign county names and filter out only positive values to avoid double counting
  diff_all_counties <- rbind(diff_all_counties, 
                             left_join(diff_cases_df, 
                                       dist_df, 
                                       by = "CountyName")
                             ) %>% filter(cases >= 0)
}

# compute mean difference in 1-lag COVID-19 ID across time for each weighting 
# scheme
mean_diff_all_counties <- diff_all_counties %>% 
  group_by(CountyName, dist, pop_product, pop_diff) %>% 
  summarise(cases = mean(cases)) %>% 
  ungroup()  

# visualise the log transformed distance between counties against the 
# difference in 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = dist %>% log(), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log distance between counties") +
  ylab("Average difference in COVID-19 ID")
ggsave("plots/spatialCor/plot_distance_cases.pdf", 
       width = 23, height = 14, unit = "cm")

# visualise the log transformed INV-D weights against the difference in 
# 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = log(1 / dist), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log inverse distance between counties") +
  ylab("Average difference in COVID-19 ID")
ggsave("plots/spatialCor/plot_inverse_distance_cases.pdf", 
       width = 23, height = 14, unit = "cm")

# visualise the log transformed PB weights against the difference in 
# 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = log(pop_product / dist), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log population-based weights") +
  ylab("Average difference in COVID-19 ID")
ggsave("plots/spatialCor/plot_pop_weighted_cases.pdf", 
       width = 23, height = 14, unit = "cm")

# visualise the log transformed PD weights against the difference in 
# 1-lag COVID-19 ID 
ggplot(mean_diff_all_counties, 
       aes(x = log(abs(pop_diff) / dist), 
           y = cases)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE, method = "lm") +
  xlab("Log abs. difference in population sizes times inverse distance") +
  ylab("Average difference in COVID-19 ID")
ggsave("plots/spatialCor/plot_pop_difference_id_weight_cases.pdf", 
       width = 23, height = 14, unit = "cm")


# compute the shortest path length between vertices and the average 
# difference in 1-lag COVID-19 ID over time for each county-county pair
spl_mean_queen <- spl_and_mean_difference(covid_net_queen_igraph, 
                                          county_index_queen, 
                                          network_name = "queen")
spl_mean_eco_hubs <- spl_and_mean_difference(covid_net_eco_hubs_igraph, 
                                             county_index_eco_hubs,
                                             numeric_vertices = TRUE, 
                                             network_name = "eco_hubs")
spl_mean_train <- spl_and_mean_difference(covid_net_train_igraph, 
                                          county_index_train, 
                                          numeric_vertices = TRUE, 
                                          network_name = "train")
spl_mean_delaunay <- spl_and_mean_difference(covid_net_delaunay_igraph, 
                                             county_index_delaunay, 
                                             network_name = "delaunay")
spl_mean_gabriel <- spl_and_mean_difference(covid_net_gabriel_igraph, 
                                            county_index_gabriel, 
                                            network_name = "gabriel")
spl_mean_relative <- spl_and_mean_difference(covid_net_relative_igraph, 
                                             county_index_relative, 
                                             network_name = "relative")
spl_mean_soi <- spl_and_mean_difference(covid_net_soi_igraph, 
                                        county_index_soi, 
                                        network_name  = "soi")

# concatenate to one data frame 
spl_mean_all <- rbind(spl_mean_queen, 
                      spl_mean_eco_hubs, 
                      spl_mean_train, 
                      spl_mean_delaunay, 
                      spl_mean_gabriel, 
                      spl_mean_relative, 
                      spl_mean_soi)

# plot shortest path length against average difference in 1-lag COVID-19 ID
# filter out negative difference to avoid overcounting 
ggplot(spl_mean_all %>% filter(mean_diff >= 0),
       aes(x = spl,
           y = mean_diff)) +
  geom_point() +
  xlab("Shortest path length") +
  ylab("Average difference in COVID-19 ID")
ggsave("plots/spatialCor/plot_shortest_path_length_mean_difference.pdf",
       width = 14, height = 14, unit = "cm")



# Benchmark ARIMA  --------------------------------------------------------
# construct ARIMA models with optimal p for each county in isolation 
results_arima <- list()
for (county in covid_cases %>% colnames()) {
  
  covid_cases_county <- covid_cases_df %>% 
    dplyr::select(county %>% all_of()) %>% 
    ts(frequency = 52, start = c(2020, 1))
  
  # lag 1 to 5 tested 
  results_arima[[county]] <- sapply(seq(1, 5), FUN = function(i) {
    arima_mod <- arima(covid_cases_county, order = c(i, 0, 0))
    return(c(arima_mod %>% BIC(), 
             arima_mod %>% AIC()))
  })
}

mean_results_arima <- do.call(cbind.data.frame, results_arima)  %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate("p" = rep(seq(1, 5), 26)) 

colnames(mean_results_arima) <- c("BIC", "AIC", "p")

# compute mean BIC values across networks 
mean_results_arima_df <- mean_results_arima %>% 
  group_by(p) %>% 
  summarise(mean_BIC = mean(BIC), 
            mean_AIC = mean(AIC))

mean_results_arima_df$mean_BIC %>% which.min()
mean_results_arima_df$mean_BIC %>% min()
# lag 2 produces minimal BIC, 1482.536
mean_results_arima_df$mean_AIC %>% which.min()
mean_results_arima_df$mean_AIC %>% min()
# lag 4 produces minimal AIC, 1469.017


# for latex
strCaption <- "ARIMA models of order $p \\in {1, ..., 5}$ fitted for each county, 
mean BIC and AIC reported"
print(xtable(mean_results_arima_df,
             digits=2,
             caption=strCaption,
             label="tab:arima", 
             align = c("", "l", "|", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(mean_results_arima_df)),
                        command = c(paste("\\toprule \n",
                                          " order p & mean BIC & mean AIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)



# predict and compute MASE for ARIMA models 
mase_arima <- fit_and_predict_arima(lag = 2, 
                                    forecast_window = 10)


# Queen's contiguity ------------------------------------------------------

# fit GNAR models for INV-D weighting 
results_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                          county_index = county_index_queen, 
                                          inverse_distance = TRUE)
return_best_model(results_queen)
# GNAR-4-2211-TRUE 195.4917

# for latex 
strCaption <- "\\code{GNAR} models for Queen's case contiguity neighbourhood 
network, INV-D weighting without classification, BIC for model selection; 
the best performing model is highlighted in bold red"
print(xtable(results_queen[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_queens", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_queen[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# fit GNAR models for PB weighting 
results_pop_weighted_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                                       county_index = county_index_queen,
                                                       inverse_distance = FALSE, 
                                                       weight_index = population_weight)
return_best_model(results_pop_weighted_queen)
# GNAR-5-11100-TRUE 196.3142


# fit GNAR model for INV-D weighting with classification 
results_class_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                                             county_index = county_index_queen,
                                                             inverse_distance = TRUE, 
                                                             weight_index = population_weight, 
                                                             weight_factor = urbanisation_factor,
                                                             globalalpha = TRUE)
return_best_model(results_class_queen)
# GNAR-1-1-TRUE 196.315

# fit GNAR model for PB weighting with classification 
results_pop_weighted_class_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                                             county_index = county_index_queen,
                                                             inverse_distance = FALSE, 
                                                             weight_index = population_weight, 
                                                             weight_factor = urbanisation_factor,
                                                             globalalpha = TRUE)
return_best_model(results_pop_weighted_class_queen)
# GNAR-1-1-TRUE 196.6433


# fit GNAR models for SPL weighting 
results_old_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                              old = TRUE)
return_best_model(results_old_queen)
# GNAR-4-2211-TRUE 195.1655

# fit GNAR models for SPL weighting with classification 
results_old_class_queen <- fit_and_predict_for_many(net = covid_net_queen_gnar, 
                                              old = TRUE, 
                                              weight_factor = urbanisation_factor,
                                              globalalpha = TRUE)
return_best_model(results_old_class_queen)
# GNAR-1-1-TRUE 195.7187


# visualise BIC
plot_BIC_for_GNAR(results_queen,
                  results_class_queen, 
                  results_pop_weighted_queen, 
                  results_pop_weighted_class_queen, 
                  results_old_queen, 
                  results_old_class_queen)
ggsave("plots/modelfit/BIC/queen_bic.pdf", 
       width = 27, height = 14, unit = "cm")


# Economic hub ------------------------------------------------------------
# fit GNAR models for INV-D weighting 
results_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                             county_index = county_index_eco_hubs, 
                                             inverse_distance = TRUE, 
                                             numeric_vertices = TRUE)
return_best_model(results_eco_hubs)
# GNAR-4-2211-TRUE 195.4185

# for latex 
strCaption <- "\\code{GNAR} models for Economic hub network, 
INV-D weighting without classification, BIC for model selection, 
the best performing model is highlighted in bold red"
print(xtable(results_eco_hubs[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_economic_hubs", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_eco_hubs[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# fit GNAR models for PB weighting 
results_pop_weighted_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                          county_index = county_index_eco_hubs, 
                                                          inverse_distance = FALSE, 
                                                          numeric_vertices = TRUE)
return_best_model(results_pop_weighted_eco_hubs)
# GNAR-4-2211-TRUE 196.1491

# fit GNAR models for INV-D weighting with classification 
results_class_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                   county_index = county_index_eco_hubs, 
                                                   inverse_distance = TRUE, 
                                                   numeric_vertices = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_class_eco_hubs)
# GNAR-1-1-TRUE 196.1043


# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                   county_index = county_index_eco_hubs, 
                                                   inverse_distance = FALSE, 
                                                   numeric_vertices = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_pop_weighted_class_eco_hubs)
# GNAR-1-1-TRUE 196.5155

# fit GNAR models for SPL weighting  
results_old_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                 old = TRUE)
return_best_model(results_old_eco_hubs)
# GNAR-5-11100-TRUE 194.7234

# fit GNAR models for SPL weighting with classification 
results_old_class_eco_hubs <- fit_and_predict_for_many(net = covid_net_eco_hubs_gnar, 
                                                 old = TRUE, 
                                                 weight_factor = urbanisation_factor, 
                                                 globalalpha = TRUE)
return_best_model(results_old_class_eco_hubs)
# GNAR-1-1-TRUE 195.2224

# visualise BIC
plot_BIC_for_GNAR(results_eco_hubs, 
                  results_class_eco_hubs, 
                  results_pop_weighted_eco_hubs, 
                  results_pop_weighted_class_eco_hubs, 
                  results_old_eco_hubs, 
                  results_old_class_eco_hubs)
ggsave("plots/modelfit/BIC/eco_hubs_bic.pdf", 
       width = 27, height = 14, unit = "cm")

# Rail-based --------------------------------------------------------------
# fit GNAR models for INV-D weighting 
results_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                          county_index = county_index_train, 
                                          inverse_distance = TRUE, 
                                          numeric_vertices = TRUE)
return_best_model(results_train)
# GNAR-5-22111-TRUE 195.8426

# for latex
strCaption <- "\\code{GNAR} models for Railway-based network, 
counties represented by their county town, 
INV-D weighting without classification, BIC for model selection; 
the best performing model is highlighted in bold red"
print(xtable(results_train[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_rail_ct", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_train[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# fit GNAR models for PB weighting 
results_pop_weighted_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                          county_index = county_index_train, 
                                          inverse_distance = FALSE, 
                                          numeric_vertices = TRUE)
return_best_model(results_pop_weighted_train)
# GNAR-4-2211-TRUE 196.1854 --> many with only degree two

# fit GNAR models for INV-D weighting with classification  
results_class_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                                county_index = county_index_train, 
                                                inverse_distance = TRUE, 
                                                numeric_vertices = TRUE, 
                                                weight_factor = urbanisation_factor, 
                                                globalalpha = TRUE)
return_best_model(results_class_train)
# GNAR-1-0-TRUE 196.6496

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                                             county_index = county_index_train, 
                                                             inverse_distance = FALSE, 
                                                             numeric_vertices = TRUE, 
                                                             weight_factor = urbanisation_factor, 
                                                             globalalpha = TRUE)
return_best_model(results_pop_weighted_class_train)
# GNAR-1-0-TRUE 196.6496

# fit GNAR models for SPL weighting 
results_old_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                              old = TRUE)
return_best_model(results_old_train)
# GNAR-4-2111-TRUE 194.9313

# fit GNAR models for SPL weighting with classification 
results_old_class_train <- fit_and_predict_for_many(net = covid_net_train_gnar, 
                                                    old = TRUE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
return_best_model(results_old_class_train)
# GNAR-1-1-TRUE 196.1892


# visualise BIC 
plot_BIC_for_GNAR(results_train, 
                  results_class_train, 
                  results_pop_weighted_train, 
                  results_pop_weighted_class_train, 
                  results_old_train, 
                  results_old_class_train)
ggsave("plots/modelfit/BIC/train_bic.pdf", 
       width = 27, height = 14, unit = "cm")

# Delaunay triangulation --------------------------------------------------
# fit GNAR models for INV-D weighting 
results_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                             county_index = county_index_delaunay, 
                                             inverse_distance = TRUE)
return_best_model(results_delaunay)
# GNAR-4-2211-TRUE 195.1285

# for latex 
strCaption <- "\\code{GNAR} models for Delaunay triangulation network, 
INV-D weighting without classification, BIC for model selection; 
the best performing model is highlighted in bold red"
print(xtable(results_delaunay[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_delaunay", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_delaunay[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# fit GNAR models for PB weighting 
results_pop_weighted_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                          county_index = county_index_delaunay, 
                                                          inverse_distance = FALSE)
return_best_model(results_pop_weighted_delaunay)
# GNAR-4-2211-TRUE 195.5891

# fit GNAR models for INV-D weighting with classification 
results_class_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                   county_index = county_index_delaunay, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
return_best_model(results_class_delaunay)
# GNAR-5-11000-TRUE 196.1148

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                                county_index = county_index_delaunay, 
                                                                inverse_distance = FALSE, 
                                                                weight_factor = urbanisation_factor,
                                                                globalalpha = TRUE)
return_best_model(results_pop_weighted_class_delaunay)
# GNAR-2-22-TRUE 196.4688


# fit GNAR models for SPL weighting 
results_old_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                              old = TRUE)
return_best_model(results_old_delaunay)
# GNAR-4-2111-TRUE 192.6884

# fit GNAR models for SPL weighting with classification 
results_old_class_delaunay <- fit_and_predict_for_many(net = covid_net_delaunay_gnar, 
                                                       old = TRUE, 
                                                       weight_factor = urbanisation_factor, 
                                                       globalalpha = TRUE)
return_best_model(results_old_class_delaunay)
# GNAR-1-1-TRUE 195.1826

# visualise BIC 
plot_BIC_for_GNAR(results_delaunay, 
                  results_class_delaunay, 
                  results_pop_weighted_delaunay, 
                  results_pop_weighted_class_delaunay, 
                  results_old_delaunay, 
                  results_old_class_delaunay)
ggsave("plots/modelfit/BIC/delaunay_bic.pdf", 
       width = 27, height = 14, unit = "cm")

# Gabriel neighbourhood  --------------------------------------------------
# fit GNAR models for INV-D weighting 
results_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                            county_index = county_index_gabriel, 
                                            inverse_distance = TRUE)
return_best_model(results_gabriel)
# GNAR-5-11000-TRUE 195.2479

# for latex
strCaption <- "\\code{GNAR} models for Gabriel network, 
INV-D weighting without classification, BIC for model selection; 
the best performing model is highlighted in bold red"
print(xtable(results_gabriel[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_gabriel", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_gabriel[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# fit GNAR models for PB weighting 
results_pop_weighted_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                            county_index = county_index_gabriel, 
                                            inverse_distance = FALSE)
return_best_model(results_pop_weighted_gabriel)
# GNAR-4-2221-TRUE 195.7398


# fit GNAR models for INV-D weighting with classification 
results_class_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                            county_index = county_index_gabriel, 
                                            inverse_distance = TRUE, 
                                            weight_factor = urbanisation_factor, 
                                            globalalpha = TRUE)
return_best_model(results_class_gabriel)
# GNAR-1-1-TRUE 195.5572

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                                  county_index = county_index_gabriel, 
                                                  inverse_distance = FALSE, 
                                                  weight_factor = urbanisation_factor, 
                                                  globalalpha = TRUE)
return_best_model(results_pop_weighted_class_gabriel)
# GNAR-1-1-TRUE 196.0854


# fit GNAR models for SPL weighting 
results_old_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                                 old = TRUE)
return_best_model(results_old_gabriel)
# GNAR-5-11100-TRUE 195.0627

# fit GNAR models for SPL weighting with classification 
results_old_class_gabriel <- fit_and_predict_for_many(net = covid_net_gabriel_gnar, 
                                                      old = TRUE, 
                                                      weight_factor = urbanisation_factor, 
                                                      globalalpha = TRUE)
return_best_model(results_old_class_gabriel)
# GNAR-1-1-TRUE 195.5408

# visualise BIC 
plot_BIC_for_GNAR(results_gabriel, 
                  results_class_gabriel, 
                  results_pop_weighted_gabriel, 
                  results_pop_weighted_class_gabriel, 
                  results_old_gabriel, 
                  results_old_class_gabriel)
ggsave("plots/modelfit/BIC/gabriel_bic.pdf", 
       width = 27, height = 14, unit = "cm")

# Relative neighbourhood --------------------------------------------------
# fit GNAR models for INV-D weighting 
results_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = TRUE)
return_best_model(results_relative)
# GNAR-4-2211-TRUE 195.9906

# for latex 
strCaption <- "\\code{GNAR} models for Relative neighbourhood network, 
INV-D weighting without classification, BIC for model selection; 
the best performing model is highlighted in bold red"
print(xtable(results_relative[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_relative", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_relative[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# fit GNAR models for PB weighting 
results_pop_weighted_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = FALSE)
return_best_model(results_pop_weighted_relative)
# GNAR-1-1-TRUE 196.4352

# fit GNAR models for INV-D weighting with classification  
results_class_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = TRUE, 
                                             weight_factor = urbanisation_factor, 
                                             globalalpha = TRUE)
return_best_model(results_class_relative)
# GNAR-1-1-TRUE 196.4462

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                             county_index = county_index_relative, 
                                             inverse_distance = FALSE, 
                                             weight_factor = urbanisation_factor, 
                                             globalalpha = TRUE)
return_best_model(results_pop_weighted_class_relative)
# GNAR-1-1-TRUE 196.498

# fit GNAR models for SPL weighting 
results_old_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                                 old = TRUE)
return_best_model(results_old_relative)
# GNAR-4-2211-TRUE 195.988

# fit GNAR models for SPL weighting with classification 
results_old_class_relative <- fit_and_predict_for_many(net = covid_net_relative_gnar, 
                                                 old = TRUE, 
                                                 weight_factor = urbanisation_factor, 
                                                 globalalpha = TRUE)
return_best_model(results_old_class_relative)
# GNAR-1-1-TRUE 196.1255

# visualise BIC
plot_BIC_for_GNAR(results_relative, 
                  results_class_relative, 
                  results_pop_weighted_relative, 
                  results_pop_weighted_class_relative, 
                  results_old_relative,
                  results_old_class_relative)
ggsave("plots/modelfit/BIC/relative_bic.pdf", 
       width = 27, height = 14, unit = "cm")

# Sphere of influence neighbourhood ---------------------------------------
# fit GNAR models for INV_D weighting 
results_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = TRUE)
return_best_model(results_soi)
# GNAR-4-2221-TRUE 194.9368

# for latex
strCaption <- "\\code{GNAR} models for SOI network, 
INV-D weighting without classification, 
BIC for model selection; the best performing model is highlighted in bold red"
print(xtable(results_soi[, c(1, 3)],
             digits=2,
             caption=strCaption,
             label="tab:gnar_soi", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(results_soi[, c(1, 3)])),
                        command = c(paste("\\toprule \n",
                                          " Model name & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# fit GNAR models for PB weighting 
results_pop_weighted_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = FALSE)
return_best_model(results_pop_weighted_soi)
# GNAR-5-11100-TRUE 196.0015


# fit GNAR models for INV-D weighting with classification  
results_class_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = TRUE, 
                                        weight_factor = urbanisation_factor, 
                                        globalalpha = TRUE)
return_best_model(results_class_soi)
# GNAR-1-1-TRUE 196.0503

# fit GNAR models for PB weighting with classification 
results_pop_weighted_class_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                        county_index = county_index_soi,
                                        inverse_distance = FALSE, 
                                        weight_factor = urbanisation_factor, 
                                        globalalpha = TRUE)
return_best_model(results_pop_weighted_class_soi)
# GNAR-1-1-TRUE 196.3516

# fit GNAR models for SPL weighting 
results_old_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                                 old = TRUE)
return_best_model(results_old_soi)
# GNAR-4-2111-TRUE 194.2081

# fit GNAR models for SPL weighting with classification 
results_old_class_soi <- fit_and_predict_for_many(net = covid_net_soi_gnar, 
                                            old = TRUE, 
                                            weight_factor = urbanisation_factor, 
                                            globalalpha = TRUE)
return_best_model(results_old_class_soi)
# GNAR-5-21111-TRUE 195.2722

# visualise BIC
plot_BIC_for_GNAR(results_soi, 
                  results_class_soi, 
                  results_pop_weighted_soi, 
                  results_pop_weighted_class_soi, 
                  results_old_soi, 
                  results_old_class_soi)
ggsave("plots/modelfit/BIC/soi_bic.pdf", 
       width = 27, height = 14, unit = "cm")


# KNN ---------------------------------------------------------------------
# up to fully connected 
max_k <- dim(covid_cases)[2] - 1

# create list to save best performing model for each neighbourhood size 
knn_best <- list()

for (k in seq(1, max_k, by = 2)) {
  # create nb list for neighbourhood size k
  nb_knn <- knearneigh(x = coord_urbanisation, 
                       k = k, 
                       longlat = TRUE) %>% 
    knn2nb(row.names = coord_urbanisation %>% row.names())
  
  # create igraph object 
  covid_net_knn_igraph <- neighborsDataFrame(nb = nb_knn) %>% 
    graph_from_data_frame(directed = FALSE) %>% 
    igraph::simplify() 
  
  # create GNAR object 
  covid_net_knn <- covid_net_knn_igraph %>% 
    igraphtoGNAR()
  
  # create ordered coutny index data frame 
  county_index_knn <- data.frame("CountyName" = covid_net_knn_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))
  # fit GNAR models for INV-D weighting 
  results_knn_id <- fit_and_predict_for_many(net = covid_net_knn, 
                                             county_index = county_index_knn, 
                                             inverse_distance = TRUE)
  # fit GNAR models for PB weighting 
  results_knn_pop <- fit_and_predict_for_many(net = covid_net_knn, 
                                              county_index = county_index_knn, 
                                              inverse_distance = FALSE)
  
  # fit GNAR models for INV-D weighting with classification
  results_knn_id_urban <- fit_and_predict_for_many(net = covid_net_knn, 
                                                   county_index = county_index_knn, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
  
  # fit GNAR models for PB weighting with classification 
  results_knn_pop_urban <- fit_and_predict_for_many(net = covid_net_knn, 
                                                    county_index = county_index_knn, 
                                                    inverse_distance = FALSE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  # # fit GNAR models for SPL weighting 
  results_knn_old <- fit_and_predict_for_many(net = covid_net_knn, 
                                              county_index = county_index_knn, 
                                              old = TRUE)
  
  # fit GNAR models for SPL weighting with classification 
  results_knn_old_class <- fit_and_predict_for_many(net = covid_net_knn, 
                                                    county_index = county_index_knn, 
                                                    old = TRUE, 
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  
  # save best performing model for every k
  knn_best[[length(knn_best) + 1]] <- data.frame("k" = k, 
                                                 "best model id" = return_best_model(results_knn_id), 
                                                 "best model pop" = return_best_model(results_knn_pop),
                                                 "best model id class" = return_best_model(results_knn_id_urban), 
                                                 "best model pop class" = return_best_model(results_knn_pop_urban), 
                                                 "best model old" = return_best_model(results_knn_old), 
                                                 "best model old class" = return_best_model(results_knn_old_class)
  )
  
}

knn_best_df <- do.call(rbind.data.frame, knn_best)

# id - id.class - pop - pop.class - old - old.class
knn_best_res <- return_best_knn_dnn(knn_best_df)
# 1 21  GNAR-5-11111-TRUE  193.224464618423
# 2 11   GNAR-4-2221-TRUE  193.332889189113
# 3 25  GNAR-5-10000-TRUE  194.118994421962
# 4 15  GNAR-5-10000-TRUE  194.377915086384
# 5 21   GNAR-4-1100-TRUE  192.533347515879
# 6 21   GNAR-4-1000-TRUE  193.415126643114

knn_best_res$best.model.id.BIC <- knn_best_res$best.model.id.BIC %>% as.numeric()
knn_best_res$weighting <- c("INV-D", 
                            "INV-D+class",
                            "PB", 
                            "PB+class", 
                            "SPL", 
                            "SPL+class")


# for latex
strCaption <- "Overview over best \\code{GNAR} model for different weighting 
schemes, k denotes size of neighbourhood, BIC for model selection; 
the best performing model is highlighted in bold red"
print(xtable(knn_best_res[, c(4, 1, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:performance_knn", 
             align = c("", "l", "|", "r", "c", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(knn_best_res[, c(4, 1, 2, 3)])),
                        command = c(paste("\\toprule \n",
                                          "Weighting & size k & \\code{GNAR} model & 
                                          BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# plot k against BIC 
colnames(knn_best_df) <- c("k", 
                           "INV-D_name", 
                           "INV-D", 
                           "PB_name",
                           "PB", 
                           "INV-D_class_name", 
                           "INV-D+class", 
                           "PB_class_name", 
                           "PB+class",
                           "SPL_name", 
                           "SPL", 
                           "SPL_class_name", 
                           "SPL+class")

knn_plot <- knn_best_df %>% 
  dplyr::select("k",
                "INV-D",
                "PB",
                "INV-D+class",
                "PB+class",
                "SPL",
                "SPL+class") %>%
  gather("weight", "BIC", -k)


ggplot(knn_plot, 
       aes(x = k, 
           y = BIC, 
           group = weight, 
           color = weight)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Neighbourhood size k") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") 
ggsave("plots/hyperparameter/knn/k_vs_bic.pdf", 
       width = 26, height = 14, unit = "cm")


# DNN ---------------------------------------------------------------------
# find thresholds for DNN network 
max_dist_threshold <- dist_urbanisation %>% colMax() %>% max()
dist_urbanisation %>% colMax() %>% which.max()

min_dist_threshold <- dist_urbanisation %>% colMin() %>% max()

dist_urbanisation %>% colMin() %>% mean() # 42km
dist_urbanisation %>% colMin() %>% sd() # 18km

# create list to save best performing model for each distance threshold d
dnn_best <- list()
for (d in seq(100, 
              350, 
              by = 25)) {
  # create nb list 
  nb_dnn <- dnearneigh(x = coord_urbanisation, 
                       d1 = 0, 
                       d2 = d,
                       row.names = coord_urbanisation %>% rownames(),
                       longlat = TRUE, 
                       use_s2 = TRUE)
  
  # create igraph object
  covid_net_dnn_igraph <- neighborsDataFrame(nb = nb_dnn) %>% 
    graph_from_data_frame(directed = FALSE) %>% 
    igraph::simplify() 
  
  # create GNAR object 
  covid_net_dnn <- covid_net_dnn_igraph %>% 
    igraphtoGNAR()
  
  # create ordered county index data frame 
  county_index_dnn <- data.frame("CountyName" = covid_net_dnn_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))
  
  # fit GNAR models for INV-D weighting 
  results_dnn_id <- fit_and_predict_for_many(net = covid_net_dnn, 
                                             county_index = county_index_dnn, 
                                             inverse_distance = TRUE)
  
  # fit GNAR models for PB weighting
  results_dnn_pop <- fit_and_predict_for_many(net = covid_net_dnn, 
                                              county_index = county_index_dnn, 
                                              inverse_distance = FALSE)
  
  # fit GNAR models for INV-D weighting with classification 
  results_dnn_id_urban <- fit_and_predict_for_many(net = covid_net_dnn, 
                                                   county_index = county_index_dnn, 
                                                   inverse_distance = TRUE, 
                                                   weight_factor = urbanisation_factor, 
                                                   globalalpha = TRUE)
  
  # fit GNAR models for PB weighting with classification 
  results_dnn_pop_urban <- fit_and_predict_for_many(net = covid_net_dnn, 
                                                    county_index = county_index_dnn, 
                                                    inverse_distance = FALSE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  # fit GNAR models for SPL weighting
  results_dnn_old <- fit_and_predict_for_many(net = covid_net_dnn, 
                                              county_index = county_index_dnn, 
                                              old = TRUE)
  
  # fit GNAR models for SPL weighting with classification 
  results_dnn_old_class <- fit_and_predict_for_many(net = covid_net_dnn, 
                                                    county_index = county_index_dnn, 
                                                    old = TRUE,
                                                    weight_factor = urbanisation_factor, 
                                                    globalalpha = TRUE)
  
  # save best performing model for every d
  dnn_best[[length(dnn_best) + 1]] <- data.frame("d in km" = d, 
                                                 "best model id" = return_best_model(results_dnn_id), 
                                                 "best model pop" = return_best_model(results_dnn_pop),
                                                 "best model id class" = return_best_model(results_dnn_id_urban), 
                                                 "best model pop class" = return_best_model(results_dnn_pop_urban),
                                                 "best model old" = return_best_model(results_dnn_old), 
                                                 "best model old class" = return_best_model(results_dnn_old_class)
  )
}

dnn_best_df <- do.call(rbind.data.frame, dnn_best)

# id - id.class - pop - pop.class - old - old.class 
dnn_best_res <- return_best_knn_dnn(dnn_best_df)
# 1     175   GNAR-4-2221-TRUE  192.658054255805
# 2     200  GNAR-5-11110-TRUE  193.145995455023
# 3     175     GNAR-2-22-TRUE  193.661264024826
# 4     200     GNAR-2-11-TRUE  193.889497530882
# 5     250  GNAR-5-11110-TRUE  192.692730321918
# 6     250  GNAR-5-10000-TRUE   193.61360573789

dnn_best_res$best.model.id.BIC <- dnn_best_res$best.model.id.BIC %>% as.numeric()

dnn_best_res$weighting <- c("INV-d", 
                            "INV-D+class",
                            "PB", 
                            "PB+class", 
                            "SPL", 
                            "SPL+class")


# for latex 
strCaption <- "Overview over best \\code{GNAR} model for different weighting schemes, 
BIC for model selection; the best performing model is highlighted in bold red"
print(xtable(dnn_best_res[, c(4, 1, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:performance_dnn", 
             align = c("", "l", "|", "r", "c", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(dnn_best_res[, c(4, 1, 2, 3)])),
                        command = c(paste("\\toprule \n",
                                          "Weighting scheme & distance d (in km) & 
                                          \\code{GNAR} model & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# plot d against BIC values 
colnames(dnn_best_df) <- c("d", 
                           "INV-D_name", 
                           "INV-D", 
                           "PB_name",
                           "PB", 
                           "INV-D_class_name", 
                           "INV-D+class", 
                           "PB_class_name", 
                           "PB+class",
                           "SPL_name", 
                           "SPL", 
                           "SPL_class_name", 
                           "SPL+class")

dnn_plot <- dnn_best_df %>% dplyr::select("d", 
                                          "INV-D", 
                                          "PB", 
                                          "INV-D+class", 
                                          "PB+class", 
                                          "SPL", 
                                          "SPL+class") %>% 
  gather("weight", "BIC", -d)


ggplot(dnn_plot, 
       aes(x = d, 
           y = BIC, 
           group = weight, 
           color = weight)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Distance threshold d") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") 
ggsave("plots/hyperparameter/dnn/d_vs_bic.pdf", 
       width = 26, height = 14, unit = "cm")


# Compare BIC for GNAR models ---------------------------------------------
# for latex 
bic_development <- cbind(results_delaunay[, c(1, 3)], 
                         results_gabriel[, 3], 
                         results_soi[, 3], 
                         results_relative[, 3], 
                         results_queen[, 3], 
                         results_eco_hubs[, 3], 
                         results_train[, 3])


strCaption <- "\\code{GNAR} models with INV-D weighting and without vertex classification
for \\textbf{Delaunay triangulation}, \\textbf{Gabriel}, \\textbf{SOI}, 
\\textbf{Relative neighbourhood}, \\textbf{Queen's contiguity}, 
\\textbf{Economic hub} and \\textbf{Railway-based} network and the corresponding
BIC for model selection; the best performing model for each network is 
highlighted in bold red"
print(xtable(bic_development,
             digits=2,
             caption=strCaption,
             label="tab:gnar_all", 
             align = c("", "l", "|", "r", "r", "r", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(bic_development)),
                        command = c(paste("\\toprule \n",
                                          " Model name & Delaunay & Gabriel & 
                                          SOI & Relative & Queen & Eco. hub & 
                                          Railway \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)



# Performance overview ----------------------------------------------------
compare_df <- data.frame("network" = c(rep("Delaunay", 6),
                                       rep("Gabriel", 6),
                                       rep("SOI", 6),
                                       rep("Relative", 6), 
                                       rep("Complete", 6), 
                                       rep("Queen", 6),
                                       rep("Eco. hub", 6), 
                                       rep("Rail", 6)
                                       ), 
                         "hyperparameters" = rep(c("INV-D", 
                                                   "INV-D+class", 
                                                   "PB", 
                                                   "PB+class", 
                                                   "SPL", 
                                                   "SPL+class"), 
                                                 8), 
                         "model" = c(
                                     return_best_model(results_delaunay)$name, 
                                     return_best_model(results_class_delaunay)$name, 
                                     return_best_model(results_pop_weighted_delaunay)$name,
                                     return_best_model(results_pop_weighted_class_delaunay)$name,
                                     return_best_model(results_old_delaunay)$name,
                                     return_best_model(results_old_class_delaunay)$name,
                                     
                                     return_best_model(results_gabriel)$name, 
                                     return_best_model(results_class_gabriel)$name, 
                                     return_best_model(results_pop_weighted_gabriel)$name,
                                     return_best_model(results_pop_weighted_class_gabriel)$name,
                                     return_best_model(results_old_gabriel)$name,
                                     return_best_model(results_old_class_gabriel)$name,
                                     
                                     return_best_model(results_soi)$name, 
                                     return_best_model(results_class_soi)$name, 
                                     return_best_model(results_pop_weighted_soi)$name,
                                     return_best_model(results_pop_weighted_class_soi)$name, 
                                     return_best_model(results_old_soi)$name,
                                     return_best_model(results_old_class_soi)$name,
                                     
                                     return_best_model(results_relative)$name, 
                                     return_best_model(results_class_relative)$name, 
                                     return_best_model(results_pop_weighted_relative)$name,
                                     return_best_model(results_pop_weighted_class_relative)$name,
                                     return_best_model(results_old_relative)$name,
                                     return_best_model(results_old_class_relative)$name,
                                     
                                     return_best_model(results_complete)$name, 
                                     return_best_model(results_class_complete)$name, 
                                     return_best_model(results_pop_weighted_complete)$name,
                                     return_best_model(results_pop_weighted_class_complete)$name,
                                     return_best_model(results_old_complete)$name, 
                                     return_best_model(results_old_class_complete)$name, 
                                     
                                     return_best_model(results_queen)$name, 
                                     return_best_model(results_class_queen)$name, 
                                     return_best_model(results_pop_weighted_queen)$name,
                                     return_best_model(results_pop_weighted_class_queen)$name,
                                     return_best_model(results_old_queen)$name,
                                     return_best_model(results_old_class_queen)$name, 
                                     
                                     return_best_model(results_eco_hubs)$name, 
                                     return_best_model(results_class_eco_hubs)$name, 
                                     return_best_model(results_pop_weighted_eco_hubs)$name,
                                     return_best_model(results_pop_weighted_class_eco_hubs)$name,
                                     return_best_model(results_old_eco_hubs)$name,
                                     return_best_model(results_old_class_eco_hubs)$name,
                                     
                                     return_best_model(results_train)$name, 
                                     return_best_model(results_class_train)$name, 
                                     return_best_model(results_pop_weighted_train)$name,
                                     return_best_model(results_pop_weighted_class_train)$name,
                                     return_best_model(results_old_train)$name,
                                     return_best_model(results_old_class_train)$name
                                     ), 
                         "BIC" = c(
                                   return_best_model(results_delaunay)$BIC, 
                                   return_best_model(results_class_delaunay)$BIC, 
                                   return_best_model(results_pop_weighted_delaunay)$BIC,
                                   return_best_model(results_pop_weighted_class_delaunay)$BIC,
                                   return_best_model(results_old_delaunay)$BIC,
                                   return_best_model(results_old_class_delaunay)$BIC,
                                   
                                   return_best_model(results_gabriel)$BIC, 
                                   return_best_model(results_class_gabriel)$BIC, 
                                   return_best_model(results_pop_weighted_gabriel)$BIC,
                                   return_best_model(results_pop_weighted_class_gabriel)$BIC,
                                   return_best_model(results_old_gabriel)$BIC,
                                   return_best_model(results_old_class_gabriel)$BIC,
                                   
                                   return_best_model(results_soi)$BIC, 
                                   return_best_model(results_class_soi)$BIC, 
                                   return_best_model(results_pop_weighted_soi)$BIC,
                                   return_best_model(results_pop_weighted_class_soi)$BIC, 
                                   return_best_model(results_old_soi)$BIC,
                                   return_best_model(results_old_class_soi)$BIC,
                                   
                                   return_best_model(results_relative)$BIC, 
                                   return_best_model(results_class_relative)$BIC, 
                                   return_best_model(results_pop_weighted_relative)$BIC,
                                   return_best_model(results_pop_weighted_class_relative)$BIC,
                                   return_best_model(results_old_relative)$BIC,
                                   return_best_model(results_old_class_relative)$BIC,
                                   
                                   return_best_model(results_complete)$BIC, 
                                   return_best_model(results_class_complete)$BIC, 
                                   return_best_model(results_pop_weighted_complete)$BIC,
                                   return_best_model(results_pop_weighted_class_complete)$BIC, 
                                   return_best_model(results_old_complete)$BIC, 
                                   return_best_model(results_old_class_complete)$BIC,
                                   
                                   return_best_model(results_queen)$BIC, 
                                   return_best_model(results_class_queen)$BIC, 
                                   return_best_model(results_pop_weighted_queen)$BIC,
                                   return_best_model(results_pop_weighted_class_queen)$BIC,
                                   return_best_model(results_old_queen)$BIC, 
                                   return_best_model(results_old_class_queen)$BIC,
                                   
                                   return_best_model(results_eco_hubs)$BIC, 
                                   return_best_model(results_class_eco_hubs)$BIC, 
                                   return_best_model(results_pop_weighted_eco_hubs)$BIC,
                                   return_best_model(results_pop_weighted_class_eco_hubs)$BIC,
                                   return_best_model(results_old_eco_hubs)$BIC, 
                                   return_best_model(results_old_class_eco_hubs)$BIC, 
                                   
                                   return_best_model(results_train)$BIC, 
                                   return_best_model(results_class_train)$BIC, 
                                   return_best_model(results_pop_weighted_train)$BIC,
                                   return_best_model(results_pop_weighted_class_train)$BIC,
                                   return_best_model(results_old_train)$BIC,
                                   return_best_model(results_old_class_train)$BIC
                         ))

# for latex 
strCaption <- "Comparison of best performing \\code{GNARfit()} models for 
COVID-19 networks, for all weighting schemes with and without classification 
(+class) according to urbanisation of county; 
the best performing model is highlighted in bold red"
print(xtable(compare_df,
             digits=2,
             caption=strCaption,
             label="tab:overview_weighting_best_model", 
             align = c("", "l", "|", "r", "c", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(compare_df)),
                        command = c(paste("\\toprule \n",
                                          " Network & hyperparameter & \\code{GNAR} model & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# for KNN and DNN
compare_knn_dnn <- data.frame("network" = c(rep("KNN", 6), 
                                            rep("DNN", 6)), 
                              "hyperparameters" = rep(c("INV-D", 
                                                        "INV-D+class", 
                                                        "PB", 
                                                        "PB+class", 
                                                        "SPL", 
                                                        "SPL+class"), 
                                                      2), 
                              "k_d" = c(knn_best_res$k, 
                                        dnn_best_res$d.in.km), 
                              "model" = c(knn_best_res$best.model.id.name, 
                                          dnn_best_res$best.model.id.name), 
                              "BIC" = c(knn_best_res$best.model.id.BIC, 
                                        dnn_best_res$best.model.id.BIC) %>% 
                                as.numeric() %>% 
                                round(2)
                              )


# for latex 
strCaption <- "Comparison of best performing \\code{GNARfit()} models for KNN 
and DNN networks, for all weighting schemes with and without classification 
(+class) according to urbanisation of county; 
the best performing model is highlighted in bold red"
print(xtable(compare_knn_dnn,
             digits=2,
             caption=strCaption,
             label="tab:overview_best_model_knn_dnn", 
             align = c("", "l", "|", "r", "r", "c", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(compare_knn_dnn)),
                        command = c(paste("\\toprule \n",
                                          " Network & hyperparameter & k / d [km] & 
                                          \\code{GNAR} model & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# best KNN / DNN 
compare_knn_dnn %>% 
  group_by(network) %>% 
  summarise(which.min(BIC))
# first model for both KNN / DNN is best performing 

# absolute best for KNN / DNN 
compare_knn_dnn[which.min(compare_knn_dnn$BIC), ]
# KNN             spl  21 GNAR-4-1100-TRUE 192.53

# absolute best 
compare_df[which.min(compare_df$BIC), ]
# Complete             spl GNAR-5-11110-TRUE 192.7288

# best sparse network 
compare_df[compare_df$network != "Complete", ][which.min(compare_df[
  compare_df$network != "Complete", ]$BIC), 
  ]
# Delaunay             spl GNAR-4-2211-TRUE 193.681
graph.density(covid_net_delaunay_igraph) # 0.2061538


# Best KNN ----------------------------------------------------------------
# construct best performing KNN network and GNAR model
# create nb list
opt_knn_net <- knearneigh(x = coord_urbanisation, 
                          k = 21, 
                          longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames())

# create igraph object
opt_knn_net_igraph <- neighborsDataFrame(nb = opt_knn_net) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify()

# create GNAR object 
opt_knn_net_gnar <- opt_knn_net_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame
county_index_opt_knn <- data.frame("CountyName" = opt_knn_net_igraph %>%
                                     V() %>% 
                                     names(), 
                                   "index" = seq(1, 26))

# compute network characteristics 
graph_char_knn <- network_characteristics(opt_knn_net_igraph, 
                                          "KNN")

# for latex 
strCaption <- paste0("Summary for the \\textbf{KNN} network, av. short for 
                     average, s.d. short for standard deviation")
print(xtable(graph_char_knn,
             digits=2,
             caption=strCaption,
             label="tab:summary_knn", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_char_knn)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# visualise network 
plot(st_geometry(ireland_shp),
     border="grey")
plot(opt_knn_net,
     coord_urbanisation,
     add = TRUE,
     pch = 19, cex = 0.6)
text(coord_urbanisation[, 1],
     coord_urbanisation[, 2],
     labels = rownames(coord_urbanisation),
     cex = 0.8, font = 2, pos = 1)

opt_knn_mod <- fit_and_predict(net = opt_knn_net_gnar, 
                               alpha = 4, 
                               beta = c(1, 1, 0, 0), 
                               globalalpha = TRUE, 
                               old = TRUE, 
                               county_index = county_index_opt_knn, 
                               return_model = TRUE)
summary(opt_knn_mod$mod)
BIC(opt_knn_mod)
opt_knn_mod %>% coef()


# compute Moran's I
moran_knn <- moran_I(data = COVID_weekly_data, 
                     nb_list = opt_knn_net)
# visualise Moran's I
ggplot(moran_knn, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"), 
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_knn.pdf", 
       width = 27, height = 14, unit = "cm")


# Best DNN ----------------------------------------------------------------
# construct best performing DNN network and GNAR model 
# create nb list
opt_dnn_net <- dnearneigh(x = coord_urbanisation, 
                          d1 = 0, 
                          d2 = 175,
                          row.names = coord_urbanisation %>% rownames(),
                          longlat = TRUE, 
                          use_s2 = TRUE) 

# create igraph object 
opt_dnn_net_igraph <- neighborsDataFrame(opt_dnn_net) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

# create GNAR object 
opt_dnn_net_gnar <- opt_dnn_net_igraph %>% 
  igraphtoGNAR()

# create ordered county index data frame
county_index_opt_dnn <- data.frame("CountyName" = opt_dnn_net_igraph %>%
                                     V() %>% 
                                     names(), 
                                   "index" = seq(1, 26))

# compute network characteristics
graph_char_dnn <- network_characteristics(opt_dnn_net_igraph, 
                                          "DNN")

# for latex 
strCaption <- paste0("Summary for the \\textbf{DNN} network, av. short for 
                     average, s.d. short for standard deviation")
print(xtable(graph_char_dnn,
             digits=2,
             caption=strCaption,
             label="tab:summary_dnn", 
             align = c("", "l", "|", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_char_dnn)),
                        command = c(paste("\\toprule \n",
                                          " metric & value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Network characteristics for KNN and DNN
graph_knn_dnn <- cbind(graph_char_knn, 
                       graph_char_dnn$DNN, 
                       graph_complete$complete)

strCaption <- "Overview of network characteristics for best performing KNN (k = 21),
DNN (d = 175) and Complete network, including average (av.) degree, density, average (av.) 
shortest path length (SPL), global and average (av.) local clustering (clust.) 
as well as average (av.) betweenness (betw.) and its standard deviation (s.d.)"
print(xtable(graph_knn_dnn,
             digits=2,
             caption=strCaption,
             label="tab:summary_knn_dnn", 
             align = c("", "l", "|", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_knn_dnn)),
                        command = c(paste("\\toprule \n",
                                          " Metric & KNN & DNN & Complete \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# visualise network 
plot(st_geometry(ireland_shp),
     border="grey")
plot(opt_dnn_net,
     coord_urbanisation,
     add = TRUE,
     pch = 19, cex = 0.6)
text(coord_urbanisation[, 1],
     coord_urbanisation[, 2],
     labels = rownames(coord_urbanisation),
     cex = 0.8, font = 2, pos = 1)


# best model 
opt_dnn_mod <- fit_and_predict(net = opt_dnn_net_gnar, 
                               alpha = 4, 
                               beta = c(2, 2, 2, 1), 
                               globalalpha = TRUE, 
                               inverse_distance = TRUE, 
                               county_index = county_index_opt_dnn, 
                               return_model = TRUE)
summary(opt_dnn_mod$mod)
BIC(opt_dnn_mod)


# compute Moran's I
moran_dnn <- moran_I(data = COVID_weekly_data, 
                     nb_list = opt_dnn_net)

# visualise Moran's I
ggplot(moran_dnn, 
       aes(x = dates, 
           y = moran)) +
  geom_line() +
  xlab("Time") +
  ylab("Moran's I") +
  geom_vline(aes(xintercept = as.Date("18.08.2020",
                                      format = "%d.%m.%Y"), 
                 color = "County-specific restrictions")) +
  geom_vline(aes(xintercept = as.Date("26.12.2020", 
                                      format = "%d.%m.%Y"), 
                 color = "Level-5 lockdown")) +
  geom_vline(aes(xintercept = as.Date("26.07.2021",
                                      format = "%d.%m.%Y"), 
                 color = "Indoor dining")) +
  geom_vline(aes(xintercept = as.Date("06.03.2022",
                                      format = "%d.%m.%Y"), 
                 color = "End")) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "None")
ggsave("plots/spatialCor/covid_moran_dnn.pdf", 
       width = 27, height = 14, unit = "cm")

# Best model for each network ---------------------------------------------
# find best model for each sparse and Complete network 
best_model <- compare_df %>% 
  group_by(network) %>% 
  summarise(BIC = min(BIC))

# filter minimum BIC for KNN and DNN network 
best_model_knn_dnn <- compare_knn_dnn %>% 
  group_by(network) %>% 
  summarise(BIC = min(BIC))


best_model_overview <- left_join(best_model, 
                                 compare_df, 
                                 by = c("BIC", "network")) %>% 
  dplyr::select(network, 
                hyperparameters, 
                model, 
                BIC) %>% 
  as.data.frame()

best_model_knn_dnn_overview <- left_join(best_model_knn_dnn, 
                                         compare_knn_dnn, 
                                         by = c("BIC", "network")) %>% 
  dplyr::select(network, 
                hyperparameters, 
                model, 
                BIC) %>% 
  as.data.frame()


# merge data frames with best models for sparse + Complete networks and 
# KNN + DNN networks 
best_model_overview_all <- rbind(best_model_overview, 
                                 best_model_knn_dnn_overview)


# Residual analysis  ------------------------------------------------------

# for each network, fit best performing model and analysis residuals 
# in scatter plot and qqplot 

# Queen 
model_queen <- fit_and_predict(alpha = 4, 
                               beta = c(2, 2, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_queen_gnar, 
                               old = TRUE, 
                               county_index = county_index_queen, 
                               return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_queen <- check_and_plot_predictions(model = model_queen, 
                                                network_name = "queen")
# compute residuals for fitted values and plot in scatter plot
residuals_queen <- check_and_plot_residuals(model = model_queen, 
                                            network_name = "queen")

# Eco hubs
model_eco_hubs <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 1, 0, 0),
                                  globalalpha = TRUE, 
                                  net = covid_net_eco_hubs_gnar, 
                                  old = TRUE, 
                                  numeric_vertices = TRUE, 
                                  county_index = county_index_eco_hubs,
                                  return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_eco_hubs <- check_and_plot_predictions(model = model_eco_hubs, 
                                                   network_name = "eco_hubs")
# compute residuals for fitted values and plot in scatter plot
residuals_eco_hubs <- check_and_plot_residuals(model = model_eco_hubs, 
                                               network_name = "eco_hubs", 
                                               alpha = 5)


# Train 
model_train <- fit_and_predict(alpha = 4, 
                               beta = c(2, 1, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_train_gnar, 
                               old = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train,
                               return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_train <- check_and_plot_predictions(model = model_train, 
                                                network_name = "train")
# compute residuals for fitted values and plot in scatter plot
residuals_train <- check_and_plot_residuals(model = model_train, 
                                            network_name = "train")


# Delaunay 
model_delaunay <- fit_and_predict(alpha = 4, 
                                  beta = c(2, 2, 1, 1),
                                  globalalpha = TRUE, 
                                  net = covid_net_delaunay_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_delaunay, 
                                  return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_delaunay <- check_and_plot_predictions(model = model_delaunay, 
                                                   network_name = "delaunay")
# compute residuals for fitted values and plot in scatter plot
residuals_delaunay <- check_and_plot_residuals(model = model_delaunay, 
                                               network_name = "delaunay")

# Gabriel
model_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(1, 1, 1, 0, 0),
                                 globalalpha = TRUE, 
                                 net = covid_net_gabriel_gnar, 
                                 old = TRUE, 
                                 county_index = county_index_gabriel,
                                 return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_gabriel <- check_and_plot_predictions(model = model_gabriel, 
                                                  network_name = "gabriel")
# compute residuals for fitted values and plot in scatter plot
residuals_gabriel <- check_and_plot_residuals(model = model_gabriel, 
                                              network_name = "gabriel", 
                                              alpha = 5)


# Relative 
model_relative <- fit_and_predict(alpha = 4, 
                                  beta = c(2, 2, 1, 1),
                                  globalalpha = TRUE, 
                                  net = covid_net_relative_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_relative,
                                  return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_relative <- check_and_plot_predictions(model = model_relative, 
                                                   network_name = "relative")
# compute residuals for fitted values and plot in scatter plot
residuals_relative <- check_and_plot_residuals(model = model_relative, 
                                               network_name = "relative")

# SOI
model_soi <- fit_and_predict(alpha = 4, 
                             beta = c(2, 1, 1, 1),
                             globalalpha = TRUE, 
                             net = covid_net_soi_gnar, 
                             old = TRUE, 
                             county_index = county_index_soi, 
                             return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_soi <- check_and_plot_predictions(model = model_soi, 
                                                network_name = "soi")
# compute residuals for fitted values and plot in scatter plot
residuals_soi <- check_and_plot_residuals(model = model_soi, 
                                         network_name = "soi")

# KNN 
# predict last 10 weeks in data set and plot true values and predicted values
predictions_knn <- check_and_plot_predictions(model = opt_knn_mod,
                                              network_name = "KNN")
# compute residuals for fitted values and plot in scatter plot
residuals_knn <- check_and_plot_residuals(model = opt_knn_mod, 
                                          network_name = "KNN")

# DNN
# predict last 10 weeks in data set and plot true values and predicted values
predictions_dnn <- check_and_plot_predictions(model = opt_dnn_mod,
                                              network_name = "DNN")
# compute residuals for fitted values and plot in scatter plot
residuals_dnn <- check_and_plot_residuals(model = opt_dnn_mod, 
                                          network_name = "DNN")

# Complete 
model_complete <- fit_and_predict(alpha = 5,
                                  beta = c(1, 1, 1, 1, 0), 
                                  globalalpha = TRUE, 
                                  net = complete_net_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_complete,
                                  return_model = TRUE)

# predict last 10 weeks in data set and plot true values and predicted values
predictions_complete <- check_and_plot_predictions(model = model_complete, 
                                                   network_name = "complete")
# compute residuals for fitted values and plot in scatter plot
residuals_complete <- check_and_plot_residuals(model = model_complete, 
                                               network_name = "complete", 
                                               alpha = 5)


# add AIC 
best_model_overview_all$AIC <- c(AIC(model_complete),
                                 AIC(model_delaunay), 
                                 AIC(model_eco_hubs), 
                                 AIC(model_gabriel), 
                                 AIC(model_queen), 
                                 AIC(model_train), 
                                 AIC(model_relative), 
                                 AIC(model_soi), 
                                 AIC(opt_dnn_mod), 
                                 AIC(opt_knn_mod)
                                 )

# for latex 
strCaption <- "Overview over best performing \\code{GNARfit()} models for 
each COVID-19 network, SPL weighting best performing for all networks 
excluding the DNN network"
print(xtable(best_model_overview_all[c(2, 4, 8, 7, 1, 10, 9, 5, 3, 6), ],
             digits=2,
             caption=strCaption,
             label="tab:best_models", 
             align = c("", "l", "|", "r", "c", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_model_overview_all)),
                        command = c(paste("\\toprule \n",
                                          " Network & weighing & 
                                          \\code{GNAR} model & BIC & 
                                          AIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# MASE --------------------------------------------------------------------
# compute MASE for the last 10 weeks 

# Queen
mase_queen <- compute_MASE(model = model_queen, 
                           network_name = "Queen")

# Eco hubs
mase_eco_hubs <- compute_MASE(model = model_eco_hubs, 
                              network_name = "Eco hubs")

# Train 
mase_train <- compute_MASE(model = model_train, 
                           network_name = "Train")

# Delaunay 
mase_delaunay <- compute_MASE(model = model_delaunay, 
                              network_name = "Delaunay")

# Gabriel 
mase_gabriel <- compute_MASE(model = model_gabriel, 
                             network_name = "Gabriel")

# Relative 
mase_relative <- compute_MASE(model = model_relative, 
                              network_name = "Relative")

# Soi
mase_soi <- compute_MASE(model = model_soi, 
                         network_name = "SOI")

# KNN
mase_knn <- compute_MASE(model = opt_knn_mod, 
                         network_name = "KNN")

# DNN
mase_dnn <- compute_MASE(model = opt_dnn_mod, 
                         network_name = "DNN")

# Complete 
mase_complete <- compute_MASE(model = model_complete, 
                              network_name = "Complete")

# summarise all MASE values for each network in a data frame  
mase_overview <- rbind.data.frame(mase_queen, 
                                  mase_eco_hubs, 
                                  mase_train, 
                                  mase_delaunay, 
                                  mase_gabriel, 
                                  mase_relative, 
                                  mase_soi, 
                                  mase_knn, 
                                  mase_dnn, 
                                  mase_complete, 
                                  mase_arima) %>% 
  na.omit()




# plot MASE for counties: Dublin, Wicklow, Kerry, Donegal
ggplot(mase_overview, 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("ARIMA" = "grey", 
                                "KNN" = "#F8766D", 
                                "DNN" = "#D89000", 
                                "Complete" = "#A3A500", 
                                "Queen" = "#E76BF3", 
                                "Eco hubs" = "#00BF7D", 
                                "Gabriel" = "#00BFC4", 
                                "Relative" = "#00B0F6", 
                                "SOI" = "#9590FF", 
                                "Delaunay" = "#39B600", 
                                "Train" = "#FF62BC"), 
                     name = "Network")
ggsave("plots/Prediction/mase.pdf", 
       width = 26, height = 13, units = "cm")

# plot MASE for the Delaunay, Gabriel, Relative and SOI as well as Rail-way 
# based network 
ggplot(mase_overview %>% filter(type %in% c("Delaunay", 
                                            "Gabriel", 
                                            "Relative", 
                                            "SOI", 
                                            "Train", 
                                            "ARIMA")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  ylim(0, 8) +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("ARIMA" = "grey", 
                                "Gabriel" = "#00BFC4", 
                                "Relative" = "#00B0F6", 
                                "SOI" = "#9590FF", 
                                "Delaunay" = "#39B600", 
                                "Train" = "#FF62BC"), 
                     name = "Network")
ggsave("plots/Prediction/mase_delaunay_etc_zoom.pdf", 
       width = 26, height = 13, units = "cm")

# plot MASE for the KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_overview %>% filter(type %in% c("DNN", 
                                            "KNN", 
                                            "Queen", 
                                            "Eco hubs", 
                                            "Complete", 
                                            "ARIMA")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  ylim(0, 8) +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("ARIMA" = "grey", 
                                "KNN" = "#F8766D", 
                                "DNN" = "#D89000", 
                                "Complete" = "#A3A500", 
                                "Queen" = "#E76BF3", 
                                "Eco hubs" = "#00BF7D"), 
                     name = "Network")
ggsave("plots/Prediction/mase_knn_etc_zoom.pdf", 
       width = 26, height = 13, units = "cm")

# Optimise beta order -----------------------------------------------------
# Queen
# compute upper limit for neighbourhood stage 
upper_limit_queen <- distance_table(covid_net_queen_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_queen <- optimise_beta(GNAR_object  = covid_net_queen_gnar, 
                                                          alpha = 4, 
                                                          beta = c(2, 2, 1, 1), 
                                                          upper_limit_beta = upper_limit_queen, 
                                                          old_best = best_model_overview_all[best_model_overview_all$network == "Queen", ]$BIC, 
                                                          network_name = "queen")

# Economic hub 
# compute upper limit for neighbourhood stage
upper_limit_eco_hub <- distance_table(covid_net_eco_hubs_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_eco_hubs <- optimise_beta(GNAR_object  = covid_net_eco_hubs_gnar, 
                                                          alpha = 5, 
                                                          beta = c(1, 1, 1, 0, 0), 
                                                          upper_limit_beta = upper_limit_eco_hub, 
                                                          old_best = best_model_overview_all[best_model_overview_all$network == "Eco. hub", ]$BIC, 
                                                          network_name = "eco_hubs")

# Train 
# compute upper limit for neighbourhood stage 
upper_limit_train <- distance_table(covid_net_train_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_train <- optimise_beta(GNAR_object  = covid_net_train_gnar, 
                                                             alpha = 4, 
                                                             beta = c(2, 1, 1, 1), 
                                                             upper_limit_beta = upper_limit_train, 
                                                             old_best = best_model_overview_all[best_model_overview_all$network == "Rail", ]$BIC, 
                                                             network_name = "train")

# Delaunay
# compute upper limit for neighbourhood stage 
upper_limit_delaunay <- distance_table(covid_net_delaunay_igraph)$res %>% length()

# optimise beta order components
neighbourhood_size_delaunay <- optimise_beta(GNAR_object  = covid_net_delaunay_gnar, 
                                                             alpha = 4, 
                                                             beta = c(2, 2, 1, 1), 
                                                             upper_limit_beta = upper_limit_delaunay, 
                                                             old_best = best_model_overview_all[best_model_overview_all$network == "Delaunay", ]$BIC, 
                                                             network_name = "delaunay")

# Gabriel
# compute upper limit for neighbourhood stage 
upper_limit_gabriel <- distance_table(covid_net_gabriel_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_gabriel <- optimise_beta(GNAR_object  = covid_net_gabriel_gnar, 
                                                             alpha = 5, 
                                                             beta = c(1, 1, 1, 0, 0), 
                                                             upper_limit_beta = upper_limit_gabriel, 
                                                             old_best = best_model_overview_all[best_model_overview_all$network == "Gabriel", ]$BIC, 
                                                             network_name = "gabriel")

# Relative
# compute upper limit for neighbourhood stage 
upper_limit_relative <- distance_table(covid_net_relative_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_relative <- optimise_beta(GNAR_object  = covid_net_relative_gnar, 
                                                             alpha = 4, 
                                                             beta = c(2, 2, 1, 1), 
                                                             upper_limit_beta = upper_limit_relative, 
                                                             old_best = best_model_overview_all[best_model_overview_all$network == "Relative", ]$BIC, 
                                                             network_name = "relative")

# SOI
# compute upper limit for neighbourhood stage 
upper_limit_soi <- distance_table(covid_net_soi_igraph)$res %>% length()

# optimise beta order components
neighbourhood_size_soi <- optimise_beta(GNAR_object  = covid_net_soi_gnar, 
                                                             alpha = 4, 
                                                             beta = c(2, 1, 1, 1), 
                                                             upper_limit_beta = upper_limit_soi, 
                                                             old_best = best_model_overview_all[best_model_overview_all$network == "SOI", ]$BIC, 
                                                             network_name = "soi")

# KNN
# compute upper limit for neighbourhood stage 
upper_limit_knn <- distance_table(opt_knn_net_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_knn <- optimise_beta(GNAR_object  = opt_knn_net_gnar, 
                                                             alpha = 4, 
                                                             beta = c(1, 1, 0, 0), 
                                                             upper_limit_beta = upper_limit_knn, 
                                                             old_best = best_model_overview_all[best_model_overview_all$network == "KNN", ]$BIC, 
                                                             network_name = "knn")

# DNN
# compute upper limit for neighbourhood stage 
upper_limit_dnn <- distance_table(opt_dnn_net_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_dnn <- optimise_beta(GNAR_object  = opt_dnn_net_gnar, 
                                                        alpha = 4, 
                                                        beta = c(2, 2, 2, 1), 
                                                        upper_limit_beta = upper_limit_dnn,
                                                        old = FALSE, 
                                                        county_index = county_index_opt_dnn, 
                                                        old_best = best_model_overview_all[best_model_overview_all$network == "DNN", ]$BIC, 
                                                        network_name = "dnn")

# Complete
# compute upper limit for neighbourhood stage 
upper_limit_complete <- distance_table(complete_net_igraph)$res %>% length()

# optimise beta order components 
neighbourhood_size_complete <- optimise_beta(GNAR_object  = complete_net_gnar, 
                                                        alpha = 5, 
                                                        beta = c(1, 1, 1, 1, 0), 
                                                        upper_limit_beta = upper_limit_complete, 
                                                        old_best = best_model_overview_all[best_model_overview_all$network == "Complete", ]$BIC, 
                                                        network_name = "complete")

# create overview data frame with optimised beta order for each network 
new_best_parameters <- data.frame(network = c("Delaunay", 
                                              "Gabriel", 
                                              "SOI",
                                              "Relative",
                                              "Complete",
                                              "KNN", 
                                              "DNN", 
                                              "Queen",
                                              "Eco. hub", 
                                              "Rail"
                                              ), 
                                  param = c(neighbourhood_size_delaunay$best_beta,
                                            neighbourhood_size_gabriel$best_beta, 
                                            neighbourhood_size_soi$best_beta, 
                                            neighbourhood_size_relative$best_beta, 
                                            neighbourhood_size_complete$best_beta, 
                                            neighbourhood_size_knn$best_beta, 
                                            neighbourhood_size_dnn$best_beta, 
                                            neighbourhood_size_queen$best_beta, 
                                            neighbourhood_size_eco_hubs$best_beta, 
                                            neighbourhood_size_train$best_beta 
                                            
                                  ))

# compare old vs. new model 
compare_parameters <- left_join(best_model_overview_all %>% 
                                  dplyr::select(network, model), 
                                new_best_parameters, 
                                by = "network")

# for latex 
strCaption <- "Comparison of initial best performing model and new $\\beta$ 
order for each COVID-19 network after optimising the components of the 
$\\beta$ order in isolation"
print(xtable(compare_parameters,
             digits=2,
             caption=strCaption,
             label="tab:best_adapted_models", 
             align = c("", "l", "|", "c", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(compare_parameters)),
                        command = c(paste("\\toprule \n",
                                          " Network & \\code{GNAR} model & $\\beta$ order \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# Adapted models ----------------------------------------------------------
# fit GNAR models with optimised beta order
# no change for KNN network 
new_model_queen <- fit_and_predict(alpha = 4, 
                                   beta = c(3, 2, 1, 1),
                                   globalalpha = TRUE, 
                                   net = covid_net_queen_gnar, 
                                   old = TRUE,  
                                   return_model = TRUE)

new_model_eco_hubs <- fit_and_predict(alpha = 5, 
                                     beta = c(3, 1, 1, 0, 0),
                                     globalalpha = TRUE, 
                                     net = covid_net_eco_hubs_gnar, 
                                     old = TRUE,  
                                     return_model = TRUE)

new_model_train <- fit_and_predict(alpha = 4, 
                                   beta = c(4, 1, 1, 1),
                                   globalalpha = TRUE, 
                                   net = covid_net_train_gnar, 
                                   old = TRUE,  
                                   return_model = TRUE)

new_model_delaunay <- fit_and_predict(alpha = 4, 
                                      beta = c(4, 2, 1, 1),
                                      globalalpha = TRUE, 
                                      net = covid_net_delaunay_gnar, 
                                      old = TRUE,  
                                      return_model = TRUE)

new_model_gabriel <- fit_and_predict(alpha = 5, 
                                     beta = c(3, 1, 1, 0, 0),
                                     globalalpha = TRUE, 
                                     net = covid_net_gabriel_gnar, 
                                     old = TRUE,  
                                     return_model = TRUE)

new_model_relative <- fit_and_predict(alpha = 4, 
                                      beta = c(5, 2, 1, 1),
                                      globalalpha = TRUE, 
                                      net = covid_net_relative_gnar, 
                                      old = TRUE,  
                                      return_model = TRUE)

new_model_soi <- fit_and_predict(alpha = 4, 
                                 beta = c(4, 1, 1, 1),
                                 globalalpha = TRUE, 
                                 net = covid_net_soi_gnar, 
                                 old = TRUE,  
                                 return_model = TRUE)

new_model_dnn <- fit_and_predict(alpha = 4, 
                                 beta = c(2, 2, 2, 0),
                                 globalalpha = TRUE, 
                                 net = opt_dnn_net_gnar, 
                                 old = FALSE,
                                 inverse_distance = TRUE, 
                                 county_index = county_index_opt_dnn,  
                                 return_model = TRUE)

new_model_complete <- fit_and_predict(alpha = 5, 
                                      beta = c(1, 1, 0, 1, 0),
                                      globalalpha = TRUE, 
                                      net = complete_net_gnar, 
                                      old = TRUE,  
                                      return_model = TRUE)


# Compare MASE ------------------------------------------------------------
# compute MASE for GNAR models with optimised beta order and compare to 
# MASE for the initial best performing GNAR model for each network 

# Queen
mase_new_queen <- compute_MASE(model = new_model_queen, 
                               network_name = "Queen")

# Eco hubs
mase_new_eco_hubs <- compute_MASE(model = new_model_eco_hubs, 
                                  network_name = "Eco hubs")

# Train 
mase_new_train <- compute_MASE(model = new_model_train, 
                               network_name = "Train")

# Delaunay 
mase_new_delaunay <- compute_MASE(model = new_model_delaunay, 
                                  network_name = "Delaunay")

# Gabriel 
mase_new_gabriel <- compute_MASE(model = new_model_gabriel, 
                                 network_name = "Gabriel")

# Relative 
mase_new_relative <- compute_MASE(model = new_model_relative, 
                                  network_name = "Relative")

# Soi
mase_new_soi <- compute_MASE(model = new_model_soi, 
                             network_name = "SOI")

# KNN: no change in parameters 

# DNN
mase_new_dnn <- compute_MASE(model = new_model_dnn, 
                             network_name = "DNN")

# Complete 
mase_new_complete <- compute_MASE(model = new_model_complete, 
                                  network_name = "Complete")


# plot MASE for initial best performing GNAR model and the GNAR model with 
# optimised beta order 
compare_MASE(mase_queen, 
             mase_new_queen, 
             network_name = "queen", 
             old_model = "2211", 
             new_model = "3211")

compare_MASE(mase_eco_hubs, 
             mase_new_eco_hubs, 
             network_name = "eco_hub", 
             old_model = "11100", 
             new_model = "31100")

compare_MASE(mase_train, 
             mase_new_train, 
             network_name = "train",
             old_model = "2111", 
             new_model = "4111")

compare_MASE(mase_delaunay, 
             mase_new_delaunay, 
             network_name = "delaunay", 
             old_model = "2211", 
             new_model = "4211")

compare_MASE(mase_gabriel, 
             mase_new_gabriel,
             network_name = "gabriel", 
             old_model = "11100", 
             new_model = "31100")

compare_MASE(mase_relative,
             mase_new_relative, 
             network_name = "relative", 
             old_model = "2211", 
             new_model = "5211")

compare_MASE(mase_soi, 
             mase_new_soi, 
             network_name = "soi", 
             old_model = "2111", 
             new_model = "4111")

compare_MASE(mase_dnn, 
             mase_new_dnn, 
             network_name = "dnn", 
             old_model = "2221", 
             new_model = "2220")

compare_MASE(mase_complete, 
             mase_new_complete, 
             network_name = "complete", 
             old_model = "11110", 
             new_model = "11010")


# Density vs.  complexity -------------------------------------------------
# plot the density vs. the complexity of model, 
# where complexity is defined as the sum of alpha and beta order
density_vs_beta <- data.frame("density" = c(complete_net_igraph %>% graph.density(), 
                                            covid_net_delaunay_igraph %>% graph.density(), 
                                            covid_net_eco_hubs_igraph %>% graph.density(), 
                                            covid_net_gabriel_igraph %>% graph.density(), 
                                            covid_net_queen_igraph %>% graph.density(), 
                                            covid_net_train_igraph %>% graph.density(), 
                                            covid_net_relative_igraph %>% graph.density(), 
                                            covid_net_soi_igraph %>% graph.density(), 
                                            opt_knn_net_igraph %>% graph.density(), 
                                            opt_dnn_net_igraph %>% graph.density()), 
                              "sum_param" = c(c(5, 1, 1, 1, 1) %>% sum(), 
                                              c(4, 2, 2, 1, 1) %>% sum(), 
                                              c(5, 1, 1, 1) %>% sum(), 
                                              c(5, 1, 1, 1) %>% sum(), 
                                              c(4, 2, 2, 1, 1) %>% sum(), 
                                              c(4, 2, 1, 1, 1) %>% sum(), 
                                              c(4, 2, 2, 1, 1) %>% sum(), 
                                              c(4, 2, 1, 1, 1) %>% sum(), 
                                              c(4, 1, 1, 0, 0) %>% sum(), 
                                              c(4, 2, 2, 2, 1) %>% sum()), 
                              "sum_param_opt" = c(c(5, 1, 1, 1) %>% sum(), 
                                                  c(4, 4, 2, 1, 1) %>% sum(),
                                                  c(5, 3, 1, 1) %>% sum(), 
                                                  c(5, 3, 1, 1) %>% sum(),
                                                  c(4, 3, 2, 1, 1) %>% sum(),
                                                  c(4, 4, 1, 1, 1) %>% sum(), 
                                                  c(4, 5, 2, 1, 1) %>% sum(), 
                                                  c(4, 4, 1, 1, 1) %>% sum(), 
                                                  c(4, 1, 1) %>% sum(), 
                                                  c(4, 2, 2, 2) %>% sum())
)

ggplot(density_vs_beta, 
       aes(x = density, 
           y = sum_param, 
           shape = "best fit", 
           color = "best fit")) +
  geom_point(size = 2) +
  geom_line(stat = "smooth", 
            linetype = "dashed", 
            method = "lm", 
            se = FALSE, 
            alpha = 0.7) +
  geom_point(aes(x = density, 
                 y = sum_param_opt, 
                 shape = "optimised", 
                 color = "optimised"), 
             size = 2) +
  geom_line(aes(x = density, 
                y = sum_param_opt, 
                shape = "optimised", 
                color = "optimised"), 
            stat = "smooth", 
            linetype = "dashed", 
            method = "lm", 
            se = FALSE, 
            alpha = 0.7) +
  xlab("Network density") +
  ylab("Sum of parameter") +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(title = ""),
         shape = guide_legend(title = "")) +
  theme(legend.position = "bottom")
ggsave("plots/modelfit/density_vs_sum_params.pdf", 
       width = 14, height = 14, units = "cm")


# Save objects -----------------------------------------------------------
# save GNAR objects for every network
save(list = c("covid_net_queen_gnar",
              "covid_net_eco_hubs_gnar",
              "covid_net_train_gnar",
              "covid_net_delaunay_gnar",
              "covid_net_gabriel_gnar",
              "covid_net_relative_gnar",
              "covid_net_soi_gnar",
              "opt_knn_net_gnar",
              "opt_dnn_net_gnar",
              "complete_net_gnar"),
     file = "data/RObjects/GNAR.RData")

# save igraph objects for every network
save(list = c("covid_net_queen_igraph",
              "covid_net_eco_hubs_igraph",
              "covid_net_train_igraph",
              "covid_net_delaunay_igraph",
              "covid_net_gabriel_igraph",
              "covid_net_relative_igraph",
              "covid_net_soi_igraph",
              "opt_knn_net_igraph",
              "opt_dnn_net_igraph",
              "complete_net_igraph"),
     file = "data/RObjects/igraph.RData")

# save county indices for every network
save(list = c("county_index_queen",
              "county_index_eco_hubs",
              "county_index_train",
              "county_index_delaunay",
              "county_index_gabriel",
              "county_index_relative",
              "county_index_soi",
              "county_index_opt_knn",
              "county_index_opt_dnn",
              "county_index_complete"),
     file = "data/RObjects/county_index.RData")

# save vector with population sizes
save(population_weight,
     file = "data/RObjects/population_weight.RData")

# save matrix of distances between counties
save(dist_urbanisation,
     file = "data/RObjects/distance_urbanisation.RData")

# save coordinates to represent counties 
# (centroids for rural, county towns for urban)
save(coord_urbanisation,
     file = "data/RObjects/coord_urbanisation.RData")
