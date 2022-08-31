## ---- Install packages -------------------------------------------------------

# install.packages("igraph)
# install.packages("GNAR)
# install.packages("sf")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# install.packages("rnaturalearthhires")
# install.packages("mapdata")
# install.packages("rgdal")
# install.packages("ggmap")
# install.packages("rgeos")
# install.packages("broom")
# install.packages("TTR")
# install.packages("sf")
# install.packages("raster")
# install.packages("terra")
# install.packages("leaflet")


## ---- Load libraries ---------------------------------------------------------
library(igraph)
library(GNAR)
library(tidyverse)
library(magrittr) # for pipping 
library(xtable) # for tables 

# for MAPS: 
library(sp)
library(sf) 
library(raster) 
library(leaflet)
library(mapview) # for saving

# Set ggplot theme 
theme_set(theme_bw(base_size = 16))

# Load functions 
source(file = "functions.R")

## ---- Data -------------------------------------------------------------------
load("data/Mumps/mumps.RData")

mumps %>% class() # matrix 
mumps %>% nrow() # duration: 52 weeks 
mumps %>% ncol() # 47 counties in England, Wales is considered as one

## ---- Population correction for mumps count ----------------------------------
# list of population sizes for each county to perform population correction
county_inhabitants <- list("Cornwall and Isles of Scilly" = 520000, 
                           "Devon" = 1110000, 
                           "Somerset" = 886000, 
                           "Dorset" = 699000, 
                           "Isle of Wight" = 137800, 
                           "West Sussex" = 770600, 
                           "East Sussex" = 753000, 
                           "Kent" = 1792000, 
                           "Surrey" = 1070900, 
                           "Hampshire" = 1680000, 
                           "Wiltshire" = 631000,
                           "Avon" = 663000, 
                           "Oxfordshire" = 627500, 
                           "Berkshire" = 809000, 
                           "London" = 7484900, 
                           "Hertfordshire" = 1051500, 
                           "Essex" = 1654000, 
                           "Buckinghamshire" = 483700, 
                           "Suffolk" = 694600, 
                           "Gloucestershire" = 576300,
                           "Warwickshire" = 521500, 
                           "Northamptonshire" = 655400, 
                           "Bedfordshire" = 583000, 
                           "Cambridgeshire" = 744000, 
                           "Norfolk" = 824400, 
                           "Wales" = 2950100, 
                           "Hereford and Worcester" = 726300, 
                           "West Midlands" = 2593600, 
                           "Staffordshire" = 1059000, 
                           "Shropshire" = 447000, 
                           "Leicestershire" = 902000, 
                           "Nottinghamshire" = 1051000, 
                           "Derbyshire" = 985000, 
                           "Lincolnshire" = 996000, 
                           "South Yorkshire" = 1288200, 
                           "Cheshire" = 996000,
                           "Merseyside" = 1358500, 
                           "Greater Manchester" = 2543900, 
                           "West Yorkshire" = 2148600,  
                           "Humberside" = 586000, 
                           "North Yorkshire" = 772000, 
                           "Lancashire" = 1444000, 
                           "Cumbria" = 496700, 
                           "Northumberland" = 309400, 
                           "Tyne and Wear" = 1088400, 
                           "Durham" = 594000, 
                           "Cleveland" = 556000
)

# convert matrix to data frame 
mumps_corr <- mumps %>% as.data.frame()

# compute population-corrected mumps incidence for 1000.000 inhabitants 
for (county in colnames(mumps_corr)) {
  mumps_corr[[county]] <- (mumps_corr[[county]] * 1000000 /
                             county_inhabitants[[county]]) %>% 
    round()
}

# find max mumps incidence for each county
max_incidence <- colMax(mumps_corr)
which.max(max_incidence)
max(max_incidence)
# overall max Warwickshire with incidence 228 in week 43

# copy for editing in the following sections  
mumps_df <- mumps_corr

## ---- Data pre-processing and central towns ----------------------------------
# list to allot each county its county town 
central_towns <- list("Cornwall and Isles of Scilly" = "Truro", 
                      "Devon" = "Exeter", 
                      "Somerset" = "Taunton", 
                      "Dorset" = "Dorchester", 
                      "Isle of Wight" = "Newport", 
                      "West Sussex" = "Chichester", 
                      "East Sussex" = "Lewes", 
                      "Kent" = "Maidstone", 
                      "Surrey" = "Guildford", 
                      "Hampshire" = "Winchester", 
                      "Wiltshire" = "Devizes",
                      "Avon" = "Bristol", 
                      "Oxfordshire" = "Oxford", 
                      "Berkshire" = "Reading", 
                      "London" = "London", 
                      "Hertfordshire" = "Hertford", 
                      "Essex" = "Chelmsford", 
                      "Buckinghamshire" = "Aylesbury", 
                      "Suffolk" = "Ipswich", 
                      "Gloucestershire" = "Gloucester",
                      "Warwickshire" = "Warwick", 
                      "Northamptonshire" = "Northampton", 
                      "Bedfordshire" = "Bedford", 
                      "Cambridgeshire" = "Cambridge", 
                      "Norfolk" = "Norwich", 
                      "Wales" = "Rhayader", 
                      "Hereford and Worcester" = "Worcester", 
                      "West Midlands" = "Birmingham", 
                      "Staffordshire" = "Stafford", 
                      "Shropshire" = "Shrewsbury", 
                      "Leicestershire" = "Leicester", 
                      "Nottinghamshire" = "Nottingham", 
                      "Derbyshire" = "Derby", 
                      "Lincolnshire" = "Lincoln", 
                      "South Yorkshire" = "Sheffield", 
                      "Cheshire" = "Chester",
                      "Merseyside" = "Liverpool", 
                      "Greater Manchester" = "Manchester", 
                      "West Yorkshire" = "Leeds",  
                      "Humberside" = "Kingston Upon Hull", 
                      "North Yorkshire" = "York", 
                      "Lancashire" = "Lancaster", 
                      "Cumbria" = "Carlisle", 
                      "Northumberland" = "Morpeth", 
                      "Tyne and Wear" = "Newcastle", 
                      "Durham" = "Durham", 
                      "Cleveland" = "Middlesborough"
                      )

# check if all counties are attributed a central town 
central_towns[!(names(central_towns) %in% colnames(mumps))]


# create a data frame with central towns as column names instead of counties
mumps_central_towns <- mumps_df %>% 
  as.data.frame()

# save previous order to ensure same ordering
ordered_counties <- colnames(mumps_central_towns)

# assign new column names
colnames(mumps_central_towns) <- unlist(unname(central_towns[ordered_counties]))

# create time column with week number 
mumps_central_towns <- mumps_central_towns %>% 
  mutate(week = seq(1, nrow(mumps_central_towns)))


## ---- Maps -------------------------------------------------------------------
# visualise the data and network 

# download UK data level 2 from GADM database 
uk <- getData('GADM', 
              country='GBR', 
              level = 2)

# filter counties in Wales and England
england_wales <- uk[uk$NAME_1 %in% c("England", "Wales"), ]

# plot maps of Wales and English counties
pal <- colorFactor("Reds", uk$NAME_1)

leaflet(england_wales) %>%
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5,
              fillColor = ~pal(NAME_1),
              highlightOptions = highlightOptions(color = "white", weight =2,
                                                  bringToFront = TRUE),
              label = ~NAME_2)


# create maps visualising the mumps incidence for week 1, 25, 52 and save 
mumps_week_1 <- map_mumps(desired_week = 1)
mapshot(mumps_week_1, 
        file = "plots/mapsreplicate/mumps_map_week_1.jpeg")

mumps_week_25 <- map_mumps(desired_week = 25)
mapshot(mumps_week_25, 
        file = "plots/mapsreplicate/mumps_map_week_25.jpeg")

mumps_week_52 <- map_mumps(desired_week = 52)
mapshot(mumps_week_52, 
        file = "plots/mapsreplicate/mumps_map_week_52.jpeg")


## ---- Plot data --------------------------------------------------------------
# plot mumps incidence across all counties and time 
mumps_central_towns_long <- mumps_central_towns %>% 
  gather("Town", "mumps",  -"week") 

ggplot(data = mumps_central_towns_long, 
       aes(x = week, 
           y = mumps, 
           color = Town)) +
  geom_line() +
  ggtitle("Population-corrected mumps cases") +
  xlab("time") +
  ylab("Mumps cases") +
  geom_vline(aes(xintercept = 20), 
             linetype = "dashed", 
             color = "grey")

# filter out the 10 towns with lowest case number in week 1
town_lowest_cases <- mumps_central_towns_long %>% 
  filter(week == 1) %>% 
  arrange(mumps) %>% 
  pull(Town) %>% 
  head(10)

# plot mumps incidence for 10 towns with lowest incidence in week 1
ggplot(data = mumps_central_towns_long %>% 
         filter(Town %in% town_lowest_cases), 
       aes(x = week, 
           y = mumps, 
           color = Town)) +
  geom_line() +
  ggtitle("Population-corrected mumps cases") +
  xlab("time") +
  ylab("Mumps cases")


## ---- Adjacency matrix -------------------------------------------------------
# manually create mumps network from Knight, Nunes and Nason 2016
adjacency_list <- list(
  "Truro" = "Exeter", 
  "Exeter" = c("Truro", "Taunton", "Dorchester"), 
  "Dorchester" = c("Exeter", "Taunton", "Devizes", "Winchester"), 
  "Newport" = "Winchester", 
  "Winchester" = c("Dorchester", "Devizes", "Guildford", "Reading", 
                   "Chichester", "Newport"), 
  "Chichester" = c("Winchester", "Guildford", "Lewes"), 
  "Lewes" = c("Chichester", "Guildford", "Maidstone"), 
  "Maidstone" = c("Lewes", "Guildford", "London", "Chelmsford"), 
  "Taunton" = c("Devizes", "Exeter", "Dorchester", "Bristol"), 
  "Devizes" = c("Taunton", "Winchester", "Dorchester", "Bristol", 
                "Oxford", "Reading", "Gloucester"), 
  "Guildford" = c("Reading", "Winchester", "Chichester", "Lewes", 
                  "Maidstone", "London"), 
  "Reading" = c("Devizes", "Oxford", "London", "Aylesbury", "Guildford"), 
  "London"= c("Guildford", "Reading", "Maidstone", "Chelmsford", "Hertford"), 
  "Chelmsford" = c("London", "Hertford", "Cambridge", "Ipswich"), 
  "Ipswich" = c("Chelmsford", "Cambridge", "Norwich"), 
  "Norwich" = c("Ipswich", "Cambridge", "Lincoln"), 
  "Cambridge" = c("Lincoln", "Norwich", "Chelmsford", "Northampton", 
                  "Bedford", "Hertford", "Ipswich"), 
  "Bedford" = c("Aylesbury", "Hertford", "Cambridge", "Northampton"), 
  "Northampton" = c("Cambridge", "Bedford", "Warwick", "Oxford", "Leicester", 
                    "Lincoln"), 
  "Oxford" = c("Gloucester", "Warwick", "Worcester", "Northampton", 
               "Aylesbury","Reading", "Winchester", "Devizes"), 
  "Bristol" = c("Taunton", "Gloucester", "Devizes"), 
  "Hertford" = c("London", "Chelmsford", "Aylesbury", "Cambridge", "Bedford"), 
  "Aylesbury" = c("Hertford", "Bedford", "Oxford", "Northampton", "Reading"),
  "Gloucester" = c("Bristol", "Devizes", "Oxford", "Warwick", "Worcester", 
                   "Rhayader"),
  "Rhayader" = c("Chester", "Shrewsbury", "Worcester", "Gloucester"), 
  "Worcester" = c("Rhayader", "Shrewsbury", "Stafford", "Birmingham", 
                  "Warwick", "Oxford", "Gloucester"), 
  "Birmingham" = c("Worcester", "Warwick", "Stafford"), 
  "Warwick" = c("Birmingham", "Stafford", "Worcester", "Leicester", 
                "Northampton", "Oxford", "Gloucester"), 
  "Stafford" = c("Worcester", "Birmingham", "Warwick", "Leicester", 
                 "Derby", "Chester", "Shrewsbury"), 
  "Shrewsbury" = c("Rhayader", "Worcester", "Stafford", "Chester"), 
  "Leicester" = c("Stafford", "Nottingham", "Derby", "Lincoln", "Northampton", 
                  "Warwick"), 
  "Nottingham" = c("Leicester", "Derby", "Kingston Upon Hull", "Lincoln", 
                   "Sheffield"), 
  "Derby" = c("Chester", "Manchester", "Sheffield", "Leeds"), 
  "Lincoln" = c("Leicester", "Nottingham", "Northampton", "Norwich", 
                "Cambridge", "Kingston Upon Hull"), 
  "Sheffield" = c("Leeds", "Derby", "Nottingham", "Kingston Upon Hull", "York"), 
  "Leeds" = c("Sheffield", "Derby", "Manchester", "Lancaster", "York"), 
  "Kingston Upon Hull" = c("York", "Sheffield", "Nottingham", "Lincoln"), 
  "York" = c("Lancaster", "Carlisle", "Durham", "Middlesborough", 
             "Kingston Upon Hull", "Sheffield", "Leeds"), 
  "Manchester" = c("Leeds", "Lancaster", "Liverpool", "Chester", "Derby"), 
  "Chester" = c("Liverpool", "Manchester", "Derby", "Stafford", "Shrewsbury", 
                "Rhayader"), 
  "Liverpool" = c("Manchester", "Chester", "Lancaster"), 
  "Lancaster" = c("Liverpool", "Carlisle", "Manchester", "Leeds", "York"), 
  "Carlisle"= c("Morpeth", "Durham", "York", "Lancaster"), 
  "Morpeth" = c("Carlisle", "Newcastle", "Durham"), 
  "Newcastle" = c("Morpeth", "Durham"), 
  "Durham" = c("Morpeth", "Carlisle", "Middlesborough", "York"), 
  "Middlesborough" = c("Durham", "York")
)

# create a data frame with 0
adjacency_df <- as.data.frame(matrix(rep(0, 47 * 47), 
                                         nrow = 47))
# assign central towns as column and row nanes 
colnames(adjacency_df) <- names(adjacency_list)
rownames(adjacency_df) <- names(adjacency_list)

# insert 1 to indicate edges between central towns 
for (name in names(adjacency_list)) {
  adjacency_df[name, 
               unlist(adjacency_list[name])] = rep(1, 
                                                   length(adjacency_list[name]))
}

# transform adjacency data frame into a matrix
adjacency_matrix <- adjacency_df %>% as.matrix()

# create igraph object from adjacency matrix  
mumps_graph <- adjacency_matrix %>% 
  igraph::graph_from_adjacency_matrix(mode = "undirected")


# substitute county towns with county index for nb-list generation
town_to_index <- data.frame("town" = rownames(adjacency_df), 
                            "index" = seq(1, 47))

rownames(adjacency_df) <- seq(1, 47)
colnames(adjacency_df) <- seq(1, 47)

# create igraph object with numeric vertex names to extract nb list 
mumps_numeric_graph <- adjacency_df %>% as.matrix() %>% 
  igraph::graph_from_adjacency_matrix(mode = "undirected")

# generate nb list 
mumps_graph_nb <- mumps_numeric_graph %>% igraph2nb()

# GNARnet object
mumps_gnar <- igraphtoGNAR(mumps_graph)

## ---- Network characteristics ------------------------------------------------
gorder(mumps_graph) # 47 vertices in graph 
gsize(mumps_graph) # 111 edges 

graph_char <- network_characteristics(mumps_graph, 
                                      "Mumps")

# for latex 
strCaption <- "Summary table with network characteristics 
for \\textbf{Mumps} network"
print(xtable(graph_char,
             digits=2,
             caption=strCaption,
             label="tab:summary_mumps_network", 
             align = c("", "l", "|", "c")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(graph_char)),
                        command = c(paste("\\toprule \n",
                                          " metric & value 
                                          \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# plot mumps network 
pdf("plots/mapsreplicate/mumps_graph.pdf")
plot(mumps_graph, 
     layout=layout_nicely, 
     vertex.size = degree(mumps_graph), 
     vertex.label.font = 1, 
     rescale = TRUE, 
     vertex.label.color = "#242526", 
     vertex.color = "#FF7F7F")
dev.off()


## ---- NARIMA -----------------------------------------------------------------
# transform mumps matrix into a ts object 
mumps_ts <- as.matrix(mumps_central_towns %>%
                        dplyr::select(-week))

# log transformation to ensure variance homogeneity 
mumps_log_ts <- log(mumps_ts + 1)

## ---- NAR(1, 1) --------------------------------------------------------------
nar_1_1 <- GNARfit(vts = mumps_ts, 
                   net = mumps_gnar, 
                   alphaOrder = 1, 
                   betaOrder = 1, 
                   globalalpha = TRUE)

summary(nar_1_1)

# PREDICTIONS
mumps_pred_nar_1_1 <- fitted(nar_1_1) %>% 
  as.data.frame()

# assign central towns as column names  
colnames(mumps_pred_nar_1_1) <- colnames(mumps_central_towns %>%
                                           dplyr::select(-week))
# generate time column 
mumps_pred_nar_1_1 <- mumps_pred_nar_1_1 %>% 
  mutate(week = seq(2, nrow(mumps_pred_nar_1_1) + 1))

# visualise predicted values for county Oxford 
ggplot(data = mumps_central_towns, 
       aes(x = week, 
           y = Oxford)) +
  geom_line() +
  geom_line(data = mumps_pred_nar_1_1, 
            aes(x = week, 
                y = Oxford), 
            color = "red")
ggsave("plots/modelreplicate/mumps_nar_1_1_oxford.pdf")

# RESIDUALS
mumps_res_nar_1_1 <- residuals(nar_1_1) %>% 
  as.data.frame()
# assign ccentrral towns as column names 
colnames(mumps_res_nar_1_1) <- colnames(mumps_central_towns %>% 
                                          dplyr::select(-week))
# generate time column 
mumps_res_nar_1_1 <- mumps_res_nar_1_1 %>% 
  mutate(week = seq(2, nrow(mumps_res_nar_1_1) + 1))

# visualise residuals for selected counties according to original paper
ggplot(data = mumps_res_nar_1_1, 
       aes(x = week, 
           y = Bedford)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Bedfordshire (Bedford)")
ggsave("plots/modelreplicate/mumps_nar_1_1_bedford_residuals.pdf", 
       width = 14, height = 14, unit = "cm")

ggplot(data = mumps_res_nar_1_1, 
       aes(x = week, 
           y = Aylesbury)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Buckinghamshire (Aylesbury)")
ggsave("plots/modelreplicate/mumps_nar_1_1_aylesbury_residuals.pdf", 
       width = 14, height = 14, unit = "cm")

ggplot(data = mumps_res_nar_1_1, 
       aes(x = week, 
           y = Cambridge)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Cambridgeshire (Cambridge)")
ggsave("plots/modelreplicate/mumps_nar_1_1_cambridge_residuals.pdf", 
       width = 14, height = 14, unit = "cm")

ggplot(data = mumps_res_nar_1_1, 
       aes(x = week, 
           y = Chester)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Cheshire (Chester)")
ggsave("plots/modelreplicate/mumps_nar_1_1_chester_residuals.pdf", 
       width = 14, height = 14, unit = "cm")

# visualise residuals as QQ plot
g1 <- ggplot(data = mumps_res_nar_1_1, 
       aes(sample = Bedford)) +
  geom_qq() +
  geom_qq_line() +
  xlab("theor. quantiles") +
  ylab("emp. quantiles for Bedfordshire")
ggsave("plots/modelreplicate/qq_nar_1_1_bedford.pdf", 
       width = 14, height = 14, unit = "cm")

g2 <- ggplot(data = mumps_res_nar_1_1, 
       aes(sample = Aylesbury)) +
  geom_qq() +
  geom_qq_line() +
  xlab("theor. quantiles") +
  ylab("emp. quantiles for Buckinghamshire")
ggsave("plots/modelreplicate/qq_nar_1_1_buckinghamshire.pdf", 
       width = 14, height = 14, unit = "cm")

g3 <- ggplot(data = mumps_res_nar_1_1, 
       aes(sample = Cambridge)) +
  geom_qq() +
  geom_qq_line()+ 
  xlab("theor. quantiles") +
  ylab("emp. quantiles for Cambridgeshire")
ggsave("plots/modelreplicate/qq_nar_1_1_cambridge.pdf", 
       width = 14, height = 14, unit = "cm")

g4 <- ggplot(data = mumps_res_nar_1_1, 
       aes(sample = Chester)) +
  geom_qq() +
  geom_qq_line() +
  xlab("theor. quantiles") +
  ylab("emp. quantiles for Cheshire")
ggsave("plots/modelreplicate/qq_nar_1_1_cheshire.pdf", 
       width = 14, height = 14, unit = "cm")



## ---- NAR(1, 1) log-transformed ----------------------------------------------
nar_log_1_1 <- GNARfit(vts = mumps_log_ts, 
                       net = mumps_gnar, 
                       alphaOrder = 1, 
                       betaOrder = 1, 
                       globalalpha = TRUE)


summary(nar_log_1_1) 


# PREDICTIONS
mumps_pred_nar_log_1_1 <- fitted(nar_log_1_1) %>% 
  as.data.frame()
# assign county towns as column names
colnames(mumps_pred_nar_log_1_1) <- colnames(mumps_central_towns %>% 
                                               dplyr::select(-week))
# generate time column 
mumps_pred_nar_log_1_1 <- mumps_pred_nar_log_1_1 %>% 
  mutate(week = seq(2, nrow(mumps_pred_nar_log_1_1) + 1))

# visualise predicted values for county Oxford 
ggplot(data = mumps_central_towns, 
       aes(x = week, 
           y = log(Oxford + 1))) +
  geom_line() +
  geom_line(data = mumps_pred_nar_log_1_1, 
            aes(x = week, 
                y = Oxford), 
            color = "red")
ggsave("plots/modelreplicate/mumps_nar_1_1_log_oxford.pdf", 
       width = 23, height = 14, unit = "cm")

# RESIDUALS
mumps_res_nar_log_1_1 <- residuals(nar_log_1_1) %>% 
  as.data.frame()

colnames(mumps_res_nar_log_1_1) <- colnames(mumps_central_towns %>% 
                                              dplyr::select(-week))

mumps_res_nar_log_1_1 <- mumps_res_nar_log_1_1 %>% 
  mutate(week = seq(2, nrow(mumps_res_nar_log_1_1) + 1))

# visualise residuals 
ggplot(data = mumps_res_nar_log_1_1, 
       aes(x = week, 
           y = Bedford)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Bedfordshire (Bedford)")
ggsave("plots/modelreplicate/mumps_nar_1_1_log_bedford_residuals.pdf", 
       width = 14, height = 14, unit = "cm")

ggplot(data = mumps_res_nar_log_1_1, 
       aes(x = week, 
           y = Aylesbury)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Buckinghamshire (Aylesbury)")
ggsave("plots/modelreplicate/mumps_nar_1_1_log_aylesbury_residuals.pdf", 
       width = 14, height = 14, unit = "cm")

ggplot(data = mumps_res_nar_log_1_1, 
       aes(x = week, 
           y = Cambridge)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Cambridgeshire (Cambridge)")
ggsave("plots/modelreplicate/mumps_nar_1_1_log_cambridge_residuals.pdf", 
       width = 14, height = 14, unit = "cm")

ggplot(data = mumps_res_nar_log_1_1, 
       aes(x = week, 
           y = Chester)) +
  geom_point() +
  xlab("Time") +
  ylab("Residuals for Cheshire (Chester)")
ggsave("plots/modelreplicate/mumps_nar_1_1_log_chester_residuals.pdf", 
       width = 14, height = 14, unit = "cm")



## ---- NAR(1, 0) --------------------------------------------------------------
nar_1_0 <- GNARfit(vts = mumps_ts, 
                   net = mumps_gnar, 
                   alphaOrder = 1, 
                   betaOrder = 0, 
                   globalalpha = TRUE)
summary(nar_1_0)

## ---- NAR(1, 0) log-transformed ----------------------------------------------
nar_log_1_0 <- GNARfit(vts = mumps_log_ts, 
                       net = mumps_gnar, 
                       alphaOrder = 1, 
                       betaOrder = 0, 
                       globalalpha = TRUE)
summary(nar_log_1_0)

## ---- NAR(1,  2) -------------------------------------------------------------
nar_1_2 <- GNARfit(vts = mumps_ts, 
                   net = mumps_gnar, 
                   alphaOrder = 1, 
                   betaOrder = 2, 
                   globalalpha = TRUE)
summary(nar_1_2)

## ---- NAR(1, 2) log-transformed ----------------------------------------------
nar_log_1_2 <- GNARfit(vts = mumps_log_ts, 
                       net = mumps_gnar, 
                       alphaOrder = 1, 
                       betaOrder = 2, 
                       globalalpha = TRUE)
summary(nar_log_1_2)

## ---- NAR(2, [1, 0]) ---------------------------------------------------------
nar_2_10 <- GNARfit(vts = mumps_ts, 
                    net = mumps_gnar, 
                    alphaOrder = 2, 
                    betaOrder = c(1, 0), 
                    globalalpha = TRUE)
summary(nar_2_10)

# RESIDUALS
mumps_res_nar_2_10 <- residuals(nar_2_10) %>% 
  as.data.frame()
# assign county towns as column names
colnames(mumps_res_nar_2_10) <- colnames(mumps_central_towns %>% 
                                          dplyr::select(-week))

## ---- NAR(2, [1, 0]) log-transformed -----------------------------------------
nar_log_2_10 <- GNARfit(vts = mumps_log_ts, 
                       net = mumps_gnar, 
                       alphaOrder = 2, 
                       betaOrder = c(1, 0), 
                       globalalpha = TRUE)
summary(nar_log_2_10)

## ---- NAR(2, [1, 1]) ---------------------------------------------------------
nar_2_11 <- GNARfit(vts = mumps_ts, 
                    net = mumps_gnar, 
                    alphaOrder = 2, 
                    betaOrder = c(1, 1), 
                    globalalpha = TRUE)
summary(nar_2_11)

## ---- NAR(2, [1, 1]) log-transformed -----------------------------------------
nar_log_2_11 <- GNARfit(vts = mumps_log_ts, 
                        net = mumps_gnar, 
                        alphaOrder = 2, 
                        betaOrder = c(1, 1), 
                        globalalpha = TRUE)
summary(nar_log_2_11)

## ---- Model comparison -------------------------------------------------------
# comparison of different models 
model_comparison <- data_frame("NAR(1, 0)" = c(BIC(nar_1_0), 
                                               AIC(nar_1_0), 
                                               RSS(nar_1_0)), 
                               "NAR(1,1)" = c(BIC(nar_1_1), 
                                              AIC(nar_1_1), 
                                              RSS(nar_1_1)), 
                               "log NAR(1, 1)" = c(BIC(nar_log_1_1), 
                                                   AIC(nar_log_1_1), 
                                                   RSS(nar_log_1_1)), 
                               "NAR(1, 2)" = c(BIC(nar_1_2), 
                                               AIC(nar_1_2), 
                                               RSS(nar_1_2)), 
                               "NAR(2, [1, 0])" = c(BIC(nar_2_10), 
                                                    AIC(nar_2_10), 
                                                    RSS(nar_2_10)), 
                               "NAR(2, [1, 1])" = c(BIC(nar_2_11),
                                                    AIC(nar_2_11),  
                                                    RSS(nar_2_11))
                               ) %>% t() %>% as.data.frame()
colnames(model_comparison) <- c("BIC", "AIC", "RSS")

# table for latex
strCaption <- "Comparison of NARIMA models through BIC, AIC and residual 
sum of squares (RSS)"
print(xtable(model_comparison,
             digits=2,
             caption=strCaption,
             label="tab:narima_model_comparison", 
             align = c("l", "|", "r", "r", "r")),
      include.rownames=TRUE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(model_comparison)),
                        command = c(paste("\\toprule \n",
                                          "Model & BIC & AIC & RSS \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

## ---- Model comparison log-transformed ---------------------------------------
# comparison of different models
model_comparison_log <- data_frame("log NAR(1, 0)" = c(BIC(nar_log_1_0), 
                                                       AIC(nar_log_1_0), 
                                                       RSS(nar_log_1_0), 
                                                       1212.2), 
                                   "log NAR(1,1)" = c(BIC(nar_log_1_1), 
                                                      AIC(nar_log_1_1), 
                                                      RSS(nar_log_1_1), 
                                                      1029.8),  
                                   "log NAR(1, 2)" = c(BIC(nar_log_1_2), 
                                                       AIC(nar_log_1_2), 
                                                       RSS(nar_log_1_2), 
                                                       NA), 
                                   "log NAR(2, [1, 0])" = c(BIC(nar_log_2_10), 
                                                            AIC(nar_log_2_10), 
                                                            RSS(nar_log_2_10), 
                                                            862.2), 
                                   "log NAR(2, [1, 1])" = c(BIC(nar_log_2_11), 
                                                            AIC(nar_log_2_11),  
                                                            RSS(nar_log_2_11), 
                                                            862.1)
                                   ) %>% 
  t() %>% 
  as.data.frame()
colnames(model_comparison_log) <- c("BIC", "AIC", "RSS", "orig. RSS")

# table for latex
strCaption <- "Comparison of NARIMA models through BIC, AIC and residual 
sum of squares (RSS) on log transformed data, compared to original RSS in 
\\cite{knight2016modelling}"

print(xtable(model_comparison_log,
             digits=2,
             caption=strCaption,
             label="tab:log_narima_model_comparison", 
             align = c("l", "|", "r", "r", "r", "r")),
      include.rownames=TRUE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(model_comparison_log)),
                        command = c(paste("\\toprule \n",
                                          "Model & BIC & AIC & RSS & 
                                          orig. RSS \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


## ---- ACF / PACF -------------------------------------------------------------
# cross-covariance analysis for residuals (Somerset Taunton, Avon Bristol)
pdf("plots/modelreplicate/mumps_nar_1_1_acf_bristol_taunton.pdf")
acf(mumps_res_nar_1_1[, c("Bristol", "Taunton")])
dev.off()

pdf("plots/modelreplicate/mumps_nar_1_1_pacf_bristol_taunton.pdf")
pacf(mumps_res_nar_1_1[, c("Bristol", "Taunton")])
dev.off()

pdf("plots/modelreplicate/mumps_nar_log_1_1_acf_bristol_taunton.pdf")
acf(mumps_res_nar_log_1_1[, c("Bristol", "Taunton")])
dev.off()

pdf("plots/modelreplicate/mumps_nar_log_1_1_pacf_bristol_taunton.pdf")
pacf(mumps_res_nar_log_1_1[, c("Bristol", "Taunton")])
dev.off()

pdf("plots/modelreplicate/mumps_nar_2_10_acf_bristol_taunton.pdf")
acf(mumps_res_nar_2_10[, c("Bristol", "Taunton")])
dev.off()

# explore other counties for NAR(1, 1)
acf(mumps_res_nar_1_1[, c("Newport", "Newcastle")])
acf(mumps_res_nar_1_1[, c("Rhayader", "Norwich")])
acf(mumps_res_nar_1_1[, c("London", "Hertford")])
acf(mumps_res_nar_1_1[, c("Liverpool", "Lancaster")])


## ---- AR models --------------------------------------------------------------
# obtain county specific data sets for Somerset and Avon
mumps_somerset <- mumps_log_ts[, colnames(mumps_log_ts) == "Taunton"] 
mumps_avon <- mumps_log_ts[, colnames(mumps_log_ts) == "Bristol"] 

# fit AR model and compute residuals, order according to original paper
mod_somerset <- arima(mumps_somerset, order = c(1, 0, 0))
res_somerset <- residuals(mod_somerset)

mod_avon <- arima(mumps_avon, order = c(2, 0, 0))
res_avon <- residuals(mod_avon)

# RESIDUALS
res_somerset_avon <- data.frame("Avon" = res_avon, 
                                "Somerset" = res_somerset)
rownames(res_somerset_avon) <- names(mumps_avon)

# plot ACF for residuals 
pdf("plots/modelreplicate/mumps_acf_ar_taunton_bristol.pdf")
acf(res_somerset_avon)
dev.off()

# data frame including mumps cases 
mumps_somerset_avon <- data.frame("Avon" = mumps_avon, 
                                  "Somerset" = mumps_somerset)

# plot ACF for mumps cases
pdf("plots/modelreplicate/mumps_acf_taunton_bristol.pdf")
acf(mumps_somerset_avon)
dev.off()


## ---- Reproducibility assessment ---------------------------------------------
original_coef <- list(
  "NAR(1, [1])" = c(0.683, 0.263), 
  "log NAR(1, [1])" = c(0.647, 0.33), 
  "NAR(2, [1, 0])" = c(0.394, 0.204, 0.380)
)

fitted_coef <- list(
  "NAR(1, [1])" = coef(nar_1_1), 
  "log NAR(1, [1])" = coef(nar_log_1_1), 
  "NAR(2, [1, 0])" = coef(nar_2_10)
)

fitted_coef_ci <- list(
  "NAR(1, [1])" = confint(nar_1_1$mod) %>% round(digits = 3), 
  "log NAR(1, [1])" = confint(nar_log_1_1$mod) %>% round(digits = 3), 
  "NAR(2, [1, 0])" = confint(nar_2_10$mod) %>% round(digits = 3)
)

## ---- Deep dive residuals ----------------------------------------------------
# generate vector incorporating all counties 
counties <- mumps_ts %>% colnames()

# plot residual scatter plot for all counties for NAR(1, 1) model 
for (county in counties) {
  res_county <- mumps_res_nar_1_1 %>% 
    dplyr::select(all_of(county), week)
  
  colnames(res_county) <- c("county", "week")
  
  ggplot(data = res_county, 
         aes(x = week, 
             y = county)) +
    geom_point() +
    xlab("Time") +
    ylab(paste0("Residuals for ", county))
  ggsave(paste0("plots/modelreplicate/NAR_1_1_residuals/", county, ".pdf"), 
         width = 14, height = 14, unit = "cm")
}

