### CHANGE IN GNAR MODELS FOR DATA SUBSETS ACCORDING TO COVID-19 
### REGULATIONS 

# set seed to guarantee reproducibility 
set.seed(1234)

# Load libraries and objects ----------------------------------------------
library(readr)
library(igraph)
library(GNAR)
library(MASS) # for box cox 
library(tidyverse)
library(magrittr) # for pipping 
library(xtable) # for tables 
library(spdep) # for neighbourhood construction 
library(expp)


# load vectors and GNAR objects 
load(file = "data/RObjects/GNAR.RData")
load(file = "data/RObjects/population_weight.RData")
load(file = "data/RObjects/distance_urbanisation.RData")
load(file = "data/RObjects/county_index.RData")
load(file = "data/RObjects/coord_urbanisation.RData")

source("functions.R")

# turn off warnings
options(warn = -1)

# set ggplot theme 
theme_set(theme_bw(base_size = 16))

# Load data ---------------------------------------------------------------
COVID_weekly_data <- read_csv(file = "data/COVID/ireland_covid_weekly.csv", 
                              show_col_types = FALSE)

# convert to correct format for GNARfit()
covid_cases <-  COVID_weekly_data %>% 
  dplyr::select(CountyName,
                yw, 
                weeklyCases) %>% 
  spread(CountyName, 
         weeklyCases) %>% 
  column_to_rownames(var = "yw") %>% 
  as.matrix()

# convert to data frame with time column for plotting 
covid_cases_df <- covid_cases %>% 
  as.data.frame() %>% 
  mutate(time = as.Date(rownames(covid_cases))) 

# Relevant restrictions ---------------------------------------------------
# create data frame for relevant COVID-19 regulations and their start date
restrictions_df <- data.frame(date = c("28.02.2020", 
                                       "18.08.2020" ,
                                       "26.12.2020" ,
                                       "26.07.2021" ,
                                       "06.03.2022", 
                                       "18.06.2022"
                                       ) %>% as.Date(format = "%d.%m.%Y"),
                              content = c( "Start", 
                                          "County-specific lockdowns", 
                                          "Level-5 lockdown", 
                                          "Indoor dining", 
                                          "Ease", 
                                          "End"
                              ), 
                              index = seq(1, 6), 
                              sentiment = c("Restriction", 
                                            "Restriction", 
                                            "Restriction",
                                            "Ease",  
                                            "Ease",  
                                            "Ease"
                                            )
                              )

# create data frame with COVID-19 virus strains and the date at which more 
# than 50% of COVID-19 infections where caused by the strain 
strains <- data.frame(variant = c("Original", 
                                  "Alpha", 
                                  "Delta", 
                                  "Omicron I", 
                                  "Omicron II", 
                                  NA), 
                      onset = c("28/02/2020", 
                                "27/12/2020", 
                                "06/06/2021", 
                                "13/12/2021", 
                                "13/03/2022", 
                                "18/06/2022") %>% 
                        as.Date(format = "%d/%m/%Y"), 
                      index = seq(1, 6))

COVID_weekly_data$restriction <- NA
COVID_weekly_data$restriction_date <- as.Date(NA)
COVID_weekly_data$sentiment <- NA
COVID_weekly_data$strain <- NA

# assign each week the predominant COVID-19 virus strain  
for (i in seq(1, nrow(COVID_weekly_data))) {
  current_date <- COVID_weekly_data$yw[i]
  
  which_strain <- strains[strains$onset > current_date, ]$index %>% min() - 1
  
  COVID_weekly_data$strain[i] <- strains[which_strain, ]$variant
  
  which_restriction <- restrictions_df[restrictions_df$date > current_date, ]$index %>% 
    min() - 1
  
  if (which_restriction == 1) {
    which_restriction <- 2
  }
  COVID_weekly_data$restriction[i] <- which_restriction
  COVID_weekly_data$restriction_date[i] <- restrictions_df[which_restriction, ]$date
  COVID_weekly_data$sentiment[i] <- restrictions_df[which_restriction, ]$sentiment
}


# plot restrictions and virus strains 
ggplot(COVID_weekly_data) +
  geom_rect(aes(xmin = yw, xmax = lead(yw),
                ymin = -Inf, ymax = Inf,
                fill = strain)) +
  geom_line(aes(x = yw, 
                y = weeklyCases, 
                group = CountyName)) +
  geom_vline(aes(xintercept = restriction_date, 
                 color = restriction %>% as.factor(), 
                 linetype = sentiment), 
             size = 1) +
  xlab("Time") +
  ylab("COVID-19 ID") + 
  scale_color_brewer(palette = "Set1") +
  scale_fill_manual(values = c("Original" = "#FAFAFA", 
                               "Alpha" = "#eeeeee", 
                               "Delta" = "#E1E5E8", 
                               "Omicron I" = "#D0D5D9", 
                               "Omicron II" = "#ABB0B8")) +
  theme(legend.position = "none")
ggsave("plots/dataVisualisation/covid_id_restrictions.pdf", 
       width = 23, height = 14, unit = "cm")



# Split dataset -----------------------------------------------------------
splits <- c("27.02.2020",  
            "18.08.2020",
            "26.12.2020",
            "26.07.2021", 
            "06.03.2022",
            "18.06.2022") %>% 
  as.Date(format = "%d.%m.%Y")

# split COVID-19 data set into separate data sets according to COVID-19 
# regulations 
datasets_list <- list()

for (i in seq(2, length(splits))) {
  datasets_list[[i-1]] <- COVID_weekly_data %>% 
    filter(splits[i-1] < yw & yw <= splits[i]) %>% 
    dplyr::select(CountyName,
                  yw, 
                  weeklyCases) %>% 
    spread(CountyName, 
           weeklyCases) %>% 
    column_to_rownames(var = "yw") %>% 
    as.matrix()
}

# Variability for each data set --------------------------------------------
# compute variance in 1-lag COVID-19 ID for each data subset
lapply(datasets_list, FUN = function(i) {
  apply(i, MARGIN = 2, FUN = function(j) {sd(j)}) %>% mean()
})

# Change in coefficients for best performing models -----------------------
# compute how the coefficients for the best performing model for each network
# change if fitted to the different data subsets

# Queen 
param_queen <- parameter_development(county_index = county_index_queen, 
                                     net = covid_net_queen_gnar, 
                                     alpha = 4, 
                                     beta = c(2, 2, 1, 1), 
                                     globalalpha = TRUE,
                                     old = TRUE, 
                                     network_name = "queen")

# Economic hubs 
param_eco_hubs <- parameter_development(county_index = county_index_eco_hubs, 
                                        net = covid_net_eco_hubs_gnar, 
                                        alpha = 5, 
                                        beta = c(1, 1, 1, 0, 0), 
                                        globalalpha = TRUE, 
                                        old = TRUE, 
                                        network_name = "eco_hubs", 
                                        numeric_vertices = TRUE)

# Railway-based
param_train <- parameter_development(county_index = county_index_train, 
                                     net = covid_net_train_gnar, 
                                     alpha = 4, 
                                     beta = c(2, 1, 1, 1), 
                                     globalalpha = TRUE, 
                                     old = TRUE, 
                                     network_name = "train", 
                                     numeric_vertices = TRUE)

# Delaunay 
param_delaunay <- parameter_development(county_index = county_index_delaunay, 
                                        net = covid_net_delaunay_gnar, 
                                        alpha = 4, 
                                        beta = c(2, 2, 1, 1), 
                                        globalalpha = TRUE, 
                                        old = TRUE, 
                                        network_name = "delaunay")


# Gabriel 
param_gabriel <- parameter_development(county_index = county_index_gabriel, 
                                       net = covid_net_gabriel_gnar, 
                                       alpha = 5, 
                                       beta = c(1, 1, 1, 0, 0), 
                                       globalalpha = TRUE, 
                                       old = TRUE, 
                                       network_name = "gabriel")

# Relative
param_relative <- parameter_development(county_index = county_index_relative, 
                                        net = covid_net_relative_gnar, 
                                        alpha = 4, 
                                        beta = c(2, 2, 1, 1), 
                                        globalalpha = TRUE,
                                        old = TRUE, 
                                        network_name = "relative")

# SOI
param_soi <- parameter_development(county_index = county_index_soi, 
                                   net = covid_net_soi_gnar, 
                                   alpha = 4, 
                                   beta = c(2, 1, 1, 1), 
                                   globalalpha = TRUE,
                                   old = TRUE, 
                                   network_name = "soi")

# KNN
param_knn <- parameter_development(county_index = county_index_opt_knn, 
                                   net = opt_knn_net_gnar, 
                                   alpha = 4, 
                                   beta = c(1, 1, 0, 0), 
                                   globalalpha = TRUE, 
                                   old = TRUE, 
                                   network_name = "KNN")


# DNN 
param_dnn <- parameter_development(county_index = county_index_opt_dnn, 
                                   net = opt_dnn_net_gnar, 
                                   alpha = 4, 
                                   beta = c(2, 2, 2, 1), 
                                   globalalpha = TRUE, 
                                   old = FALSE, 
                                   inverse_distance = TRUE, 
                                   network_name = "DNN")



# Complete 
param_complete <- parameter_development(county_index = county_index_complete, 
                                        net = complete_net_gnar, 
                                        alpha = 5, 
                                        beta = c(1, 1, 1, 1, 0), 
                                        globalalpha = TRUE, 
                                        network_name = "complete")


# Best performing model for every subset ----------------------------------
# compute GNAR models for each data subsets and select the best performing one 
# based on the BIC 

# Queen 
best_for_subset_queen <- fit_and_predict_for_restrictions(net = covid_net_queen_gnar)
best_for_subset_queen$network <- "Queen"

# Economic hub
best_for_subset_eco_hub <- fit_and_predict_for_restrictions(net = covid_net_eco_hubs_gnar, 
                                                            numeric_vertices = TRUE, 
                                                            county_index = county_index_eco_hubs, 
                                                            upper_limit = 4)
best_for_subset_eco_hub$network <- "Eco. hub"

# Railway-based 
best_for_subset_train <- fit_and_predict_for_restrictions(net = covid_net_train_gnar, 
                                                          numeric_vertices = TRUE, 
                                                          county_index = county_index_train)
best_for_subset_train$network <- "Train"

# Delaunay triangulation 
best_for_subset_delaunay <- fit_and_predict_for_restrictions(net = covid_net_delaunay_gnar)
best_for_subset_delaunay$network <- "Delaunay"

# Gabriel 
best_for_subset_gabriel <- fit_and_predict_for_restrictions(net = covid_net_gabriel_gnar)
best_for_subset_gabriel$network <- "Gabriel"

# Relative neighbourhood
best_for_subset_relative <- fit_and_predict_for_restrictions(net = covid_net_relative_gnar)
best_for_subset_relative$network <- "Relative"

# SOI 
best_for_subset_soi <- fit_and_predict_for_restrictions(net = covid_net_soi_gnar)
best_for_subset_soi$network <- "SOI"

# Complete
best_for_subset_complete <- fit_and_predict_for_restrictions(net = complete_net_gnar, 
                                                             upper_limit = 1)
best_for_subset_complete$network <- "Complete"

# compare best-performing models across networks
best_for_subset <- rbind(best_for_subset_queen,
                         best_for_subset_eco_hub, 
                         best_for_subset_train,
                         best_for_subset_delaunay, 
                         best_for_subset_gabriel, 
                         best_for_subset_relative, 
                         best_for_subset_soi, 
                         best_for_subset_complete)

# for latex 
strCaption <- "Overview over the best performing model for every COVID-19 data 
subset for every COVID-19 network, excluding the KNN and DNN network"
print(xtable(best_for_subset[, c(4, 1, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_subsets", 
             align = c("", "l", "|", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_for_subset[, c(4, 1, 2, 3)])),
                        command = c(paste("\\toprule \n",
                                          "Network & data subset & best model &
                                          BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# construct KNN networks with different neighbourhood size k, fit GNAR models 
# and select best performing model via the BIC
knn_best <- list()

for (k in seq(1, 26, by = 2)) {
  # create nb list
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
  
  # create ordered county index data frame 
  county_index_knn <- data.frame("CountyName" = covid_net_knn_igraph %>%
                                   V() %>% 
                                   names(), 
                                 "index" = seq(1, 26))
  
  # compute an upper limit for neighbourhood stage 
  upper_limit_knn <- distance_table(covid_net_knn_igraph)$res %>% length()
  
  # fit GNAR models and select the best performing one for each data subset 
  res <- fit_and_predict_for_restrictions(net = covid_net_knn, 
                                          upper_limit = upper_limit_knn)
  
  res$hyperparam <- k
  
  # save best performing model for every k across all data subsets  
  knn_best[[length(knn_best) + 1]] <- res
  
}

# filter the best performing GNAR model for each data subset across all 
# neighbourhood sizes 
knn_best_df <- do.call(rbind.data.frame, knn_best) %>% 
  group_by(data_subset) %>% 
  filter(BIC == min(BIC)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  arrange(data_subset)
knn_best_df$network <-  "KNN"



# construct DNN networks for different distance thresholds d, fit GNAR models 
# and select the best performing one via the BIC 
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
  
  # compute an upper limit for neighbourhood stage 
  upper_limit_dnn <- distance_table(covid_net_dnn_igraph)$res %>% length()
  
  # fit GNAR models and select the best performing one for each data subset
  res <- fit_and_predict_for_restrictions(net = covid_net_dnn, 
                                          upper_limit = upper_limit_dnn)
  
  res$hyperparam <- d
  
  # save best performing model for every distance threshold d
  dnn_best[[length(dnn_best) + 1]] <- res
  
}

# filter best performing model across distance threshold for each data subset
dnn_best_df <- do.call(rbind.data.frame, dnn_best) %>% 
  group_by(data_subset) %>% 
  filter(BIC == min(BIC)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  arrange(data_subset)
dnn_best_df$network <- "DNN"


best_subset_knn_dnn <- rbind(knn_best_df,
                             dnn_best_df)

# for latex 
strCaption <- "Overview over the best performing model and optimal 
neighbourhood size $k$ / distance threshold $d$ for the KNN and DNN network"
print(xtable(best_subset_knn_dnn[, c(5, 1, 4, 2, 3)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_knn_dnn_subsets", 
             align = c("", "l", "|", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_subset_knn_dnn[, c(5, 4, 1, 2, 3)])),
                        command = c(paste("\\toprule \n",
                                          "Network & data subset & k / d [in km] &
                                          \\code{GNAR} model & BIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)



# Construct networks ------------------------------------------------------
# construct the networks for the best performing GNAR models 

# DNN d = 300
dnn_300 <- dnearneigh(x = coord_urbanisation, 
                     d1 = 0, 
                     d2 = 300,
                     row.names = coord_urbanisation %>% rownames(),
                     longlat = TRUE, 
                     use_s2 = TRUE)

dnn_300_igraph<- neighborsDataFrame(nb = dnn_300) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_300_gnar <- dnn_300_igraph %>% 
  igraphtoGNAR()

# DNN d = 175
dnn_175 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 175,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_175_igraph<- neighborsDataFrame(nb = dnn_175) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_175_gnar <- dnn_175_igraph %>% 
  igraphtoGNAR()

# DNN d = 225 
dnn_225 <- dnearneigh(x = coord_urbanisation, 
                      d1 = 0, 
                      d2 = 225,
                      row.names = coord_urbanisation %>% rownames(),
                      longlat = TRUE, 
                      use_s2 = TRUE)

dnn_225_igraph<- neighborsDataFrame(nb = dnn_225) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

dnn_225_gnar <- dnn_225_igraph %>% 
  igraphtoGNAR()

# KNN k = 7
knn_7 <- knearneigh(x = coord_urbanisation, 
                    k = 7, 
                    longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_7_igraph<- neighborsDataFrame(nb = knn_7) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_7_gnar <- knn_7_igraph %>% 
  igraphtoGNAR()

# KNN k = 9
knn_9 <- knearneigh(x = coord_urbanisation, 
                    k = 9, 
                    longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_9_igraph<- neighborsDataFrame(nb = knn_9) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_9_gnar <- knn_9_igraph %>% 
  igraphtoGNAR()

# KNN k = 11
knn_11 <- knearneigh(x = coord_urbanisation, 
                    k = 11, 
                    longlat = TRUE) %>% 
  knn2nb(row.names = coord_urbanisation %>% rownames(),)

knn_11_igraph<- neighborsDataFrame(nb = knn_11) %>% 
  graph_from_data_frame(directed = FALSE) %>% 
  igraph::simplify() 

knn_11_gnar <- knn_11_igraph %>% 
  igraphtoGNAR()


# Compare networks --------------------------------------------------------
# compute network characteristics

# DNN 175
char_dnn_175 <- network_characteristics(dnn_175_igraph, 
                                        network_name = "dnn_175")

# DNN 225
char_dnn_225 <- network_characteristics(dnn_225_igraph, 
                                        network_name = "dnn_225")

# DNN 300
char_dnn_300 <- network_characteristics(dnn_300_igraph, 
                                        network_name = "dnn_300")

# KNN 9
char_knn_9 <- network_characteristics(knn_9_igraph, 
                                      network_name = "knn_9")

# Delaunay 
char_delaunay <- network_characteristics(covid_net_delaunay_igraph, 
                                         network_name = "delaunay")

# compare network characteristics 
cbind(char_dnn_175, 
      char_dnn_225[, 2], 
      char_dnn_300[, 2], 
      char_knn_9[, 2], 
      char_delaunay[, 2])



# Best GNAR model for each data subset ------------------------------------
# create a data frame with the best performing GNAR model and network for 
# each data subset 
best_subset_knn_dnn_final <- best_subset_knn_dnn %>% 
  mutate(network = paste(network, hyperparam, sep = "-")) %>% 
  dplyr::select(-hyperparam)

best_for_subset_large_df <- rbind(best_for_subset, 
                                  best_subset_knn_dnn_final)

best_for_subset_all <-  best_for_subset_large_df %>% 
  group_by(data_subset) %>% 
  filter(BIC == min(BIC)) %>% 
  arrange(data_subset)

# fit models 
# data set 1
best_model_1 <- fit_and_predict(alpha = 5, 
                                beta = c(1, 0, 0, 0, 0), 
                                net = dnn_300_gnar, 
                                vts = datasets_list[[1]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 2
best_model_2 <- fit_and_predict(alpha = 5, 
                                beta = c(1, 1, 1, 1, 1), 
                                net = dnn_175_gnar, 
                                vts = datasets_list[[2]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 3
best_model_3 <- fit_and_predict(alpha = 5, 
                                beta = c(3, 1, 1, 0, 0), 
                                net = knn_9_gnar, 
                                vts = datasets_list[[3]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 4
best_model_4 <- fit_and_predict(alpha = 5, 
                                beta = c(1, 1, 1, 0, 0), 
                                net = dnn_225_gnar, 
                                vts = datasets_list[[4]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)
# data set 5
best_model_5 <- fit_and_predict(alpha = 5, 
                                beta = c(4, 2, 2, 1, 1), 
                                net = covid_net_eco_hubs_gnar, 
                                vts = datasets_list[[5]], 
                                globalalpha = TRUE, 
                                old = TRUE,
                                forecast_window = 5, 
                                return_model = TRUE)

best_for_subset_all$AIC <- c(AIC(best_model_1), 
                             AIC(best_model_2), 
                             AIC(best_model_3), 
                             AIC(best_model_4), 
                             AIC(best_model_5))


# for latex 
strCaption <- "Overview over the best performing model and network for every 
COVID-19 data subset"
print(xtable(best_for_subset_all[, c(1, 4, 2, 3, 5)],
             digits=2,
             caption=strCaption,
             label="tab:best_model_datasets", 
             align = c("", "l", "|", "r", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(best_for_subset_all[, c(1, 4, 2, 3, 5)])),
                        command = c(paste("\\toprule \n",
                                          "Data subset & network & \\code{GNAR} model & 
                                          BIC & AIC \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# Data formatting ----------------------------------------------------------
# format data subsets as data frame with time column 
data_1 <- datasets_list[[1]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[1]]) %>% as.Date())

data_2 <- datasets_list[[2]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[2]])  %>% as.Date())

data_3 <- datasets_list[[3]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[3]])  %>% as.Date())

data_4 <- datasets_list[[4]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[4]])  %>% as.Date())

data_5 <- datasets_list[[5]] %>% 
  as.data.frame() %>% 
  mutate(time = rownames(datasets_list[[5]])  %>% as.Date())



# MASE for data subset 1 --------------------------------------------------
# fit best performing GNAR model for each network for data subset 1

# Queen
mod_1_queen <- fit_and_predict(alpha = 5, 
                               beta = c(3, 2, 2, 1, 1), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list[[1]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)

# Eco hub 
mod_1_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(4, 1, 1, 0, 0), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list[[1]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_1_train <- fit_and_predict(alpha = 5, 
                               beta = c(1, 0, 0, 0, 0), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list[[1]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_1_delaunay <- fit_and_predict(alpha = 5, 
                                  beta = c(2, 2, 0, 0, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_1_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(3, 1, 1, 0, 0), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list[[1]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_1_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(4, 0, 0, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_1_soi <- fit_and_predict(alpha = 5, 
                             beta = c(4, 1, 1, 1, 0), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_1_knn <- fit_and_predict(alpha = 5, 
                             beta = c(4, 1, 1, 0, 0), 
                             net = knn_7_gnar, 
                             vts = datasets_list[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_1_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 0, 0, 0, 0), 
                             net = dnn_300_gnar, 
                             vts = datasets_list[[1]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_1_complete <- fit_and_predict(alpha = 4, 
                                  beta = c(1, 1, 0, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list[[1]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)


# compute MASE for best performing GNAR model for each network
mase_1_queen <- compute_MASE(model = mod_1_queen, 
                             network_name = "subset_1_queen", 
                             n_ahead = 5, 
                             data_df = data_1)

mase_1_eco_hub <- compute_MASE(model = mod_1_eco_hub, 
                               network_name = "subset_1_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_1)

mase_1_train <- compute_MASE(model = mod_1_train, 
                             network_name = "subset_1_train", 
                             n_ahead = 5, 
                             data_df = data_1)

mase_1_delaunay <- compute_MASE(model = mod_1_delaunay, 
                                network_name = "subset_1_delaunay", 
                                n_ahead = 5, 
                                data_df = data_1)

mase_1_gabriel <- compute_MASE(model = mod_1_gabriel, 
                               network_name = "subset_1_gabriel", 
                               n_ahead = 5, 
                               data_df = data_1)

mase_1_relative <- compute_MASE(model = mod_1_relative, 
                                network_name = "subset_1_relative", 
                                n_ahead = 5, 
                                data_df = data_1)

mase_1_soi <- compute_MASE(model = mod_1_soi, 
                           network_name = "subset_1_soi", 
                           n_ahead = 5, 
                           data_df = data_1)

mase_1_knn <- compute_MASE(model = mod_1_knn, 
                           network_name = "subset_1_knn", 
                           n_ahead = 5, 
                           data_df = data_1)

mase_1_dnn <- compute_MASE(model = mod_1_dnn, 
                           network_name = "subset_1_dnn", 
                           n_ahead = 5, 
                           data_df = data_1)

mase_1_complete <- compute_MASE(model = mod_1_complete, 
                                network_name = "subset_1_complete", 
                                n_ahead = 5, 
                                data_df = data_1)


mase_1_overview <- rbind(mase_1_queen, 
                         mase_1_eco_hub, 
                         mase_1_train, 
                         mase_1_delaunay, 
                         mase_1_gabriel, 
                         mase_1_relative, 
                         mase_1_soi, 
                         mase_1_knn, 
                         mase_1_dnn, 
                         mase_1_complete)


# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network
ggplot(mase_1_overview %>% filter(type %in% c("subset_1_delaunay", 
                                              "subset_1_gabriel", 
                                              "subset_1_relative", 
                                              "subset_1_soi", 
                                              "subset_1_train")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  ylim(0, 3.5) +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_1_gabriel" = "#00BFC4", 
                                "subset_1_relative" = "#00B0F6", 
                                "subset_1_soi" = "#9590FF", 
                                "subset_1_delaunay" = "#E76BF3", 
                                "subset_1_train" = "#FF62BC"),
                     labels = c("Gabriel", 
                                "Relative", 
                                "SOI", 
                                "Delaunay", 
                                "Train"), 
                     name = "Network")
ggsave("plots/Prediction/mase_delaunay_etc_subset_1.pdf", 
       width = 26, height = 13, units = "cm")



# plot MASE for KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_1_overview %>% filter(type %in% c("subset_1_dnn", 
                                              "subset_1_knn", 
                                              "subset_1_queen", 
                                              "subset_1_eco_hub", 
                                              "subset_1_complete")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  ylim(0, 3.5) +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_1_knn" = "#F8766D", 
                                "subset_1_dnn" = "#D89000", 
                                "subset_1_complete" = "#A3A500", 
                                "subset_1_queen" = "#39B600", 
                                "subset_1_eco_hub" = "#00BF7D"),
                     label = c("KNN", 
                               "DNN", 
                               "Complete",
                               "Queen", 
                               "Eco hub"), 
                     name = "Network")
ggsave("plots/Prediction/mase_knn_etc_subset_1.pdf", 
       width = 26, height = 13, units = "cm")




# MASE for data subset 4 --------------------------------------------------
# fit best performing GNAR model for each network for data subset 4

# Queen
mod_4_queen <- fit_and_predict(alpha = 5, 
                               beta = c(4, 0, 0, 0, 0), 
                               net = covid_net_queen_gnar, 
                               vts = datasets_list[[4]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE)


# Eco hub 
mod_4_eco_hub <- fit_and_predict(alpha = 5, 
                                 beta = c(3, 0, 0, 0, 0), 
                                 net = covid_net_eco_hubs_gnar, 
                                 vts = datasets_list[[4]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE, 
                                 numeric_vertices = TRUE, 
                                 county_index = county_index_eco_hubs)

# Railway 
mod_4_train <- fit_and_predict(alpha = 5, 
                               beta = c(5, 1, 0, 0, 0), 
                               net = covid_net_train_gnar, 
                               vts = datasets_list[[4]], 
                               globalalpha = TRUE, 
                               old = TRUE,
                               forecast_window = 5, 
                               return_model = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train)
# Delaunay
mod_4_delaunay <- fit_and_predict(alpha = 5, 
                                  beta = c(2, 1, 0, 0, 0), 
                                  net = covid_net_delaunay_gnar, 
                                  vts = datasets_list[[4]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# Gabriel
mod_4_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(2, 1, 1, 0, 0), 
                                 net = covid_net_gabriel_gnar, 
                                 vts = datasets_list[[4]], 
                                 globalalpha = TRUE, 
                                 old = TRUE,
                                 forecast_window = 5, 
                                 return_model = TRUE)

# Relative 
mod_4_relative <- fit_and_predict(alpha = 5, 
                                  beta = c(5, 1, 0, 0, 0), 
                                  net = covid_net_relative_gnar, 
                                  vts = datasets_list[[4]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# SOI
mod_4_soi <- fit_and_predict(alpha = 5, 
                             beta = c(5, 1, 1, 0, 0), 
                             net = covid_net_soi_gnar, 
                             vts = datasets_list[[4]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# KNN
mod_4_knn <- fit_and_predict(alpha = 5, 
                             beta = c(2, 2, 1, 0, 0), 
                             net = knn_11_gnar, 
                             vts = datasets_list[[4]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# DNN
mod_4_dnn <- fit_and_predict(alpha = 5, 
                             beta = c(1, 1, 1, 0, 0), 
                             net = dnn_225_gnar, 
                             vts = datasets_list[[4]], 
                             globalalpha = TRUE, 
                             old = TRUE,
                             forecast_window = 5, 
                             return_model = TRUE)

# Complete
mod_4_complete <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 1, 0, 0), 
                                  net = complete_net_gnar, 
                                  vts = datasets_list[[4]], 
                                  globalalpha = TRUE, 
                                  old = TRUE,
                                  forecast_window = 5, 
                                  return_model = TRUE)

# compute MASE for best performing GNAR model for each network 
mase_4_queen <- compute_MASE(model = mod_4_queen, 
                             network_name = "subset_4_queen", 
                             n_ahead = 5, 
                             data_df = data_4)

mase_4_eco_hub <- compute_MASE(model = mod_4_eco_hub, 
                               network_name = "subset_4_eco_hub", 
                               n_ahead = 5, 
                               data_df = data_4)

mase_4_train <- compute_MASE(model = mod_4_train, 
                             network_name = "subset_4_train", 
                             n_ahead = 5, 
                             data_df = data_4)

mase_4_delaunay <- compute_MASE(model = mod_4_delaunay, 
                                network_name = "subset_4_delaunay", 
                                n_ahead = 5, 
                                data_df = data_4)

mase_4_gabriel <- compute_MASE(model = mod_4_gabriel, 
                               network_name = "subset_4_gabriel", 
                               n_ahead = 5, 
                               data_df = data_4)

mase_4_relative <- compute_MASE(model = mod_4_relative, 
                                network_name = "subset_4_relative", 
                                n_ahead = 5, 
                                data_df = data_4)

mase_4_soi <- compute_MASE(model = mod_4_soi, 
                           network_name = "subset_4_soi", 
                           n_ahead = 5, 
                           data_df = data_4)

mase_4_knn <- compute_MASE(model = mod_4_knn, 
                           network_name = "subset_4_knn", 
                           n_ahead = 5, 
                           data_df = data_4)

mase_4_dnn <- compute_MASE(model = mod_4_dnn, 
                           network_name = "subset_4_dnn", 
                           n_ahead = 5, 
                           data_df = data_4)

mase_4_complete <- compute_MASE(model = mod_4_complete, 
                                network_name = "subset_4_complete", 
                                n_ahead = 5, 
                                data_df = data_4)


mase_4_overview <- rbind(mase_4_queen, 
                         mase_4_eco_hub, 
                         mase_4_train, 
                         mase_4_delaunay, 
                         mase_4_gabriel, 
                         mase_4_relative, 
                         mase_4_soi, 
                         mase_4_knn, 
                         mase_4_dnn, 
                         mase_4_complete)

# plot MASE for Delaunay, Gabriel, Relative and SOI as well as Railway-based 
# network 
ggplot(mase_4_overview %>% filter(type %in% c("subset_4_delaunay", 
                                              "subset_4_gabriel", 
                                              "subset_4_relative", 
                                              "subset_4_soi", 
                                              "subset_4_train")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  ylim(0, 10) +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_4_gabriel" = "#00BFC4", 
                                "subset_4_relative" = "#00B0F6", 
                                "subset_4_soi" = "#9590FF", 
                                "subset_4_delaunay" = "#E76BF3", 
                                "subset_4_train" = "#FF62BC"),
                     labels = c("Gabriel", 
                                "Relative", 
                                "SOI", 
                                "Delaunay", 
                                "Train"), 
                     name = "Network")
ggsave("plots/Prediction/mase_delaunay_etc_subset_4.pdf", 
       width = 26, height = 13, units = "cm")



# plot MASE for KNN, DNN, Queen, Eco hub, Rail and Complete network
ggplot(mase_4_overview %>% filter(type %in% c("subset_4_dnn", 
                                              "subset_4_knn", 
                                              "subset_4_queen", 
                                              "subset_4_eco_hub", 
                                              "subset_4_complete")), 
       aes(x = time, 
           y = mase, 
           color = type)) +
  geom_point() +
  geom_line(linetype = "dashed") +
  xlab("Time") + 
  ylab("MASE") +
  ylim(0, 10) +
  facet_grid(~ CountyName) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  scale_color_manual(values = c("subset_4_knn" = "#F8766D", 
                                "subset_4_dnn" = "#D89000", 
                                "subset_4_complete" = "#A3A500", 
                                "subset_4_queen" = "#39B600", 
                                "subset_4_eco_hub" = "#00BF7D"),
                     label = c("KNN", 
                               "DNN", 
                               "Complete",
                               "Queen", 
                               "Eco hub"), 
                     name = "Network")
ggsave("plots/Prediction/mase_knn_etc_subset_4.pdf", 
       width = 26, height = 13, units = "cm")




# Save data subsets -------------------------------------------------------
save(datasets_list,
     file = "data/RObjects/data_subsets_regulations.RData")
