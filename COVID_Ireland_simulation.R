# SIMULATION OF DATA AND VERIFICATION OF GNAR MODELS

set.seed(1234)
# Load libraries and objects ----------------------------------------------

library(readr)
library(igraph)
library(GNAR)
library(MASS) # for box cox 
library(tidyverse)
library(magrittr) # for pipping 
library(xtable) # for tables 
library(spdep) # for neighborhood construction 
library(expp)


# Load vectors, GNAR objects and data subsets 
load(file = "data/RObjects/GNAR.RData")
load(file = "data/RObjects/county_index.RData")

load(file = "data/RObjects/population_weight.RData")
load(file = "data/RObjects/distance_urbanisation.RData")
load(file = "data/RObjects/coord_urbanisation.RData")

load(file = "data/RObjects/data_subsets_regulations.RData")

source("functions.R")

# Turn off warnings
options(warn = -1)

# Set ggplot theme 
theme_set(theme_bw(base_size = 16))

# Load data ---------------------------------------------------------------

COVID_weekly_data <- read_csv(file = "data/COVID/ireland_covid_weekly.csv", 
                              show_col_types = FALSE)

# correct format for GNAR
covid_cases <-  COVID_weekly_data %>% 
  dplyr::select(CountyName,
                yw, 
                weeklyCases) %>% 
  spread(CountyName, 
         weeklyCases) %>% 
  column_to_rownames(var = "yw") %>% 
  as.matrix()

# as data frame with time for plotting 
covid_cases_df <- covid_cases %>% 
  as.data.frame() %>% 
  mutate(time = as.Date(rownames(covid_cases))) 


# Fit models --------------------------------------------------------------

# Queen 
model_queen <- fit_and_predict(alpha = 4, 
                               beta = c(2, 2, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_queen_gnar, 
                               old = TRUE, 
                               county_index = county_index_queen, 
                               return_model = TRUE, 
                               forecast_window = 0)

# Eco hubs
model_eco_hubs <- fit_and_predict(alpha = 5, 
                                  beta = c(1, 1, 1, 0, 0),
                                  globalalpha = TRUE, 
                                  net = covid_net_eco_hubs_gnar, 
                                  old = TRUE, 
                                  numeric_vertices = TRUE, 
                                  county_index = county_index_eco_hubs,
                                  return_model = TRUE, 
                                  forecast_window = 0)

# Train 
model_train <- fit_and_predict(alpha = 4, 
                               beta = c(2, 1, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_train_gnar, 
                               old = TRUE, 
                               numeric_vertices = TRUE, 
                               county_index = county_index_train,
                               return_model = TRUE, 
                               forecast_window = 0)

# Delaunay 
model_delaunay <- fit_and_predict(alpha = 4, 
                                  beta = c(2, 2, 1, 1),
                                  globalalpha = TRUE, 
                                  net = covid_net_delaunay_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_delaunay, 
                                  return_model = TRUE, 
                                  forecast_window = 0)

# Gabriel
model_gabriel <- fit_and_predict(alpha = 5, 
                                 beta = c(1, 1, 1, 0, 0),
                                 globalalpha = TRUE, 
                                 net = covid_net_gabriel_gnar, 
                                 old = TRUE, 
                                 county_index = county_index_gabriel,
                                 return_model = TRUE, 
                                 forecast_window = 0)

# Relative 
model_relative <- fit_and_predict(alpha = 4, 
                                  beta = c(2, 2, 1, 1),
                                  globalalpha = TRUE, 
                                  net = covid_net_relative_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_relative,
                                  return_model = TRUE,
                                  forecast_window = 0)

# SOI
model_soi <- fit_and_predict(alpha = 4, 
                             beta = c(2, 1, 1, 1),
                             globalalpha = TRUE, 
                             net = covid_net_soi_gnar, 
                             old = TRUE, 
                             county_index = county_index_soi, 
                             return_model = TRUE, 
                             forecast_window = 0)

# KNN 
model_knn <- fit_and_predict(net = opt_knn_net_gnar, 
                             alpha = 4, 
                             beta = c(1, 1, 0, 0), 
                             globalalpha = TRUE, 
                             old = TRUE, 
                             county_index = county_index_opt_knn, 
                             return_model = TRUE, 
                             forecast_window = 0)

# DNN
model_dnn <- fit_and_predict(net = opt_dnn_net_gnar, 
                             alpha = 4, 
                             beta = c(2, 2, 2, 1), 
                             globalalpha = TRUE, 
                             inverse_distance = TRUE, 
                             county_index = county_index_opt_dnn, 
                             return_model = TRUE, 
                             forecast_window = 0)

# Complete 
model_complete <- fit_and_predict(alpha = 5,
                                  beta = c(1, 1, 1, 1, 0), 
                                  globalalpha = TRUE, 
                                  net = complete_net_gnar, 
                                  old = TRUE, 
                                  county_index = county_index_complete,
                                  return_model = TRUE, 
                                  forecast_window = 0)


# Fit models for subsets -------------------------------------------------
# Queen 
beta_queen <- list(c(3, 2, 2, 1, 1), 
                   c(2, 2, 1, 1, 0), 
                   c(1, 1, 1, 1, 0), 
                   c(4, 0, 0, 0, 0), 
                   c(2, 2, 1, 1, 1))
models_queen <- fit_subset_models(net = covid_net_queen_gnar, 
                                  beta_list = beta_queen, 
                                  alpha_vector = rep(5, 5))





# Average variance --------------------------------------------------------

# average variance across counties to model homogeneity 
av_var_queen <- compute_variance(model_queen)
av_var_eco_hubs <- compute_variance(model_eco_hubs)
av_var_train <- compute_variance(model_train)
av_var_delaunay <- compute_variance(model_delaunay)
av_var_gabriel <- compute_variance(model_gabriel)
av_var_relative <- compute_variance(model_relative)
av_var_soi <- compute_variance(model_soi)
av_var_knn <- compute_variance(model_knn)
av_var_dnn <- compute_variance(model_dnn)
av_var_complete <- compute_variance(model_complete)

# average variance for data subsets 
var_queen <- compute_variance_subsets(models_queen)


# Simulation --------------------------------------------------------------
# all counties for which to simulate data 
counties_v <- covid_cases %>% colnames()


# homogeneity in variance across counties modeled for simplification 

# data simulation 
simulation_queen <- simulate_time_series(gnar_object = covid_net_queen_gnar, 
                                         GNARmodel = model_queen, 
                                         beta_order =  c(2, 2, 1, 1), 
                                         av_var = av_var_queen, 
                                         county_index = county_index_queen)

simulation_eco_hubs <- simulate_time_series(gnar_object = covid_net_eco_hubs_gnar, 
                                         GNARmodel = model_eco_hubs, 
                                         beta_order =  c(1, 1, 1, 0, 0), 
                                         av_var = av_var_eco_hubs, 
                                         county_index = county_index_eco_hubs)

simulation_train <- simulate_time_series(gnar_object = covid_net_train_gnar, 
                                         GNARmodel = model_train, 
                                         beta_order =  c(2, 1, 1, 1), 
                                         av_var = av_var_train, 
                                         county_index = county_index_train)

simulation_delaunay <- simulate_time_series(gnar_object = covid_net_delaunay_gnar, 
                                            GNARmodel = model_delaunay, 
                                            beta_order =  c(2, 2, 1, 1), 
                                            av_var = av_var_delaunay, 
                                            county_index = county_index_delaunay)

simulation_gabriel <- simulate_time_series(gnar_object = covid_net_gabriel_gnar, 
                                            GNARmodel = model_gabriel, 
                                            beta_order =  c(1, 1, 1, 0, 0), 
                                            av_var = av_var_gabriel, 
                                            county_index = county_index_gabriel)

simulation_relative <- simulate_time_series(gnar_object = covid_net_relative_gnar, 
                                            GNARmodel = model_relative, 
                                            beta_order =  c(2, 2, 1, 1), 
                                            av_var = av_var_relative, 
                                            county_index = county_index_relative)

simulation_soi <- simulate_time_series(gnar_object = covid_net_soi_gnar, 
                                            GNARmodel = model_soi, 
                                            beta_order =  c(2, 1, 1, 1), 
                                            av_var = av_var_soi, 
                                            county_index = county_index_soi)

simulation_knn <- simulate_time_series(gnar_object = opt_knn_net_gnar, 
                                            GNARmodel = model_knn, 
                                            beta_order =  c(1, 1, 0, 0), 
                                            av_var = av_var_knn, 
                                            county_index = county_index_opt_knn)

simulation_dnn <- simulate_time_series(gnar_object = opt_dnn_net_gnar, 
                                            GNARmodel = model_dnn, 
                                            beta_order =  c(2, 2, 2, 1), 
                                            av_var = av_var_dnn, 
                                            county_index = county_index_opt_dnn)

simulation_complete <- simulate_time_series(gnar_object = complete_net_gnar, 
                                            GNARmodel = model_complete, 
                                            beta_order =  c(1, 1, 1, 1, 0), 
                                            av_var = av_var_complete, 
                                            county_index = county_index_complete)


# Simulation for data sets ------------------------------------------------

simulation_subsets_queen <- simulate_time_series_subsets(gnar_object = covid_net_queen_gnar, 
                                                         GNARmodel_list = models_queen, 
                                                         beta_list = beta_queen, 
                                                         var_vector = var_queen, 
                                                         county_index = county_index_queen)



ggplot(simulation_subsets_queen %>% filter(CountyName == "Dublin"), 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = subset %>% as.factor())) +
  geom_line()

# Visualization -----------------------------------------------------------

# plot simulated data 
ggplot(simulation_queen, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/queen_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_eco_hubs, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/eco_hubs_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_train, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom") 
ggsave("plots/simulations/train_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_delaunay, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/delaunay_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_gabriel, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/gabriel_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_relative, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/relative_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")

ggplot(simulation_soi, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/soi_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_knn, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/knn_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_dnn, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/dnn_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


ggplot(simulation_complete, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/simulations/complete_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


# Reconstruct coefficients ------------------------------------------------
coef_queen <- model_queen %>% coef() %>% as.data.frame()
colnames(coef_queen) <- "real_param"

simulation_queen_short <- simulation_queen %>% 
  spread(CountyName, ID) %>% 
  dplyr::select(-time) %>% 
  as.matrix()

gnar_model_queen <- GNARfit(vts =  simulation_queen_short[1:20, ], 
        net = covid_net_queen_gnar, 
        alphaOrder = 4, 
        betaOrder = c(2, 2, 1, 1), 
        globalalpha = TRUE)

coef_queen$recomp_param <- gnar_model_queen %>% coef()
coef_queen$type <- rownames(coef_queen)

strCaption <- "Simulation parameter values and re-computed parameter values for \\code{GNAR} model for \\textbf{Delaunay triangulation} network"
print(xtable(coef_queen[, c(3, 1, 2)],
             digits=2,
             caption=strCaption,
             label="tab:coef_queen", 
             align = c("", "l", "|", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(coef_queen)),
                        command = c(paste("\\toprule \n",
                                          "Parameter & real value & re-computed value \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)








