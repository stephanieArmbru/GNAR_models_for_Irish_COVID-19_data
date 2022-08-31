# SIMULATION OF DATA AND VERIFICATION OF GNAR MODELS

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
library(spdep) # for neighborhood construction 
library(expp)


# load vectors, GNAR objects and data subsets 
load(file = "data/RObjects/GNAR.RData")
load(file = "data/RObjects/county_index.RData")

load(file = "data/RObjects/population_weight.RData")
load(file = "data/RObjects/distance_urbanisation.RData")
load(file = "data/RObjects/coord_urbanisation.RData")

load(file = "data/RObjects/data_subsets_regulations.RData")

# load functions 
source("functions.R")

# turn off warnings
options(warn = -1)

# set ggplot theme 
theme_set(theme_bw(base_size = 16))


# Load data ---------------------------------------------------------------
# load pre-processed data
COVID_weekly_data <- read_csv(file = "data/COVID/ireland_covid_weekly.csv", 
                              show_col_types = FALSE)

# transform into correct format for GNAR
covid_cases <-  COVID_weekly_data %>% 
  dplyr::select(CountyName,
                yw, 
                weeklyCases) %>% 
  spread(CountyName, 
         weeklyCases) %>% 
  column_to_rownames(var = "yw") %>% 
  as.matrix()

# transform into data frame with time column
covid_cases_df <- covid_cases %>% 
  as.data.frame() %>% 
  mutate(time = as.Date(rownames(covid_cases))) 


# Fit models --------------------------------------------------------------
# fit best performing GNAR model for each network and return GNAR model 
# Queen 
model_queen <- fit_and_predict(alpha = 4, 
                               beta = c(2, 2, 1, 1),
                               globalalpha = TRUE, 
                               net = covid_net_queen_gnar, 
                               old = TRUE, 
                               county_index = county_index_queen, 
                               return_model = TRUE, 
                               forecast_window = 0)

# Fit models for subsets -------------------------------------------------
# fit best performing GNAR model for each data subsets and return model 
# for each subset 
# exemplarily for Queen's contiguity network 

# Queen 
alpha_queen <- rep(5, 5)
beta_queen <- list(c(3, 2, 2, 1, 1), 
                   c(2, 2, 1, 1, 0), 
                   c(1, 1, 1, 1, 0), 
                   c(4, 0, 0, 0, 0), 
                   c(2, 2, 1, 1, 1))
models_queen <- fit_subset_models(net = covid_net_queen_gnar, 
                                  beta_list = beta_queen, 
                                  alpha_vector = alpha_queen)

# Average variance --------------------------------------------------------
# compute average variance across counties to model homogeneity 
av_var_queen <- compute_variance(model_queen)

# compute average variance for data subsets individually 
var_queen <- compute_variance_subsets(models_queen)


# Simulation --------------------------------------------------------------
# generate vector containing all counties for which to simulate data 
counties_v <- covid_cases %>% colnames()

# homogeneity in variance across counties modeled for simplification 
simulation_queen <- simulate_time_series(gnar_object = covid_net_queen_gnar, 
                                         GNARmodel = model_queen, 
                                         beta_order =  c(2, 2, 1, 1), 
                                         av_var = av_var_queen, 
                                         county_index = county_index_queen, 
                                         timeframe = 120)

# plot simulated data and true 1-lag COVID-19 ID for each network to assess how 
# well the true data is predicted 

# Queen
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
ggsave("plots/simulations/full_dataset/queen_20.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")

# Simulation for data sets ------------------------------------------------
# simulate data according to variance in data subsets
# 5 simulations with different variance 
# time period of simulation equal to the length of the corresponding data subset 

# exemplarily for Queen's contiguity network 
simulation_subsets_queen <- simulate_time_series_subsets(gnar_object = covid_net_queen_gnar, 
                                                         GNARmodel_list = models_queen, 
                                                         beta_list = beta_queen, 
                                                         var_vector = var_queen, 
                                                         county_index = county_index_queen)


# visualise simulate data for all counties 
ggplot(simulation_subsets_queen, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Subset")) 

ggsave("plots/simulations/subsets/queen_all_counties.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")


# Example: Queen + Dublin -------------------------------------------------
# for GNAR model fit on entire data set
# simulated data 
simulation_queen_dublin <- simulation_queen %>% 
  filter(CountyName == "Dublin")

# real 1-lag COVID-19 ID
simulation_queen_dublin$real <- COVID_weekly_data %>% 
  filter(CountyName == "Dublin") %>% 
  pull(weeklyCases)

# visualise predictive performance 
ggplot(simulation_queen_dublin, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  geom_line(aes(x = time, 
                y = real), 
            color = "grey", 
            linetype = "dashed") +
xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  theme(legend.position = "none")
ggsave("plots/simulations/full_dataset/queen_dublin_120.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")

# for GNAR models fit on data subsets 
# simulated data
simulation_subsets_queen_dublin <- simulation_subsets_queen %>% 
  filter(CountyName == "Dublin")

# real 1-lag COVID-19 ID 
simulation_subsets_queen_dublin$real <- COVID_weekly_data %>% 
  filter(CountyName == "Dublin") %>% 
  pull(weeklyCases)

# visualise predictive performance 
ggplot(simulation_subsets_queen_dublin, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = subset %>% as.factor())) +
  geom_line() +
  geom_line(aes(x = time, 
                y = real), 
            color = "grey", 
            linetype = "dashed") +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Subset")) 
ggsave("plots/simulations/subsets/queen_dublin.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")

# Reconstruct coefficients ------------------------------------------------
# re-fit best performing GNAR model based on the simulation data
# compare original coefficient estimates which determine the simulation 
# with their new estimate to assess how well GNAR models detect underlying
# temporal and spatial dependencies 

# recompute GNAR model coefficients 
coef_queen <- reconstruct_coefficients(model_queen, 
                                       simulation_queen, 
                                       alphaOrder = 4, 
                                       betaOrder = c(2, 2, 1, 1), 
                                       covid_net_queen_gnar)

# for latex 
strCaption <- "Coefficients for the \\code{GNAR} model serving the data simulation 
and their re-computed values after fitting the \\code{GNAR} model to the simulated data 
and their 95\\% confidence interval (CI) for the \\textbf{Queen's contiguity} network"
print(xtable(coef_queen,
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
                                          "Coefficient & real value & re-computed value [95\\% CI] \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)

# Reconstruct coefficients for subsets ------------------------------------
# re-fit best performing GNAR model based on the simulation data
# compare original coefficient estimates which determine the simulation 
# with their new estimate to assess how well GNAR models detect underlying
# temporal and spatial dependencies 

# extract coefficients for best performing GNAR model for each coefficients 
real_coef_queen <- extract_coef(models_queen)

# re-compute GNAR model coefficients
coef_subsets_queen <- reconstruct_coefficients_subsets(GNAR_model_list = models_queen, 
                                                       simulation_subsets_df = simulation_subsets_queen, 
                                                       alphaOrder_vector = alpha_queen, 
                                                       betaOrder_vector = beta_queen, 
                                                       net = covid_net_queen_gnar)

# for latex
strCaption <- "Coefficients for the \\code{GNAR} model serving the data subset simulation 
and their re-computed values after fitting the \\code{GNAR} model to the simulated data subset
and their 95\\% confidence interval (CI) for the \\textbf{Queen's contiguity} network"
print(xtable(coef_subsets_queen,
             digits=2,
             caption=strCaption,
             label="tab:coef_queen", 
             align = c("", "l", "|", "r", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(coef_subsets_queen)),
                        command = c(paste("\\toprule \n",
                                          "Subset & coefficient & real value & re-computed value [95\\% CI] \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)



# Simulation Test ---------------------------------------------------------
# assess if bad performance of GNAR model in estimating correct coefficients
# can be traced back to high variance of the error term by simulating data 
# with a standard normal error (1) and re-computing model coefficients 

# exemplarily for Queen's contiguity network

# simulate 120 weeks 
simulation_test_queen <- simulate_time_series(gnar_object = covid_net_queen_gnar, 
                                              GNARmodel = model_queen, 
                                              beta_order =  c(2, 2, 1, 1), 
                                              av_var = 1, 
                                              county_index = county_index_queen, 
                                              timeframe = 120)

# re-compute GNA model coefficients 
coef_test_queen <- reconstruct_coefficients(model_queen, 
                                       simulation_test_queen, 
                                       alphaOrder = 4, 
                                       betaOrder = c(2, 2, 1, 1), 
                                       covid_net_queen_gnar)

# for latex 
strCaption <- "Coefficients for the \\code{GNAR} model serving the data simulation 
and their re-computed values after fitting the \\code{GNAR} model to the simulated data 
and their 95\\% confidence interval (CI) for the \\textbf{Queen's contiguity} network; 
the variance is set to 1 for testing."
print(xtable(coef_test_queen,
             digits=2,
             caption=strCaption,
             label="tab:coef_testing_queen", 
             align = c("", "l", "|", "r", "r")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(coef_test_queen)),
                        command = c(paste("\\toprule \n",
                                          "Coefficient & real value & re-computed value [95\\% CI] \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)


# visualise predictive performance for Dublin 

# simulated data 
simulation_test_queen_dublin <- simulation_test_queen %>% 
  filter(CountyName == "Dublin")

# true data
simulation_test_queen_dublin$real <- COVID_weekly_data %>% 
  filter(CountyName == "Dublin") %>% 
  pull(weeklyCases)


# visualise predictive performance 
ggplot(simulation_test_queen_dublin, 
       aes(x = time, 
           y = ID, 
           group = CountyName, 
           color = CountyName)) +
  geom_line() +
  geom_line(aes(x = time, 
                y = real), 
            color = "grey", 
            linetype = "dashed") +
  xlab("Time") +
  ylab("Simulated COVID-19 ID") +
  theme(legend.position = "none")
ggsave("plots/simulations/simulation_test_dublin.pdf", 
       width = 26, 
       height = 14, 
       unit = "cm")

