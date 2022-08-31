## SIMULATE NETWORKS 

library(readr)
library(tidyverse)

library(igraph)

source("functions.R")

set.seed(1234)

# Load pre-processed data -------------------------------------------------
COVID_weekly_data <- read_csv(file = "data/COVID/ireland_covid_weekly.csv", 
                              show_col_types = FALSE)


# Erdös-Renyi -------------------------------------------------------------

# no structure, edges at random not based on any prior knowledge / geography 

# number of vertices 
N <- COVID_weekly_data %>% 
  dplyr::select(CountyName) %>% 
  unique() %>% 
  nrow()

# one dense 
p_dense <- 0.95

# one sparse
p_sparse <- 0.21

# 1. Dense network
# reproducible through seed values saved in seed.nos
net1 <- seedToNet(seed.no = seed.nos[1], 
                  nnodes = N, 
                  graph.prob = p_dense)

net1_igraph <- GNARtoigraph(net1)

network_characteristics(net1_igraph)

# very comparable to best KNN graph 
# metric     values
# 1           av. degree 23.5384615
# 2              density  0.9415385
# 3      av. path length  1.0584615
# 4    global clustering  0.9463203
# 5 av. local clustering  0.9479992
# 6      av. betweenness  0.7307692
# 7     s.d. betweenness  0.2665347

pdf("plots/simulations/ER_covid.pdf")
plot(net1, 
     layout=layout_nicely, 
     vertex.label = unique(COVID_weekly_data$CountyName),
     vertex.label.font = 1, 
     rescale = TRUE, 
     vertex.label.color = "#242526", 
     vertex.color = "#FF7F7F")
dev.off()



# all models to investigate
alpha_options <- seq(1, 5)
beta_options <- list(0, 1, 
                    c(1, 0), 
                    c(1, 1), 
                    c(2, 0), 
                    c(2, 1), 
                    c(2, 2), 
                    c(1, 0, 0), 
                    c(1, 1, 0), 
                    c(1, 1, 1), 
                    c(2, 1, 1),
                    c(1, 0, 0, 0), 
                    c(1, 1, 0, 0), 
                    c(1, 1, 1, 0), 
                    c(1, 1, 1, 1), 
                    c(2, 1, 1, 1),
                    c(2, 2, 1, 1),
                    c(1, 0, 0, 0, 0), 
                    c(1, 1, 0, 0, 0), 
                    c(1, 1, 1, 0, 0), 
                    c(1, 1, 1, 1, 0), 
                    c(1, 1, 1, 1, 1),
                    c(2, 1, 1, 1, 1), 
                    c(2, 2, 1, 1, 1))
globalalpha <- c("TRUE", "FALSE")
# create all possible combinations of model parameters 
model_options <-  expand.grid(alpha_options, 
                              beta_options, 
                              globalalpha)

# filter out valid parameter combinations 
model_options$valid <- (model_options$Var1 == model_options$Var2 %>% 
                          lapply(FUN = function(i) {length(i)}) %>% 
                          unlist())
model_options_valid <- model_options %>% filter(valid == TRUE)

# extract alphas an betas 
alphas <- model_options_valid$Var1
betas <- model_options_valid$Var2


# simulate ER graph for each model n_rep times with global alpha and
# non-global alpha  
n_rep <- 10 

sim_res <- simNetER_pre_rep(alphas = alphas, 
                            betas = betas, 
                            data = covid_cases, 
                            n_rep = n_rep, 
                            ER_p = p_dense)

sim_res_df <- do.call(cbind.data.frame, sim_res) %>% 
  gather(key = "model", 
         value = "pred_error")


# prediction error boxplots for population proportional cases
ggplot(data = sim_res_df,
       aes(x = model, 
           y = pred_error)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1)) +
  xlab("") +
  ylab("prediction error")
ggsave("plots/simulations/covid_pop_cases_pred_error_boxplots.pdf")





