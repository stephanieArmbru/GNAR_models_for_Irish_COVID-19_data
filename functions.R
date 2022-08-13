### FUNCTIONS


# Load libraries ----------------------------------------------------------
library(ade4) # igraph to neighbourhood list object
library(Hmisc) # for weighted variance 
library(Metrics) # for MASE computation 

# source adapted functions from GNAR package 
source("GNAR/GNARdesign.R")
source("GNAR/GNARfit.R")
source("GNAR/NofNeighbours.R")

# Dataframe functions -----------------------------------------------------
# find maximum value for each column
colMax <- function(data) {
  return(sapply(data, max, na.rm = TRUE))
}

# find minimum value for each column 
colMin <- function(data) {
  return(sapply(data, min, na.rm = TRUE))
}


# Maps COVID --------------------------------------------------------------
# plot county map of Ireland, coloured indicating the COVID-19 incidence
map_covid <- function(desired_date # form "2020-03-01", has to be Monday
) {
  # filter data for desired date
  COVID_week<- COVID_weekly_data %>% 
    filter(yw == desired_date) %>% 
    dplyr::select(weeklyCases, CountyName) 
  
  ireland$COVID <- COVID_week[match(ireland$NAME_1, 
                                    COVID_week$CountyName), ] %>% 
    pull(weeklyCases)
  
  # define colour palette 
  pal <- colorNumeric(palette ="YlGnBu", 
                      domain = ireland$COVID)
  
  # generate plot 
  map <- leaflet(ireland) %>%  
    addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
                opacity = 1.0, fillOpacity = 0.5,
                fillColor = ~pal(COVID),
                highlightOptions = highlightOptions(color = "white", weight =2,
                                                    bringToFront = TRUE),
                label = ~NAME_1) %>% 
    addLegend("bottomright", 
              pal = pal, 
              values = ~COVID,
              title = paste0("COVID cases on ", desired_date),
              opacity = 1
    )
  
  return(map)
}


# Maps mumps --------------------------------------------------------------
# plot county map for England and Wales, coloured according to mumps incidence
map_mumps <- function(desired_week # integer
)  {
  # filter data for desired week
  mumps_week <- mumps_df[desired_week, ] %>% gather("County", "mumps")
  
  # assign mumps cases in desired week to counties 
  england_wales$MUMPS <- mumps_week[match(england_wales$NAME_2, mumps_week$County), ] %>% 
    pull(mumps)
  
  # assign Wales data 
  england_wales$MUMPS[england_wales$NAME_1 == "Wales"] <- mumps_week$mumps[mumps_week$County == "Wales"]
  
  # define colour palette 
  pal <- colorNumeric(palette ="YlGnBu", 
                      domain = england_wales$MUMPS)
  
  # generate map 
  map <- leaflet(england_wales) %>%  
    addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
                opacity = 1.0, fillOpacity = 0.5,
                fillColor = ~pal(MUMPS),
                highlightOptions = highlightOptions(color = "white", weight =2,
                                                    bringToFront = TRUE),
                label = ~NAME_2) %>% 
    addLegend("bottomright", 
              pal = pal, 
              values = ~MUMPS,
              title = paste0("Mumps cases (2005, week ", desired_week, ")"),
              opacity = 1
    )
  
  return(map)
}


# Network characteristics -------------------------------------------------
# compute network characteristics
network_characteristics <- function(igraph_obj, 
                                    network_name) {
  
  density <- igraph_obj %>% graph.density() 
  apl <- igraph_obj %>%  average.path.length(directed = FALSE) 
  
  global_clust <- igraph_obj %>% transitivity(type = "global") 
  mean_local_clust <- igraph_obj %>% 
    transitivity(type = "local",
                 isolates = "zero") %>% 
    mean()

  
  degree_v <- igraph_obj %>% degree()
  av_degree <- degree_v %>% mean() 
  
  
  max_degree <- degree_v %>% max() 
  which(degree_v == max_degree) 
  
  min_degree <- degree_v %>% min() 
  which(degree_v == min_degree) 
  
  # betweenness
  bet <- betweenness(igraph_obj, 
                     v=V(igraph_obj), 
                     directed = FALSE)
  
  min_bet <- bet %>% min() 
  bet[which(bet == min_bet)] 
  
  max_bet <- bet %>% max()
  bet[which(bet == max_bet)] 
  
  
  graph_char <- data.frame("metric" = c("av. degree", 
                                        "density", 
                                        "av. SPL", 
                                        "global clust.", 
                                        "av. local clust.", 
                                        "av. betw.", 
                                        "s.d. betw."), 
                           "values" = c(av_degree, 
                                        density, 
                                        apl, 
                                        global_clust, 
                                        mean_local_clust, 
                                        mean(bet), 
                                        sd(bet))) 
  colnames(graph_char) <- c("metric", network_name)
  return(graph_char)
}

# check if network is scale-free via log-log plot and R squared 
is_scale_free <- function(igraph_net, # igraph object required
                          network_name, 
                          Newman = TRUE) { # less-noise tail for formula acc. to Newman 2005
  
  
  if (!Newman) {
    df <- data.frame("X" = igraph_net %>% degree() %>% log()) %>% 
      group_by(X) %>% 
      mutate(Y = log(n() / gorder(igraph_net))) %>% 
      ungroup()
    
    name_y <- "log. emp. distribution"
  }
  
  if (Newman) {
    df <- data.frame("X" = igraph_net %>% 
                            degree() %>% 
                            unique() %>% 
                            log())
    
    ecdf_fun <- igraph_net %>% 
      degree() %>% 
      ecdf()
    
    df$Y <- igraph_net %>% 
               degree() %>% 
               unique() %>% 
               ecdf_fun() %>% log()
    
    name_y <- "log. emp. cum. distribution"
  }
  
  
  g <- ggplot(df, 
              aes(x = X, y = Y)) +
    geom_point() +
    geom_smooth(method = 'lm', 
                formula = y~x, 
                se = FALSE) +
    xlab("log(d)") +
    ylab(name_y)
  
  ggsave(file = paste0("plots/scaleFree/log_log_", network_name, ".pdf"), 
         plot = g, 
         width = 13, 
         height = 13, 
         unit = "cm")
  
  lin_reg_mod <- lm(Y ~ X, data = df)
  s <- summary(lin_reg_mod)
  
  # return plot and R squared 
  return(list("graph" = g, 
              "R_squared" = s$r.squared))
}


# Correlation and Performance ---------------------------------------------
# compute Residual Sum of Squares
RSS <- function(GNARfit_object) {
  return(sum(residuals(GNARfit_object)^2))
}

# compute Moran's I
moran_I <- function(data = COVID_weekly_data,
                    coords = coord_urbanisation, 
                    nb_list, 
                    inverse_distance = FALSE) {
  moran_list <- list()
  
  # compute Great Circle distances 
  geoms <- st_as_sf(coords %>% as.data.frame(), 
                    coords = c("X", "Y"), 
                    remove = FALSE)
  
  if (!inverse_distance) {
    dist_weights <- nb2listwdist(nb_list, 
                                 as(geoms,
                                    "Spatial"), 
                                 type = "idw", 
                                 alpha = 1, # inverse distance
                                 style = "C",
                                 longlat = TRUE)
  }
  if (inverse_distance) {
    dist_weights <- nb2listw(nb_list, 
                             style = "C")
  }
  
  # for each date, compute Moran's I
  for (date in data$yw %>% unique() %>% as.character()) {
    moran_list[[date]] <- moran(data[data$yw == date, ]$weeklyCases,
                                listw = dist_weights, # row-wise normalisation
                                n = length(nb_list), 
                                Szero(dist_weights) # global sum of weights
                                )$I
  }
  
  moran_df <- do.call(rbind.data.frame, moran_list)
  colnames(moran_df) <- "moran"
  moran_df$dates <- as.Date(names(moran_list))
  
  return(moran_df)
}


# Classical decomposition - R ---------------------------------------------
# perform decomposition into trend-season-random component
tsr_decomposition <- function(county_name, 
                              column_name = "weeklyCases",
                              type_decomposition = "additive", 
                              data = COVID_data_week_smooth) {
  covid_county <- data %>%
    filter(CountyName == county_name) %>% 
    dplyr::select(yw, 
                  column_name) 
  
  covid_county_ts <- covid_county %>% 
    dplyr::select(column_name) %>% 
    ts(frequency = 52, start = c(2020, 1))
  
  COVID_ts_decomposed <- stats::decompose(covid_county_ts, 
                                          type = type_decomposition)
  
  return(COVID_ts_decomposed)
}


# Visualise classical decomposition - manual ------------------------------
# plot trend component 
plot_trend <- function(county_name, 
                       data = COVID_data_week_smooth) {
  ggplot(data %>% filter(CountyName == county_name), 
         aes(x = yw, 
             y = weeklyCases, 
             color = "weekly cases")) +
    geom_line(linetype = "dashed") +
    geom_line(aes(x = yw, 
                  y = trendCases, 
                  color = "trend"),
              size = 0.8) +
    xlab("time") +
    ylab("COVID cases") +
    ggtitle(county_name) +
    theme(legend.title = element_blank(), 
          legend.position = "bottom") +
    scale_color_manual(values = c("weekly cases" = "#096192", 
                                  "trend" = "#24aae2"))
  
  file_name <- paste0("plots/stationarity/trend_", county_name, ".pdf")
  ggsave(file_name, 
         width = 23, height = 14, unit = "cm")
}

# plot seasonal effect 
plot_season <- function(county_name,
                        data = COVID_data_week_smooth) {
  ggplot(data %>% filter(CountyName == county_name), 
         aes(x = yw, 
             y = weeklyCases, 
             color = "weekly cases")) +
    geom_line(linetype = "dashed") +
    geom_line(aes(x = yw, 
                  y = seasonCases, 
                  color = "seasonality"), 
              size = 0.8) +
    xlab("time") +
    ylab("COVID cases") +
    ggtitle(county_name) +
    theme(legend.title = element_blank(), 
          legend.position = "bottom") +
    scale_color_manual(values = c("weekly cases" = "#096192", 
                                  "seasonality" = "#0e8c7f"))
  
  file_name <- paste0("plots/stationarity/season_", county_name, ".pdf")
  ggsave(file_name, 
         width = 23, height = 14, unit = "cm")
}

# plot residual / random error 
plot_residuals <- function(county_name, 
                           data = COVID_data_week_smooth) {
  ggplot(data %>% filter(CountyName == county_name), 
         aes(x = yw, 
             y = weeklyCases, 
             color = "weekly cases")) +
    geom_line(linetype = "dashed") +
    geom_line(aes(x = yw, 
                  y = randomError, 
                  color = "residual"), 
              size = 0.8) +
    xlab("time") +
    ylab("Weekly COVID cases") +
    ggtitle(county_name) + 
    theme(legend.title = element_blank(), 
          legend.position = "bottom") +
    scale_color_manual(values = c("weekly cases" = "#096192", 
                                  "residual" = "black"))
  
  
  file_name <- paste0("plots/stationarity/residual_", county_name, ".pdf")
  ggsave(file_name, 
         width = 23, height = 14, unit = "cm")
}

# plot entire decomposition 
plot_tsr <- function(county_name, 
                     data = COVID_data_week_smooth) {
  ggplot(data %>% filter(CountyName == county_name), 
         aes(x = yw, 
             y = weeklyCases, 
             color = "weekly cases")) +
    geom_line(linetype = "dashed") +
    # trend 
    geom_line(aes(x = yw, 
                  y = trendCases, 
                  color = "trend")) + 
    # season 
    geom_line(aes(x = yw, 
                  y = seasonCases, 
                  color = "seasonality")) +
    # residual 
    geom_line(aes(x = yw, 
                  y = randomError, 
                  color = "residual")) +
    xlab("time") +
    ylab("Weekly COVID cases") +
    theme(legend.title = element_blank(), 
          legend.position = "bottom") +
    scale_color_manual(values = c("weekly cases" = "#096192", 
                                  "trend" = "#24aae2", 
                                  "seasonality" = "#0e8c7f", 
                                  "residual" = "black"))
  
  
  file_name <- paste0("plots/stationarity/all_", county_name, ".pdf")
  ggsave(file_name, 
         width = 27, height = 14, unit = "cm")
}


# Shortest path length ----------------------------------------------------
# compute and plot shortest path length between all vertices based on igraph 
# object and the difference in 1-lag COVID-19 ID between vertex pairs 
# average the difference in 1-lag COVID-19 ID over time 
spl_and_mean_difference <- function(igraph_object, 
                                    county_index, 
                                    numeric_vertices = FALSE, 
                                    network_name) {
  
  spl_df <- data.frame(mean_diff = NA,
                       CountyName = NA, 
                       spl = NA)
  
  for (county in counties) {
    # current county we are looking at 
    cases_county <- covid_cases[, county]
    
    # difference in 1-lag COVID-19 ID with all the other counties 
    diff_cases <- covid_cases[, counties != county] - cases_county
    diff_cases <- diff_cases %>% as.data.frame()
    
    # compute average difference across time 
    diff_cases_mean <- data.frame("mean_diff" = diff_cases %>% colMeans())
    
    diff_cases_mean$CountyName <- rownames(diff_cases_mean)
    diff_cases_mean$spl <- NA
    
    # obtain vertex id number for county in question 
    county_to_vertex <- county_index[county_index$CountyName == county, ]$index
    
    # compute Great Circle distance between county in questions and all other 
    # counties 
    for (other_vertex in seq(1, 26)[-county_to_vertex]) {
      dist <- distances(graph = igraph_object,
                        v = county_to_vertex, 
                        to = other_vertex)
      
      if (numeric_vertices) {
        vertex_to_county <- county_index[county_index$index == other_vertex, ]$CountyName
      } else {
        vertex_to_county <- colnames(dist)
      }
      diff_cases_mean[diff_cases_mean$CountyName == vertex_to_county, ]$spl <- dist
    }
    spl_df <- rbind(spl_df, diff_cases_mean)
  }
  # plot average difference against shortest path length 
  # filter out positive average differences to avoid double counting 
  ggplot(spl_df %>% filter(mean_diff >= 0),
         aes(x = spl,
             y = mean_diff)) +
    geom_point() +
    xlab("Shortest path length") +
    ylab("Average difference in COVID-19 ID")
  ggsave(paste0("plots/spatialCor/shortest_path_length_mean_covid_diff_", 
                network_name, ".pdf"), 
         width = 14, height = 14, unit = "cm")
  
  return(spl_df)
}

# Great Circle distance  --------------------------------------------------
# compute Great Circle distance between all counties 
# requires data as n x 2 matrix with rownames, columns are longitude and 
# latitude coordinates 
circle_distance <- function(data) { 
  pairwise_dist <- matrix(nrow = 26, ncol = 26)
  for (i in seq(1, dim(data)[1])) {
    for (j in seq(1, dim(data)[1])) {
      long1 <- data[i, 1]
      long2 <- data[j, 1]
      lat1 <- data[i, 2]
      lat2 <- data[j, 2]
      pairwise_dist[i, j] <- distm(x = c(long1, lat1), 
                                   y = c(long2, lat2), 
                                   fun = distHaversine)
    }
  }
  # transform into data frame 
  dist_df <- pairwise_dist %>% as.data.frame()
  rownames(dist_df) <- rownames(data)
  colnames(dist_df) <- rownames(data)
  
  # substitute zero diagonals with NA  
  dist_df[dist_df == 0] <- NA
  
  return(dist_df)
}


# Epidemiology ------------------------------------------------------------
# compute reproduction number according to Pastor-Satorras and Vespignani 2002 
reproduction_rate_psv <- function(igraph_object, 
                                  beta = beta_covid, 
                                  gamma = gamma_covid) {
  # degree distribution 
  degree_vector <- igraph_object %>% degree()
  
  return(beta / gamma * (mean(degree_vector) + var(degree_vector) / mean(degree_vector)))
}


# compute reproduction number according to Lloyd and Valeika 2007 
# adapted to include weighted mean and variance  
reproduction_rate_lv <- function(igraph_object, 
                                 beta = beta_covid, 
                                 gamma = gamma_covid, 
                                 
                                 # if population_weighting, weighted mean and variance 
                                 # computed to account for population size of counties
                                 weight_df = population_weight, 
                                 
                                 numeric_vertices = FALSE, 
                                 county_index = NULL) {  
  
  # degree distribution 
  degree_vector <- igraph_object %>% degree()
  
  # weighting according to population size 
  if (!numeric_vertices) {
    weights <- weight_df[match(names(degree_vector), weight_df$CountyName), ]$weight
  } else  {
    weights <- weight_df[match(county_index$CountyName, weight_df$CountyName), ]$weight
  }
  
  mean(degree_vector) - 1 + mean(degree_vector) / var(degree_vector)
  weighted_mean_degree <- weighted.mean(degree_vector, weights)
  weighted_var_degree <- wtd.var(degree_vector, weights)
  
  return(beta / (gamma + beta) * (weighted_mean_degree - 1 + weighted_var_degree / weighted_mean_degree))
}


# compute reproduction number according to slides 
reproduction_rate_t <- function(igraph_object, 
                                beta = beta_covid, 
                                gamma = gamma_covid) {
  # degree distribution 
  degree_vector <- igraph_object %>% degree()
  degree_vector_unique <- degree_vector %>% unique()
  
  edf_degree <- table(degree_vector) / length(degree_vector)
  
  # first and second derivative of probability generating function for empirical probability density 
  g_2 <- degree_vector_unique * (degree_vector_unique - 1) * edf_degree[match(degree_vector_unique, 
                                                                              names(edf_degree))]
  g_1 <- degree_vector_unique * edf_degree[match(degree_vector_unique, 
                                                 names(edf_degree))]
  
  return(beta / (gamma + beta) * (sum(g_2) / sum(g_1)))
}

# assess the sensitivity of the reproduction number for random attacks
# with absolute number of edge deletions 
reproduction_sensitivity <- function(igraph_object, 
                                     n_deletions = seq(1, 5), 
                                     n_rep = 100, 
                                     numeric_vertices = FALSE, 
                                     county_index = NULL) {
  
  reproduction_list <- data.frame()
  
  for(i in n_deletions) {
    reproduction_vector <- c()
    
    for (n in seq(1, n_rep)) {
      # randomly delete a certain number of edges 
      random_edge_sample <- sample(igraph_object %>% E(),
                                   size = i)
      
      igraph_attacked <- igraph_object %>% 
        delete_edges(edges = random_edge_sample)
      
      # re-compute reproducton number 
      reproduction_vector <- c(reproduction_vector, 
                               reproduction_rate_lv(igraph_attacked, 
                                                    numeric_vertices = numeric_vertices, 
                                                    county_index = county_index))
    }
    reproduction_add <- list(reproduction_vector)
    names(reproduction_add) <- i 
    reproduction_list <- c(reproduction_list, reproduction_add)
  }
  
  reproduction_df <- do.call(cbind.data.frame, reproduction_list)
  colnames(reproduction_df) <- n_deletions
  
  reproduction_df_long <- reproduction_df %>% 
    gather("deletions", 
           "reproduction_rate")
  
  return(reproduction_df_long)
}

# assess the sensitivity of the reproduction number for random attacks
# with number of edge deletions proportional to all edges in network 
reproduction_sensitivity_percentage <- function(igraph_object, 
                                                percent_deletions = c(0.01, 0.1, 0.2, 0.25), 
                                                n_rep = 100, 
                                                numeric_vertices = FALSE, 
                                                county_index = NULL) {
  
  reproduction_list <- data.frame()
  
  for(i in percent_deletions) {
    reproduction_vector <- c()
    
    for (n in seq(1, n_rep)) {
      
      # compute the number of edges to delete corresponding to the deletion 
      # percentage, round up if necessary 
      n_deletions <- (igraph_object %>% V() %>% length() * i
      ) %>% 
        ceiling()
      
      # randomly delete a certain number of edges 
      random_edge_sample <- sample(igraph_object %>% E(),
                                   size = n_deletions)
      
      igraph_attacked <- igraph_object %>% 
        delete_edges(edges = random_edge_sample)
      
      # re-compute reproduction number 
      reproduction_vector <- c(reproduction_vector, 
                               reproduction_rate_lv(igraph_attacked, 
                                                    numeric_vertices = numeric_vertices, 
                                                    county_index = county_index))
    }
    reproduction_add <- list(reproduction_vector)
    names(reproduction_add) <- i 
    reproduction_list <- c(reproduction_list, reproduction_add)
  }
  
  reproduction_df <- do.call(cbind.data.frame, reproduction_list)
  
  reproduction_df_long <- reproduction_df %>% 
    gather("deletions", 
           "reproduction_rate")
  
  return(reproduction_df_long)
}


# ARIMA -------------------------------------------------------------------
# fit ARIMA model, predict the last weeks according to the forecasting 
# window and compute MASE
fit_and_predict_arima <- function(counties = c("Dublin",
                                               "Wicklow", 
                                               "Kerry", 
                                               "Donegal"), 
                                  lag = 2, 
                                  forecast_window = 10) {
  
  arima_mase <- data.frame("time" = as.Date(NA), 
                           "CountyName" = NA, 
                           "true" = NA, 
                           "predicted" = NA, 
                           "res" = NA, 
                           "mase" = NA, 
                           "type" = NA)
  
  for (county in counties) {
    covid_cases_county <- covid_cases_df %>% 
      dplyr::select(county %>% all_of()) %>% 
      ts(frequency = 52, start = c(2020, 1))
    
    # fit ARIMA model 
    arima_mod <- arima(covid_cases_county[-((120 - forecast_window + 1) : 120), ],
                       order = c(lag, 0, 0))
    
    # predict based on ARIMA model 
    prediction <- data.frame(time = rownames(covid_cases_df)[((120 - forecast_window + 1) : 120)])
    prediction$CountyName <- county
    prediction$true <- covid_cases_df %>% 
      dplyr::pull(county %>% all_of()) %>% 
      tail(forecast_window)
    prediction$predicted <- predict(arima_mod, n.ahead = forecast_window)$pred
    prediction$res <- prediction$true - prediction$predicted
    
    prediction$mase <- 0
    prediction$type <- "ARIMA"
    
    
    # compute denominator for MASE 
    denominator <- diff(prediction$true, lag = 1) %>% 
      abs() %>% 
      mean()
    
    # compute MASE 
    for (i in seq(1, nrow(prediction))) {
      prediction[i, ]$mase <- abs(prediction[i, ]$res) / denominator
    }
    arima_mase <- rbind(arima_mase, prediction)
  }
  return(arima_mase)
}


# Economic hub network ----------------------------------------------------
# create economic hub network based on Queen's contiguity network and edges to
# nearest economic hub 
create_eco_hub_igraph <- function(dist_df, # pairwise distance between vertices in data frame
                                  coord) { # matrix of vertex  coordinates
  queen_adj <- covid_net_queen_igraph %>% 
    as_adjacency_matrix()
  
  for (county in queen_adj %>% rownames()) {
    # select nearest economic hub for county in question 
    nearest_hub <- dist_df[county, hubs] %>% which.min() %>% names()
    
    queen_adj[county, nearest_hub] <- 1
    queen_adj[nearest_hub, county] <- 1
  }
  
  index_list <- data.frame("county" = queen_adj %>% rownames(), 
                           "index" = seq(1, 26))
  
  # assign numbers for row names and column names  
  rownames(queen_adj) <- seq(1, 26)
  colnames(queen_adj) <- seq(1, 26)
  
  # create igraph object from adjacency matrix
  covid_net_eco_hubs_igraph <- graph_from_adjacency_matrix(adjmatrix = queen_adj, 
                                                           mode = "undirected") 
  
  # order coordinates in the same order as adjacency matrix
  coord_ord <- coord[match(index_list$county, 
                           rownames(coord)), ] %>% 
    as.data.frame()
  coord_hubs <- coord_ord %>% filter(rownames(coord_ord) %in% hubs)
  
  return(list("igraph_net" = covid_net_eco_hubs_igraph, 
              "ordered_coord" = coord_ord, 
              "ordered_hubs" = coord_hubs))
}


# igraph ------------------------------------------------------------------
# generate nb list from igraph object
# requires numeric vertices (i.e. no county names for vertex names)
igraph2nb <- function(gr) {
  edges <- get.edgelist(gr)
  mode(edges) <- "integer"
  return(neig(edges = edges) %>% neig2nb())
}



# GNAR models --------------------------------------------------------------
# fit GNAR model with alpha and beta order according to input arguments  
fit_and_predict <- function(alpha, beta, 
                            globalalpha, 
                            net,
                            vts = covid_cases, 
                            numeric_vertices = FALSE,
                            
                            # if not NULL, coefficients are computed for 
                            # vertex classes  
                            weight_factor = NULL, 
                            
                            # if TRUE, fits INV-D weighting 
                            inverse_distance = TRUE, 
                            
                            # Great circle distance matrix with county names as row- / colnames 
                            distance_matrix = dist_urbanisation %>% as.matrix(), 
                            
                            # data frame with column CountyName and column weight
                            # for population_weight, computes PB weighting 
                            weight_index = population_weight, 
                            
                            # data frame with column CountyName and its numerical encoding 
                            county_index = NULL, 
                            
                            # if TRUE, the original GNARfit() function is applied
                            old = FALSE, 
                            
                            return_model = FALSE, 
                            forecast_window = 10
                            ) {
  
  train_window <- dim(vts)[1] - forecast_window
  
  # fit model according to given settings 
  if (weight_factor %>% is.null()) {
    if (!old) {
      model <- GNARfit_weighting(vts = vts[1:train_window, ], 
                                 net = net, 
                                 numeric_vertices = numeric_vertices, 
                                 alphaOrder = alpha, 
                                 betaOrder = beta, 
                                 globalalpha = globalalpha, 
                                 inverse_distance = inverse_distance,
                                 distance_matrix = distance_matrix,
                                 weight_index = weight_index, 
                                 county_index = county_index
      )
    } 
    if (old) {
      model <- GNARfit(vts = vts[1:train_window, ], 
                    net = net,
                    alphaOrder = alpha, 
                    betaOrder = beta, 
                    globalalpha = globalalpha
                    )
    }
    
  } else {
    if (!old) {
      model <- GNARfit_weighting(vts = vts[1:train_window, ],
                                 net = net, 
                                 numeric_vertices = numeric_vertices, 
                                 alphaOrder = alpha, 
                                 betaOrder = beta, 
                                 globalalpha = globalalpha, 
                                 fact.var = weight_factor, 
                                 inverse_distance = inverse_distance,
                                 distance_matrix = distance_matrix,
                                 weight_index = weight_index, 
                                 county_index = county_index
      )
    }
    if (old) {
      model <- GNARfit(vts = vts[1:train_window, ], 
                       net = net,
                       alphaOrder = alpha, 
                       betaOrder = beta, 
                       globalalpha = globalalpha, 
                       fact.var = weight_factor
      )
    }
  }
  
  if (!return_model) {
    # return data frame with RSS and BIC value for model 
    return(data.frame("RSS" = model$mod$residuals^2 %>% sum(), 
                      "BIC" = BIC(model)))
  } 
  if (return_model) {
    # return model 
    return(model)
  }
}

# circle through different alpha and beta order combinations and fit 
# GNAR models with the function "fit_and_predict" 
fit_and_predict_for_many <- function(alpha_options = seq(1, 5), 
                                     beta_options = list(0, 1, 
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
                                                         c(2, 2, 2, 1), 
                                                         c(1, 0, 0, 0, 0), 
                                                         c(1, 1, 0, 0, 0), 
                                                         c(1, 1, 1, 0, 0), 
                                                         c(1, 1, 1, 1, 0), 
                                                         c(1, 1, 1, 1, 1),
                                                         c(2, 1, 1, 1, 1), 
                                                         c(2, 2, 1, 1, 1), 
                                                         c(2, 2, 2, 1, 1)
                                                         ), # in form of list
                                     globalalpha = c("TRUE", "FALSE"), 
                                     net, vts = covid_cases,  
                                     numeric_vertices = FALSE, 
                                     
                                     # if not NULL, coefficients are computed for 
                                     # vertex classes
                                     weight_factor = NULL, 
                                     
                                     # if TRUE, fits INV-D weighting 
                                     inverse_distance = TRUE, 
                                     
                                     # Great circle distance matrix with county names as row- / colnames 
                                     distance_matrix = dist_urbanisation %>% as.matrix(),
                                     
                                     # data frame with column CountyName and column weight
                                     # if population_weight, fits PB weighting 
                                     weight_index = population_weight, 
                                     
                                     # data frame with column CountyName and its numerical encoding 
                                     county_index = NULL, 
                                     
                                     # if TRUE, the original GNARfit() function is applied
                                     old = FALSE,
                                     
                                     forecast_window = 10) {
  
  train_window <- dim(vts)[1] - forecast_window
  
  # create all possible combinations of alpha and beta order 
  model_options <-  expand.grid(alpha_options, 
                                beta_options, 
                                globalalpha)
  
  # filter out valid parameter combinations 
  model_options$valid <- (model_options$Var1 == model_options$Var2 %>% 
                            lapply(FUN = function(i) {length(i)}) %>% 
                            unlist())
  model_options_valid <- model_options %>% filter(valid == TRUE)
  
  # for fully connected graph, only the first stage neighbourhood can be 
  # constructed 
  if (any(net %>% 
          GNARtoigraph() %>% 
          degree() == net$edges %>% 
          length() - 1)) {
    only_1 <- model_options_valid %>% 
      pull(Var2) %>% 
      lapply(FUN = function(i) {not(2 %in% i)}) %>% 
      unlist() 
    
    model_options_valid <- model_options_valid[only_1, ]
  }
  
  BIC_RSS <- data.frame()
  
  for (i in seq(1, nrow(model_options_valid))) {
    # select model settings 
    model_setting <- model_options_valid[i, ]
    
    # create name 
    name <- paste("GNAR", 
                  model_setting$Var1 %>% as.character(), 
                  model_setting$Var2[[1]] %>% paste0(collapse = ""), 
                  model_setting$Var3 %>% as.character(), 
                  sep = "-")
    
    # fit model
    results <- fit_and_predict(alpha = model_setting$Var1, 
                               beta = model_setting$Var2[[1]], 
                               globalalpha = model_setting$Var3 %>% as.logical(), 
                               net = net, vts = vts, 
                               numeric_vertices = numeric_vertices, 
                               weight_factor = weight_factor, 
                               inverse_distance = inverse_distance,
                               distance_matrix = distance_matrix,
                               weight_index = weight_index, 
                               county_index = county_index, 
                               old = old, 
                               forecast_window = forecast_window)
    results$name <- name
    BIC_RSS <- rbind(BIC_RSS, 
                     results)
    
  }
  
  return(BIC_RSS[, c(3, 1, 2)])
}


# fit GNAR models according to certain beta- and alpha-order for all data subsets
fit_subset_models <- function(net, 
                              beta_list, 
                              alpha_vector) {
  
  mod_list <- list()
  
  for (i in seq(1, 5)) {
    mod <- fit_and_predict(alpha = alpha_vector[i], 
                           beta = beta_list[[i]], 
                           net = net, 
                           vts = datasets_list[[i]], 
                           globalalpha = TRUE, 
                           old = TRUE,
                           forecast_window = 0, 
                           return_model = TRUE)
    mod_list[[i]] <- mod
  } 
  return(mod_list)
}


# return model with lowest BIC from a generated results list 
return_best_model <- function(results_list) {
  return(results_list[which.min(results_list$BIC), c(1, 3)])
}

# return list of best performing models for KNN / DNN
return_best_knn_dnn <- function(df) {
  
  res_id <- df[which.min(df$best.model.id.BIC), ][, c(1, 2, 3)]
  res_pop <-df[which.min(df$best.model.pop.BIC), ][, c(1, 4, 5)]
  res_id_class <- df[which.min(df$best.model.id.class.BIC), ][, c(1, 6, 7)]
  res_pop_class <- df[which.min(df$best.model.pop.class.BIC), ][, c(1, 8, 9)]
  res_old <- df[which.min(df$best.model.old.BIC), ][, c(1, 10, 11)]
  res_old_class <- df[which.min(df$best.model.old.class.BIC), ][, c(1, 12, 13)]
  
  return(mapply(c, 
                res_id, 
                res_id_class,
                res_pop,  
                res_pop_class, 
                res_old, 
                res_old_class
                ) %>% data.frame())
}


# create plot visualising how BIC changes across weighting schemes
plot_BIC_for_GNAR <- function(results_id, 
                              results_id_class, 
                              results_pop_weighted, 
                              results_pop_weighted_class, 
                              results_old, 
                              results_old_class) {
  results_id$type <- "INV-D"
  results_id_class$type <- "INV-D+class"
  results_pop_weighted$type <- "PB"
  results_pop_weighted_class$type <- "PB+class"
  results_old$type <- "SPL"
  results_old_class$type <- "SPL+class"
  
  results_id$global <- regmatches(results_id$name, 
                                  regexpr('(?:TRUE|FALSE)',
                                          results_id$name))
  results_id_class$global <- regmatches(results_id_class$name, 
                                        regexpr('(?:TRUE|FALSE)',
                                                results_id_class$name))
  results_pop_weighted$global <- regmatches(results_pop_weighted$name, 
                                            regexpr('(?:TRUE|FALSE)',
                                                    results_pop_weighted$name))
  results_pop_weighted_class$global <- regmatches(results_pop_weighted_class$name, 
                                                  regexpr('(?:TRUE|FALSE)',
                                                          results_pop_weighted_class$name))
  results_old$global <- regmatches(results_old$name, 
                                   regexpr('(?:TRUE|FALSE)',
                                           results_old$name))
  results_old_class$global <- regmatches(results_old_class$name, 
                                   regexpr('(?:TRUE|FALSE)',
                                           results_old_class$name))
  
  
  results_id$model <- lapply(regmatches(results_id$name, 
                                        regexpr('-(?:TRUE|FALSE)',
                                                results_id$name), 
                                        invert = TRUE), FUN = function(i) i[[1]]) %>% unlist()
  results_id_class$model <- lapply(regmatches(results_id_class$name, 
                                              regexpr('-(?:TRUE|FALSE)',
                                                      results_id_class$name), 
                                              invert = TRUE), FUN = function(i) i[[1]]) %>% unlist()
  results_pop_weighted$model <- lapply(regmatches(results_pop_weighted$name, 
                                                  regexpr('-(?:TRUE|FALSE)',
                                                          results_pop_weighted$name), 
                                                  invert = TRUE), FUN = function(i) i[[1]]) %>% unlist()
  results_pop_weighted_class$model <- lapply(regmatches(results_pop_weighted_class$name, 
                                                        regexpr('-(?:TRUE|FALSE)',
                                                                results_pop_weighted_class$name), 
                                                        invert = TRUE), FUN = function(i) i[[1]]) %>% unlist()
  results_old$model <- lapply(regmatches(results_old$name, 
                                         regexpr('-(?:TRUE|FALSE)',
                                                 results_old$name), 
                                         invert = TRUE), FUN = function(i) i[[1]]) %>% unlist()
  results_old_class$model <- lapply(regmatches(results_old_class$name, 
                                               regexpr('-(?:TRUE|FALSE)',
                                                       results_old_class$name), 
                                               invert = TRUE), FUN = function(i) i[[1]]) %>% unlist()
  
  
  
  res_complete <- rbind.data.frame(results_id, 
                                   results_id_class, 
                                   results_pop_weighted,
                                   results_pop_weighted_class, 
                                   results_old, 
                                   results_old_class)
  
  global.label <- c("non-global", 
                   "global")
  names(global.label) <- c("FALSE", "TRUE")
  
  
  g <- ggplot(res_complete, 
              aes(x = model, 
                  y = BIC, 
                  color = type, 
                  group = interaction(global, type))) +
    geom_point() +
    geom_line(linetype = "dashed") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1), 
          legend.position = "bottom") + 
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(title = "Weighting scheme")) +
    facet_grid( ~ global, 
                labeller = labeller(global = global.label))
  
  return(g)
}

# compute and plot prediction for the n_ahead last observed weeks for 
# certain counties vs. their real values 
check_and_plot_predictions <- function(model, 
                                       network_name, 
                                       n_ahead = 10, 
                                       county_vector = c("Dublin",  
                                                         "Wicklow", 
                                                         "Kerry", 
                                                         "Donegal")) {
  predicted_df <- predict(model,  
                          n.ahead = n_ahead) %>% 
    as.data.frame()
  colnames(predicted_df) <- colnames(covid_cases)
  
  # start date for period to predict 
  prediction_time <- covid_cases[1:(120 - n_ahead), ] %>% 
    rownames() %>% 
    tail(1)
  
  predicted_df$time <- seq(as.Date(prediction_time) + 7, 
                           as.Date("2022/06/12"), 
                           by = 7)
  
  
  # compare true 1-lag COVID-19 ID and predicted values 
  true <- covid_cases_df[(120 - n_ahead):120, ] %>% 
    gather(key = "CountyName",
           value = "true", 
           -time)
  pred <- predicted_df %>% 
    gather(key = "CountyName",
           value = "predicted", 
           -time)
  
  # create data frame to assess predictive performance 
  check_predictions_df <- left_join(true, pred, 
                                    by = c("CountyName", "time")) %>% 
    mutate(res = true - predicted)
  
  # scatter plot for predictions vs. true values 
  for (county in county_vector) {
    g <- ggplot(data = check_predictions_df %>% 
                  filter(CountyName == county), 
                aes(x = as.Date(time), 
                    y = true,
                    color = "Real", 
                    shape = "Real")) +
      geom_point() +
      geom_point(aes(x = as.Date(time), 
                     y = predicted, 
                     color = "Predicted", 
                     shape = "Predicted")) +
      xlab("Time") +
      ylab("Weekly COVID-19 ID") +
      scale_color_brewer(palette = "Set1") +
      guides(color = guide_legend(title = ""), 
             shape = guide_legend(title = "")) +
      theme(legend.position = "bottom")
    ggsave(filename = paste0("plots/modelfit/Prediction/prediction_", 
                             network_name, "_", county, ".pdf"), 
           plot = g, 
           width = 23, height = 14, unit = "cm")
  }

  return(list("dataset" = check_predictions_df))
}

# compute and plot residuals for GNAR model 
check_and_plot_residuals <- function(model, 
                                     network_name, 
                                     alpha = 4, 
                                     n_ahead = 10, 
                                     county_vector = c("Dublin",
                                                       "Wicklow", 
                                                       "Kerry", 
                                                       "Donegal")) {
  fitted_df <- residuals(model) %>% as.data.frame()
  colnames(fitted_df) <- colnames(covid_cases)
  # assign time column 
  fitted_df$time <- rownames(covid_cases)[-c(1 : alpha, (120 - n_ahead + 1) : 120)]
  
  fitted_df <- fitted_df %>% gather("CountyName", "residuals", -time)
  
  for (county in county_vector) {
    # residual scatter plot for selected counties
    g5 <- ggplot(data = fitted_df %>% 
                   filter(CountyName == county),
                 aes(x = as.Date(time), 
                     y = residuals, 
                     group = CountyName)) +
      geom_point() +
      xlab("Time") +
      ylab("residuals") +
      theme(axis.text.x = element_text(angle = 0, 
                                       vjust = 0.5, 
                                       hjust=1))
    
    ggsave(filename = paste0("plots/modelfit/Residuals/residuals_", 
                             network_name, "_county_",  
                             county, ".pdf"),
           plot = g5, 
           width = 10, height = 10, unit = "cm")
    
    # residual QQ plot for selected counties
    g6 <- ggplot(data = fitted_df %>% 
                   filter(CountyName == county),
                 aes(sample = residuals)) +
      geom_qq() +
      geom_qq_line() +
      xlab("theor. quantiles") +
      ylab("emp. quantiles")
    ggsave(filename = paste0("plots/modelfit/QQ/qq_", 
                             network_name, "_county_",  
                             county, ".pdf"), 
           plot = g6, 
           width = 10, height = 10, unit = "cm")
  }
  return(fitted_df)
}

# Optimise beta -----------------------------------------------------------
# for best performing model, explore influence of a change in beta order at 
# different time lags by fitting GNAR models with different beta orders, 
# selecting the best performing one based on the BIC and plotting the 
# changes in BIC values 

optimise_beta <- function(GNAR_object, 
                          alpha = 4, 
                          beta, 
                          upper_limit_beta, 
                          global_alpha = TRUE, 
                          old = TRUE, 
                          inverse_distance = FALSE, 
                          forecast_window = 10, 
                          
                          # insert beta order for the best performing GNAR model 
                          old_best, 
                          
                          network_name, 
                          county_index = NULL) {
  beta_BIC <- data.frame(beta_vector = NA, 
                         beta = NA, 
                         BIC = NA, 
                         position = NA)
  
  # circle through all time lags and possible neighbourhood stages
  for (position in seq(1, alpha)) {
    for (i in seq(0, upper_limit_beta)) {
      beta_v <- beta
      beta_v[position] <- i
      
      # fit GNAR model 
      model_beta <- fit_and_predict(alpha = alpha,
                                    beta = beta_v, 
                                    net = GNAR_object, 
                                    old = old, 
                                    inverse_distance = inverse_distance, 
                                    globalalpha = global_alpha, 
                                    return_model = TRUE, 
                                    forecast_window = forecast_window, 
                                    county_index = county_index)
      beta_BIC <- rbind(beta_BIC, 
                        data.frame(beta_vector = paste0(beta_v, collapse = ""), 
                                   beta = i, 
                                   BIC = BIC(model_beta), 
                                   position = position))
      
    }
  }
  
  beta_BIC <- beta_BIC %>% na.omit()
  
  # plot change in BIC value for change in beta order across time lags 
  g <- ggplot(beta_BIC, 
              aes(x = beta, 
                  y = BIC, 
                  group = position, 
                  color = position %>% as.character())) +
    geom_point() +
    geom_line(linetype = "dashed") +
    geom_hline(aes(yintercept = old_best), 
               color = "grey", 
               linetype = "dashed") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Time lag")) +
    scale_color_brewer(palette = "Set1") +
    xlab(TeX("$\\beta$-order component"))
  
  ggsave(paste0("plots/modelfit/NeighbourhoodSize/neighbourhood_size_", network_name, ".pdf"), 
         plot = g, 
         width = 23, height = 14, unit = "cm")
  
  return(list("plot" = g, 
              "best_beta" = beta_BIC[which.min(beta_BIC$BIC), ]$beta_vector))
}

# Restrictions - best model -----------------------------------------------
# analyse the change in coefficient values for best performing GNAR models 
# if fitted to different data subsets 
parameter_development <- function(data_list = datasets_list, 
                                  net, 
                                  network_name,
                                  numeric_vertices = FALSE, 
                                  county_index, 
                                  old = TRUE, 
                                  inverse_distance = FALSE, 
                                  forecast_window = 0, 
                                  alpha = 4, 
                                  beta = c(2, 2, 1, 1), 
                                  globalalpha = TRUE,
                                  return_model = TRUE) {
  param <- list()
  
  # circle through all data subsets and fit the best performing GNAR model
  for (i in seq(1, length(data_list))) {
    mod <- fit_and_predict(vts = data_list[[i]], 
                           net = net, 
                           county_index = county_index, 
                           old = old, 
                           inverse_distance = inverse_distance, 
                           forecast_window = forecast_window, 
                           alpha = alpha, 
                           beta = beta, 
                           globalalpha = globalalpha,
                           return_model = return_model, 
                           numeric_vertices = numeric_vertices)
    
    param[[i]] <- mod %>% coef()
  }
  
  param_df <- do.call(rbind.data.frame, param)
  colnames(param_df) <- param[[1]] %>% names()
  param_df$restriction <-  factor(c("Start",
                                    "County-sp. lockdowns", 
                                    "Level-5 lockdown", 
                                    "Indoor dining", 
                                    "End"), 
                                  levels = c("Start", 
                                             "County-sp. lockdowns", 
                                             "Level-5 lockdown", 
                                             "Indoor dining", 
                                             "End"
                                  ))
  
  param_long <- param_df %>% gather("param", 
                                    "values", 
                                    -restriction)
  
  # visualise the change in coefficient values across data subsets  
  g_alpha <- ggplot(param_long  %>% 
                      filter(grepl("alpha", param)), 
                    aes(x = restriction, 
                        y = values, 
                        group = param, 
                        color = param)) +
    geom_point() +
    geom_line(linetype  = "dashed") +
    ylab("coefficient values") +
    xlab("Restrictions") + 
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = paste0("plots/coefficientDevelopment/alpha_order_", network_name, ".pdf"), 
         plot = g_alpha, 
         width = 26, height = 14, unit = "cm")
  
  
  g_beta <- ggplot(param_long  %>% 
                     filter(grepl("beta", param)), 
                   aes(x = restriction, 
                       y = values, 
                       group = param, 
                       color = param)) +
    geom_point() +
    geom_line(linetype = "dashed") +
    ylab("coefficient values") +
    xlab("Restrictions") + 
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "GNAR coefficient"))
  
  ggsave(file = paste0("plots/coefficientDevelopment/beta_order_", network_name, ".pdf"), 
         plot = g_beta, 
         width = 26, height = 14, unit = "cm")
  
  return(list("alpha" = g_alpha, 
              "beta" = g_beta, 
              "dataset" = param_long))
  
}


# Restrictions - tuning parameters ----------------------------------------
# circle through different alpha and beta order combinations and fit 
# GNAR models with the function "fit_and_predict" for each data subset 
fit_and_predict_for_restrictions <- function(alpha_options = seq(1, 5), 
                                             beta_options = list(0, 1, 
                                                                 c(1, 0), 
                                                                 c(1, 1), 
                                                                 c(2, 0), 
                                                                 c(2, 1), 
                                                                 c(2, 2), 
                                                                 c(1, 0, 0), 
                                                                 c(2, 0, 0),
                                                                 c(3, 0, 0),
                                                                 c(4, 0, 0),
                                                                 c(5, 0, 0),
                                                                 c(1, 1, 0),
                                                                 c(2, 1, 0),
                                                                 c(3, 1, 0),
                                                                 c(4, 1, 0),
                                                                 c(5, 1, 0),
                                                                 c(2, 2, 0),
                                                                 c(1, 1, 1), 
                                                                 c(2, 1, 1),
                                                                 c(3, 1, 1),
                                                                 c(4, 1, 1),
                                                                 c(5, 1, 1),
                                                                 c(2, 2, 1), 
                                                                 c(2, 2, 2),
                                                                 c(3, 2, 2), 
                                                                 c(4, 2, 2),
                                                                 c(5, 2, 2), 
                                                                 c(1, 0, 0, 0), 
                                                                 c(2, 0, 0, 0),
                                                                 c(3, 0, 0, 0),
                                                                 c(4, 0, 0, 0),
                                                                 c(5, 0, 0, 0),
                                                                 c(1, 1, 0, 0), 
                                                                 c(2, 1, 0, 0),
                                                                 c(3, 1, 0, 0),
                                                                 c(4, 1, 0, 0),
                                                                 c(5, 1, 0, 0),
                                                                 c(2, 2, 0, 0),
                                                                 c(1, 1, 1, 0), 
                                                                 c(2, 1, 1, 0),
                                                                 c(3, 1, 1, 0), 
                                                                 c(4, 1, 1, 0),
                                                                 c(5, 1, 1, 0),
                                                                 c(2, 2, 1, 0), 
                                                                 c(2, 2, 2, 0),
                                                                 c(3, 2, 2, 0),
                                                                 c(4, 2, 2, 0),
                                                                 c(5, 2, 2, 0),
                                                                 c(1, 1, 1, 1), 
                                                                 c(2, 1, 1, 1),
                                                                 c(3, 1, 1, 1), 
                                                                 c(4, 1, 1, 1),
                                                                 c(5, 1, 1, 1),
                                                                 c(2, 2, 1, 1),
                                                                 c(2, 2, 2, 1), 
                                                                 c(3, 2, 2, 1), 
                                                                 c(4, 2, 2, 1),
                                                                 c(5, 2, 2, 1),
                                                                 c(2, 2, 2, 2),
                                                                 c(1, 0, 0, 0, 0),
                                                                 c(2, 0, 0, 0, 0),
                                                                 c(3, 0, 0, 0, 0), 
                                                                 c(4, 0, 0, 0, 0),
                                                                 c(5, 0, 0, 0, 0), 
                                                                 c(6, 0, 0, 0, 0), 
                                                                 c(7, 0, 0, 0, 0), 
                                                                 c(1, 1, 0, 0, 0),
                                                                 c(2, 1, 0, 0, 0), 
                                                                 c(3, 1, 0, 0, 0),
                                                                 c(4, 1, 0, 0, 0),
                                                                 c(5, 1, 0, 0, 0),
                                                                 c(6, 1, 0, 0, 0),
                                                                 c(7, 1, 0, 0, 0),
                                                                 c(2, 2, 0, 0, 0),
                                                                 c(1, 1, 1, 0, 0), 
                                                                 c(2, 1, 1, 0, 0),
                                                                 c(3, 1, 1, 0, 0), 
                                                                 c(4, 1, 1, 0, 0),
                                                                 c(5, 1, 1, 0, 0),
                                                                 c(6, 1, 1, 0, 0),
                                                                 c(7, 1, 1, 0, 0),
                                                                 c(2, 2, 1, 0, 0), 
                                                                 c(2, 2, 2, 0, 0),
                                                                 c(1, 1, 1, 1, 0), 
                                                                 c(2, 1, 1, 1, 0),
                                                                 c(3, 1, 1, 1, 0), 
                                                                 c(4, 1, 1, 1, 0),
                                                                 c(5, 1, 1, 1, 0),
                                                                 c(6, 1, 1, 1, 0),
                                                                 c(7, 1, 1, 1, 0),
                                                                 c(2, 2, 1, 1, 0), 
                                                                 c(2, 2, 2, 1, 0),
                                                                 c(2, 2, 2, 2, 0),
                                                                 c(3, 2, 2, 2, 0), 
                                                                 c(4, 2, 2, 2, 0),
                                                                 c(5, 2, 2, 2, 0),
                                                                 c(6, 2, 2, 2, 0),
                                                                 c(7, 2, 2, 2, 0),
                                                                 c(1, 1, 1, 1, 1),
                                                                 c(2, 1, 1, 1, 1),
                                                                 c(3, 1, 1, 1, 1), 
                                                                 c(4, 1, 1, 1, 1),
                                                                 c(5, 1, 1, 1, 1),
                                                                 c(6, 1, 1, 1, 1),
                                                                 c(7, 1, 1, 1, 1),
                                                                 c(2, 2, 1, 1, 1), 
                                                                 c(2, 2, 2, 1, 1),
                                                                 c(3, 2, 2, 1, 1),
                                                                 c(4, 2, 2, 1, 1),
                                                                 c(5, 2, 2, 1, 1),
                                                                 c(6, 2, 2, 1, 1),
                                                                 c(7, 2, 2, 1, 1),
                                                                 c(2, 2, 2, 2, 1), 
                                                                 c(2, 2, 2, 2, 2)
                                             ), # in form of list
                                             globalalpha = TRUE, 
                                             net, 
                                             data_list = datasets_list,  
                                             numeric_vertices = FALSE,
                                             
                                             # if not NULL, fit coefficients for 
                                             # vertex classes 
                                             weight_factor = NULL, 
                                             
                                             # if TRUE, INV-D weighting 
                                             inverse_distance = FALSE, 
                                             
                                             # Great circle distance matrix with county names as row- / colnames 
                                             distance_matrix = dist_urbanisation %>% as.matrix(), 
                                             
                                             # data frame with column CountyName and column weight
                                             weight_index = population_weight, 
                                             
                                             # data frame with column CountyName and its numerical encoding 
                                             county_index = NULL, 
                                             
                                             # if TRUE, the original GNARfit() function is applied
                                             old = TRUE,
                                             
                                             forecast_window = 5,
                                             upper_limit = 5) {
  
  # construct valid beta options list 
  exclude <- lapply(beta_options, 
                    FUN = function(i) any(i > upper_limit)) %>% 
    unlist()
  
  beta_valid_options <- beta_options[!exclude]
  
  param <- list()
  
  for (i in seq(1, length(data_list))) {
    res <- fit_and_predict_for_many(vts = data_list[[i]], 
                                    net = net, 
                                    alpha_options = alpha_options, 
                                    beta_options = beta_valid_options, 
                                    globalalpha = globalalpha,
                                    old = old,
                                    forecast_window = forecast_window, 
                                    numeric_vertices = numeric_vertices, 
                                    weight_factor = weight_factor, 
                                    inverse_distance = inverse_distance, 
                                    distance_matrix = distance_matrix, 
                                    weight_index = weight_index,
                                    county_index = county_index)
    
    param[[i]] <- return_best_model(res)
  }
  param_df <- do.call(rbind.data.frame, param)
  
  param_df$data_subset <- seq(1, length(data_list))
  
  return(param_df[, c(3, 1, 2)])
}


# MASE --------------------------------------------------------------------
# compute MASE for certain counties 
compute_MASE <- function(model, 
                         network_name, 
                         n_ahead = 10, 
                         counties = c("Dublin", 
                                      "Wicklow", 
                                      "Kerry", 
                                      "Donegal"), 
                         data_df = covid_cases_df) {
  predicted_df <- predict(model,  
                          n.ahead = n_ahead) %>% 
    as.data.frame()
  
  colnames(predicted_df) <- colnames(data_df %>% dplyr::select(-time))
  
  length_data_df <- nrow(data_df)
  
  # start date for the period of prediction 
  prediction_time <- data_df[1:(length_data_df - n_ahead), ] %>% 
    rownames() %>% 
    tail(1)
  
  end_date <- data_df[length_data_df, ] %>% 
    rownames()
  
  predicted_df$time <- seq(as.Date(prediction_time) + 7, 
                           as.Date(end_date), 
                           by = 7)
  
  # compare true 1-lag COVID-19 ID and predicted values 
  true <- data_df[(length_data_df - n_ahead + 1):length_data_df, ] %>% 
    gather(key = "CountyName",
           value = "true", 
           -time)
  pred <- predicted_df %>% 
    gather(key = "CountyName",
           value = "predicted", 
           -time)
  # create data frame to assess predictive performance 
  check_predictions_df <- left_join(true, pred, 
                                    by = c("CountyName", "time")) %>% 
    mutate(res = true - predicted)
  
  mase_df <- data.frame("time" = as.Date(NA), 
                        "CountyName" = NA, 
                        "true" = NA, 
                        "predicted" = NA, 
                        "res" = NA, 
                        "mase" = NA, 
                        "type" = NA)
  
  for (county in counties) {
    check_predictions_county <- check_predictions_df %>% 
      filter(CountyName == county) 
    
    check_predictions_county$mase <- 0
    
    # compute denominator for MASE
    denominator <- diff(check_predictions_county$true, lag = 1) %>% 
      abs() %>% 
      mean()
    
    for (i in seq(1, nrow(check_predictions_county))) {
      # compute MASE values 
      check_predictions_county[i, ]$mase <- abs(check_predictions_county[i, ]$res) / denominator
      }
    
    check_predictions_county$type <- network_name
    
    mase_df <- rbind.data.frame(mase_df, 
                                check_predictions_county)
    
  }
  
  return(mase_df %>% na.omit())
}


# compare MASE for two GNAR models in plot 
compare_MASE <- function(mase_old,
                         mase_new, 
                         network_name, 
                         old_model, 
                         new_model
) {
  mase_old$model_version <- "old"
  mase_new$model_version <- "new"
  
  mase_comparison <- rbind(mase_old, 
                           mase_new)
  
  ggplot(mase_comparison, 
         aes(x = time, 
             y = mase,
             group = interaction(CountyName, model_version), 
             color = model_version)) +
    geom_point() +
    geom_line(linetype = "dashed") +
    xlab("Time") + 
    ylab("MASE") +
    facet_grid(~ CountyName) +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1)) +
    scale_color_manual(values = c("old" = "#F8766D", 
                                  "new" = "#A3A500"), 
                       labels = c("old" = paste0("old: ", old_model), 
                                  "new" = paste0("new: ", new_model)), 
                       name = "Model version")
  
  ggsave(paste0("plots/Prediction/Comparison/old_vs_new_", network_name, ".pdf"), 
         width = 23, height = 13, unit = "cm")
  
}


# Simulation data  --------------------------------------------------------
# compute average variance of residuals (for white noise error term)
compute_variance <- function(GNARmodel) {
  res <- GNARmodel %>% 
    residuals() %>% 
    apply(MARGIN = 2, FUN = function(i) var(i)) %>% 
    mean()
  
  return(res)
}

# compute average variance of residuals for each data subset (for white noise error term)
compute_variance_subsets <- function(GNARmodel_list) {
  
  res_vector <- c()
  
  for (i in seq(1, 5)) {
    res <- GNARmodel_list[[i]] %>% 
      residuals() %>% 
      apply(MARGIN = 2, FUN = function(i) var(i)) %>% 
      mean()
    
    res_vector <- c(res_vector, res)
  }
  return(res_vector)
}



# generate simulated data
simulate_time_series <- function(gnar_object, 
                                 GNARmodel, # GNARfit model
                                 beta_order, # best performing beta order 
                                 av_var, # average variance in COVID-19 ID
                                 county_index, 
                                 counties = counties_v, 
                                 timeframe = 20) {
  # determine max stage across lags 
  max_stage <- beta_order %>% max()
  
  # number of lags 
  n_alpha <- beta_order %>% length()
  
  # generate white noise 1-lag COVID-19 ID values for first n days
  initial_df <- rnorm(n = n_alpha  * 26, 
                      mean = 0, 
                      sd = sqrt(av_var)) %>%
    matrix(ncol = 26, 
           nrow = n_alpha) %>% 
    as.data.frame()
  
  # assign counties to columns  
  colnames(initial_df) <- counties
  
  # initialize data frame to store simulated data
  simulated_df <- initial_df
  
  # current time 
  current_end <- initial_df %>% nrow()
  
  # predict for weeks according to time frame
  while (current_end < timeframe) {
    simulated_values <- c()
    
    for (county in counties) {
      
      # create data frame with only relevant county as column
      initial_df_county <- simulated_df %>% dplyr::select(all_of(county))
      
      # numeric index for vertex according to county index
      county_to_vertex <- county_index[county_index$CountyName == county, ]$index
      
      # obtain model coefficients 
      model_coef <- GNARmodel %>% 
        coef() %>%
        as.data.frame()
      colnames(model_coef) <- "param"
      
      model_coef$type <- rownames(model_coef)
      
      # separate alpha and beta coefficients 
      alpha_coef <- model_coef  %>% 
        filter(grepl("alpha", type))
      
      beta_coef <- model_coef  %>% 
        filter(grepl("beta", type))
      
      # compute shortest path length between county and all remaining counties
      spl_county <- distances(graph = gnar_object %>% 
                                GNARtoigraph(), 
                              v = county_to_vertex) %>% 
        t() %>% 
        as.data.frame()
      colnames(spl_county) <- "spl"
      spl_county$index <- seq(1, 26)
      
      # add shortest path length to county index data frame 
      spl_county_names <- left_join(county_index, 
                                    spl_county, 
                                    by = "index")
      
      # create vector for beta terms to be saved in 
      beta_part_v <- c()
      
      for (lagged in seq(1, n_alpha)) {
        # determine neighborhood stage for certain lag
        beta_order_component <- beta_order[lagged]
        
        # determine beta coefficients for certain lag 
        beta_coef_lagged <- beta_coef %>%
          filter(grepl(paste0("beta", lagged), type))
        
        # for every stage, identify neighbors and compute term based on beta coefficients and lagged values 
        if (beta_order_component != 0) {
          for (stage in seq(1, beta_order_component)) {
            stage_neighborhood <- spl_county_names %>% 
              filter(spl == stage) %>% 
              pull(CountyName)
            
            initial_df_neighbors <- simulated_df %>% dplyr::select(all_of(stage_neighborhood))
            
            beta_part <- sum(beta_coef_lagged$param[stage] * (1 / stage) * initial_df_neighbors[current_end + 1 - lagged, ])
            
            beta_part_v <- c(beta_part_v, beta_part)
          }
        }
        if (beta_order_component == 0) {
          beta_part_v <- c(beta_part_v, 0)
        }
        
      }
      
      # compute alpha term 
      alpha_term <- (alpha_coef$param * initial_df_county[(current_end - n_alpha) : (current_end), 1]) %>% 
        sum()
      
      # compute beta term 
      beta_term <- beta_part_v %>% sum()
      
      # compute simulated value plus random error
      simulated_value <- alpha_term + beta_term + rnorm(n = 1, mean = 1, sd = sqrt(av_var))
      
      # add simulated data point to existing vector of simulated data points  
      simulated_values <- c(simulated_values, simulated_value)
      
    }
    # add new row with simulated data points 
    simulated_df <- rbind(simulated_df, simulated_values)
    
    # compute current time end point 
    current_end <- simulated_df %>% nrow()
  }
  
  # assign fictional "time" 
  simulated_df$time <- seq(1, timeframe)
  
  # transform into long data table 
  simulated_df_long <- simulated_df %>% 
    gather("CountyName", 
           "ID", 
           -time)
  
  return(simulated_df_long)
}


# generate simulated data with data subset specific variance 
simulate_time_series_subsets <- function(gnar_object, 
                                 GNARmodel_list, # GNARfit model list for each subset
                                 beta_list, # best performing beta order for each subset
                                 var_vector, # variance vector for subsets
                                 county_index, 
                                 counties = counties_v, 
                                 data_list = datasets_list) {
  
  n_beginning <- beta_list[[1]] %>% length()
  
  # generate white noise 1-lag COVID-19 ID values for first n days
  initial_df <- rnorm(n =  n_beginning * 26, 
                      mean = 0, 
                      sd = sqrt(var_vector[[1]])) %>%
    matrix(ncol = 26, 
           nrow = n_beginning) %>% 
    as.data.frame()
  
  # assign counties to columns  
  colnames(initial_df) <- counties
  
  # initialize data frame to store simulated data
  simulated_df <- initial_df
  simulated_df$subset <- NA
  
  # current time 
  current_end <- initial_df %>% nrow()
  
  
  for (i in seq(1, 5)) {
    # determine max stage across lags 
    max_stage <- beta_list[[i]] %>% max()
    
    # number of lags 
    n_alpha <- beta_list[[i]] %>% length()
    
    # set prediction time 
    timeframe <- datasets_list[[i]] %>% nrow()
    
    # index to control number of time steps to simulate 
    simulated_time_steps <- 1
    
    # predict for weeks according to time frame
    while (simulated_time_steps < timeframe) {
      simulated_time_steps <- simulated_time_steps + 1
      
      simulated_values <- c()
      
      for (county in counties) {
        
        # create data frame with only relevant county as column
        initial_df_county <- simulated_df %>% dplyr::select(all_of(county))
        
        # numeric index for vertex according to county index
        county_to_vertex <- county_index[county_index$CountyName == county, ]$index
        
        # obtain model coefficients 
        model_coef <- GNARmodel_list[[i]] %>% 
          coef() %>%
          as.data.frame()
        colnames(model_coef) <- "param"
        
        model_coef$type <- rownames(model_coef)
        
        # separate alpha and beta coefficients 
        alpha_coef <- model_coef  %>% 
          filter(grepl("alpha", type))
        
        beta_coef <- model_coef  %>% 
          filter(grepl("beta", type))
        
        # compute shortest path length between county and all remaining counties
        spl_county <- distances(graph = gnar_object %>% 
                                  GNARtoigraph(), 
                                v = county_to_vertex) %>% 
          t() %>% 
          as.data.frame()
        colnames(spl_county) <- "spl"
        spl_county$index <- seq(1, 26)
        
        # add shortest path length to county index data frame 
        spl_county_names <- left_join(county_index, 
                                      spl_county, 
                                      by = "index")
        
        # create vector for beta terms to be saved in 
        beta_part_v <- c()
        
        for (lagged in seq(1, n_alpha)) {
          # determine neighborhood stage for certain lag
          beta_order_component <- beta_list[[i]][lagged]
          
          # determine beta coefficients for certain lag 
          beta_coef_lagged <- beta_coef %>%
            filter(grepl(paste0("beta", lagged), type))
          
          # for every stage, identify neighbors and compute term based on beta coefficients and lagged values 
          if (beta_order_component != 0) {
            for (stage in seq(1, beta_order_component)) {
              stage_neighborhood <- spl_county_names %>% 
                filter(spl == stage) %>% 
                pull(CountyName)
              
              initial_df_neighbors <- simulated_df %>% dplyr::select(all_of(stage_neighborhood))
              
              beta_part <- sum(beta_coef_lagged$param[stage] * (1 / stage) * initial_df_neighbors[current_end + 1 - lagged, ])
              
              beta_part_v <- c(beta_part_v, beta_part)
            }
          }
          if (beta_order_component == 0) {
            beta_part_v <- c(beta_part_v, 0)
          }
          
        }
        
        # compute alpha term 
        alpha_term <- (alpha_coef$param * initial_df_county[(current_end - n_alpha) : (current_end), 1]) %>% 
          sum()
        
        # compute beta term 
        beta_term <- beta_part_v %>% sum()
        
        # compute simulated value plus random error
        simulated_value <- alpha_term + beta_term + rnorm(n = 1, mean = 1, sd = sqrt(var_vector[[i]]))
        
        # add simulated data point to existing vector of simulated data points  
        simulated_values <- c(simulated_values, simulated_value)
        
      }
      # assign data subset allocation 
      simulated_values <- c(simulated_values, i)
      
      # add new row with simulated data points 
      simulated_df <- rbind(simulated_df, simulated_values)
      
      # compute current time end point 
      current_end <- simulated_df %>% nrow()
    } # end while for time
  }
  
  # assign fictional "time" 
  simulated_df$time <- seq(1, simulated_df %>% nrow())
  
  # transform into long data table 
  simulated_df_long <- simulated_df %>% 
    gather("CountyName", 
           "ID", 
           -c("time", "subset"))
  return(simulated_df_long)
}


