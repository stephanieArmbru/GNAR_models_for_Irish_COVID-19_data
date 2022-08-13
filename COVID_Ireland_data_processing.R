### PRE-PROCESS COVID-19 DATA FOR IRELAND 

# Install packages --------------------------------------------------------
install.packages("ISOweek")

# set seed to guarantee reproduciblity 
set.seed(1234)

# Load libraries ----------------------------------------------------------
library(readr)
library(magrittr)
library(tidyverse)
library(ISOweek)

library(xtable) # for tables

library(lubridate) # for dates
library(zoo)

library(tseries) # for augmented Dickey-Fuller test
library(latex2exp) # for latex in axis labels 


# Load Data ---------------------------------------------------------------
COVID_data_orig <- read_csv("data/COVID-19_HPSC_full_week.csv", 
                            show_col_types = FALSE)


# retain the original dataset
# limit to data until 18.06.2022, full epidemic week 
COVID_data <- COVID_data_orig %>% as.data.frame() %>% 
  filter(as.Date(TimeStamp) <= as.Date("2022-06-18"))

# Data exploration --------------------------------------------------------
COVID_data %>% nrow() # 21866 rows
COVID_data %>% ncol() # 16 columns 
COVID_data %>% colnames()

# ObjectID: unique id for each row 
# OrigID: id for county, 26 in total
# CountyName 
# PopulationCensus16: 
# TimeStamp: date and time (separated in preprocessing)
# IGEasting / IGNorthing: geographic identifier
# Lat / Long: gives the "middle" of the county (not its county town)
# UGI: url link
# ConfirmedCovidCases
# PopulationProportionCovidCases
# ConfirmedCovidDeaths
# ConfirmedCovidRecovered
# SHAPE_Length/Area

# check if ObjectID and OrigID identical
COVID_data %>% filter(COVID_data$OBJECTID != COVID_data$ORIGID)
COVID_data %>% dplyr::select(OBJECTID) %>% unique() %>% nrow()
COVID_data %>% dplyr::select(ORIGID) %>% unique() %>% nrow()

# compute relative size of counties  
counties_size <- COVID_data %>% 
  dplyr::select(CountyName, 
                PopulationCensus16) %>% 
  unique() %>% 
  mutate(PopulationChar = as.character(PopulationCensus16), 
         PopulationRel = round(100 * PopulationCensus16 / sum(PopulationCensus16), 
                               digits = 2)) %>% 
  dplyr::select(CountyName, 
                PopulationChar, 
                PopulationRel) %>% 
  arrange(desc(PopulationRel))

# for latex 
strCaption <- paste0("Republic of Ireland: counties, size in number of inhabitants 
                     according to the census 2016, relative size in relative number 
                     of inhabitants compared to total population size")

print(xtable(counties_size,
             digits=2,
             caption=strCaption,
             label="tab:ireland_counties", 
             align = c("", "r", "|", "l", "l")),
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      caption.placement="bottom",
      hline.after=NULL,
      add.to.row = list(pos = list(-1,
                                   nrow(counties_size)),
                        command = c(paste("\\toprule \n",
                                          " County & size & rel. size (in \\%)
                                          \\\\\n",
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)




# extract time period for which data was recorded  
COVID_data %>% dplyr::select(TimeStamp) %>% head(1) # from 27.02.2020
COVID_data %>% dplyr::select(TimeStamp) %>% tail(1) # to 18.06.2022


# check for missing dates 
all_dates <- seq(as.Date("2020/02/27 00:00:00+00"), 
                 as.Date("2022/06/18 00:00:00+00"),
                 by = "day") 
missing_dates <- all_dates %in% (COVID_data %>% 
                                   pull(TimeStamp) %>% 
                                   unique() %>% 
                                   as.Date())

all_dates[!missing_dates] # only date "2020-02-28", "2020-02-29" missing 

# COVID-19 -------------------------------------------------------------------

# Confirmed COVID-19 Deaths / Recovered
COVID_data %>% 
  dplyr::select(ConfirmedCovidCases) %>% 
  is.na() %>% 
  sum()

COVID_data %>% 
  dplyr::select(ConfirmedCovidDeaths) %>% 
  is.na() %>% 
  sum()
# ConfirmedCovidDeaths is empty column 

COVID_data %>% 
  dplyr::select(ConfirmedCovidRecovered) %>% 
  is.na() %>% 
  sum() 
# ConfirmedCovidRecovered is empty column


COVID_data %>% 
  dplyr::select(PopulationProportionCovidCases) %>% 
  is.na() %>% 
  sum()
# 52 missing values 

COVID_data %>% 
  filter(is.na(PopulationProportionCovidCases)) %>% 
  dplyr::select(TimeStamp) %>% 
  table()
# NA for beginning of pandemic, i.e. between 27.02.2020 and 01.03.2020
# for every county, the data for 27.02 and 01.03 is recorded but NA

COVID_data %>% 
  filter(TimeStamp == "2020/02/28 00:00:00+00")
COVID_data %>% 
  filter(TimeStamp == "2020/02/29 00:00:00+00")
# data for 28.-29.02.2020 is missing  

COVID_data %>% 
  dplyr::select(TimeStamp) %>% 
  unique() %>% 
  head(1)
# first recording on the 27.02.2020

COVID_data %>% 
  filter(PopulationProportionCovidCases != 0, 
         TimeStamp == "2020/03/02 00:00:00+00") %>% 
  dplyr::select(CountyName)
# on 02.03.2020 the first COVID-19 case was recorded in Dublin  


COVID_data %>% 
  group_by(TimeStamp) %>% 
  summarise(count_counties = sum(PopulationProportionCovidCases != 0)) %>% 
  filter(count_counties == 26) %>% 
  head(1)
# on 20.03.2020 COVID-19 was recorded in every county


# Data Pre-processing -----------------------------------------------------
# make object ID and origine ID discrete and 
# remove NA-filled Deaths and Recoveries columns 
COVID_data <- COVID_data %>% 
  mutate(id = as.factor(OBJECTID), 
         origin = as.factor(ORIGID)) %>% 
  dplyr::select(-c(ConfirmedCovidDeaths, ConfirmedCovidRecovered))

# separate time into date and time
COVID_data <- COVID_data %>% 
  separate(col = "TimeStamp", 
           into = c("date", "time"), 
           sep = " ") %>% 
  mutate("date_f" = as.Date(date))

# rank for size 
size_rank <- COVID_data %>%  
  dplyr::select(CountyName, PopulationCensus16) %>% 
  unique() %>%  
  mutate(SizeRank = rank(PopulationCensus16)) %>% 
  dplyr::select(-PopulationCensus16)

COVID_data <- left_join(COVID_data, size_rank, by = "CountyName")

# remove missing dates
COVID_data <- COVID_data %>% 
  drop_na()


# compute daily COVID-19 incidence for 100,000 inhabitants
COVID_data <- COVID_data %>% 
  group_by(CountyName) %>% 
  mutate(dailyCovid = c(0, 
                        diff(PopulationProportionCovidCases, 
                             lag = 1))) %>% 
  ungroup()



# Weekly aggregation ------------------------------------------------------
# aggregate the daily COVID-19 incidence to a weekly COVID-19 incidence 
# for 100,000 inhabitants 
COVID_data_week_agg <- COVID_data %>% 
  dplyr::select(CountyName, 
         date_f, 
         PopulationProportionCovidCases) %>% 
  group_by(CountyName, 
           yw = floor_date(date_f, unit = "week")) %>% 
  summarise(weeklyCasesSum = mean(PopulationProportionCovidCases)) %>% 
  mutate(weeklyCases = c(0, 
                         diff(weeklyCasesSum, 
                              lag = 1))) %>% 
  ungroup()


census <- COVID_data %>% 
  dplyr::select(CountyName, 
         PopulationCensus16) %>% 
  unique()

COVID_data_week_agg <- left_join(x = COVID_data_week_agg, 
                                 y = census, 
                                 by = "CountyName")

COVID_data_week_agg$yw %>% unique %>% length()
# 120 weekly COVID-19 observations for each county 

# plot weekly COVID-19 incidence for each county 
ggplot(COVID_data_week_agg, 
       aes(x = yw, 
           y = weeklyCases, 
           color = CountyName)) + 
  geom_line() +
  xlab("Time") +
  ylab("COVID-19 cases") +
  guides(color = guide_legend(title = "County"))
ggsave("plots/dataVisualisation/covid_ireland_weekly_cases.pdf", 
       width = 23, height = 14, unit = "cm")


# Exploration -------------------------------------------------------------
# highest incidence (cumulative)
COVID_data_week_agg[which.max(COVID_data_week_agg$weeklyCasesSum), ] 
# Louth highest cumulative incidence

# highest incidence (daily)
COVID_data_week_agg[which.max(COVID_data_week_agg$weeklyCases), ] 
# Westmeath in the week of 2022-01-23 with week average 3446


# Moving average winter 21/22 ----------------------------------------------
# isolate peak in weekly COVID-19 incidence during winter 21/22 and 
# smooth by computing the moving average incidence over 4 consecutive 
# weeks 
COVID_data_week_smooth <- COVID_data_week_agg %>% 
  mutate(weeklyCases = ifelse(yw %in% seq(as.Date("2021-12-12"),
                                          as.Date("2022-02-27"), 
                                          by = 7), # 9 weeks over Christmas
                              rollmean(weeklyCases,
                                        k = 4,
                                        align = 'right',
                                        fill = 0), 
                              weeklyCases))

# plot weekly COVID-19 incidence with smoothed winter season 
ggplot(COVID_data_week_smooth, 
       aes(x = yw, 
           y = weeklyCases, 
           color = CountyName)) + 
  geom_line() +
  xlab("Time") +
  ylab("COVID-19-19 cases") +
  guides(color = guide_legend(title = "County"))
ggsave("plots/dataVisualisation/covid_ireland_weekly_cases_winter_smoothed.pdf", 
       width = 23, height = 14, unit = "cm")



# Classical decomposition -------------------------------------------------
# R function decompose()
# decomposition for county Dublin 
COVID_dublin_ts_decomposed <- tsr_decomposition(county_name = "Dublin", 
                                                type_decomposition = "additive")

COVID_dublin_ts_decomposed_mult <- tsr_decomposition("Dublin", 
                                                     type_decomposition = "multiplicative")

pdf("plots/stationarity/decomposition_dublin.pdf")
plot(COVID_dublin_ts_decomposed)
dev.off()

pdf("plots/stationarity/mult_decomposition_dublin.pdf")
plot(COVID_dublin_ts_decomposed_mult)
dev.off()

# Other counties 
COVID_cork_ts_decomposed <- tsr_decomposition("Leitrim")
plot(COVID_cork_ts_decomposed)

COVID_cork_ts_decomposed <- tsr_decomposition("Cork")
plot(COVID_cork_ts_decomposed)

COVID_galway_ts_decomposed <- tsr_decomposition("Galway")
plot(COVID_galway_ts_decomposed)

COVID_carlow_ts_decomposed <- tsr_decomposition("Carlow")
plot(COVID_carlow_ts_decomposed)

COVID_longford_ts_decomposed <- tsr_decomposition("Longford")
plot(COVID_longford_ts_decomposed)


# Manual trend estimation through Moving Average smoothing
COVID_data_week_smooth <- COVID_data_week_smooth %>% 
  group_by(CountyName) %>% 
  mutate(trendCases = rollapply(weeklyCases,
                                24,
                                mean,
                                align = 'right',
                                fill = 0)) %>% 
  ungroup()

plot_trend("Dublin")


## Manual season estimation through Averaging over quarters
COVID_data_week_smooth <- COVID_data_week_smooth %>%
  mutate(detrendCases = weeklyCases - trendCases,
         quarter = lubridate::quarter(yw)) %>%
  group_by(quarter) %>%
  mutate(seasonCases = mean(detrendCases)) %>%
  ungroup()

plot_season("Dublin")


# compute Residuals = True value - trend - season
COVID_data_week_smooth <- COVID_data_week_smooth %>% 
  mutate(randomError = weeklyCases - (trendCases + seasonCases))

plot_residuals("Dublin")


# plot entire decomposition 
plot_tsr("Dublin")



# Transformation ----------------------------------------------------------
# plot log transformed weekly COVID-19 incidence  
ggplot(COVID_data_week_smooth, 
       aes(x = yw, 
           y = log(weeklyCases + 1),
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("log. COVID-19-19 cases") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/stationarity/log_weekly_cases.pdf", 
       width = 27, height = 18, unit = "cm")

# classical decomposition for log(cases)
COVID_data_week_smooth <- COVID_data_week_smooth  %>% 
  mutate(logWeeklyCases = log(weeklyCases + 1))


COVID_dublin_log_ts_decomp <- tsr_decomposition(column_name = "logWeeklyCases", 
                                                county_name = "Dublin", 
                                                type_decomposition = "multiplicative")

pdf("plots/stationarity/decomposition_log_dublin.pdf")
plot(COVID_dublin_log_ts_decomp)
dev.off()

# plot square root transformed weekly COVID-19 incidence
ggplot(COVID_data_week_smooth, 
       aes(x = yw, 
           y = sqrt(weeklyCases),
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("squ.r. COVID-19-19 cases") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
ggsave("plots/stationarity/square_root_weekly_cases.pdf", 
       width = 27, height = 18, unit = "cm")

# classical decomposition for sqrt(cases)
COVID_data_week_smooth <- COVID_data_week_smooth %>% 
  mutate(sqrtWeeklyCases = sqrt(weeklyCases))

COVID_dublin_sqrt_ts_decomp <- tsr_decomposition(column_name = "sqrtWeeklyCases", 
                                                 county_name = "Dublin", 
                                                 type_decomposition = "multiplicative")
pdf("plots/stationarity/decomposition_sqrt_dublin.pdf")
plot(COVID_dublin_sqrt_ts_decomp)
dev.off()

# compute integrated weekly COVID-19 incidence through taking differences
COVID_data_week_smooth <- COVID_data_week_smooth %>%
  group_by(CountyName) %>% 
  mutate(weeklyCasesDiff = c(0, 
                             diff(weeklyCases, 
                                  lag = 1)
  )) %>% 
  ungroup()

# identify smallest difference in COVID-19 incidence 
COVID_data_week_smooth[which.min(COVID_data_week_smooth$weeklyCasesDiff), ]$weeklyCasesDiff


# plot difference in weekly COVID-19 incidence 
# stationary and hence used to fit GNAR models 
ggplot(COVID_data_week_smooth, 
       aes(x = yw, 
           y = weeklyCasesDiff,
           color = CountyName)) +
  geom_line() +
  xlab("Time") +
  ylab("1-lag diff. weekly cases") +
  guides(color = guide_legend(title = "County")) +
  theme(legend.position = "bottom")
# almost heteroscedastic, close enough
ggsave("plots/dataVisualisation/weekly_cases_diff_smoothed.pdf", 
       width = 27, height = 18, unit = "cm")


# Box-Cox transformation --------------------------------------------------
# compute Box-Cox transformation for weekly COVID-19 incidence
pdf("plots/stationarity/boxcox_weekly_cases.pdf")
boxcox(lm(weeklyCases + 1 ~ 1 + CountyName, 
          data = COVID_data_week_smooth))
# 0 : log transformation 
dev.off()

# compute Box-Cox transformation for log transformed COVID-19 incidence 
boxcox(lm(logWeeklyCases + 1 ~ 1 + CountyName, 
          data = COVID_data_week_smooth))
# 1 : no transformation

# compute Box-Cox transformation for 1-lag COVID-19 ID
pdf("plots/stationarity/boxcox_weekly_cases_diff.pdf")
boxcox(lm(weeklyCasesDiff + 693 ~ 1 + CountyName, 
          data = COVID_data_week_smooth))
# 1: no transformation 
dev.off()


# Temporal dependency ----------------------------------------------------
# compute ACF and PACF plots for 1-lag COVID-19 ID for county Dublin 
COVID_county <- COVID_data_week_smooth %>% 
  filter(CountyName == "Dublin")

acf_dublin <- acf(COVID_county$weeklyCasesDiff,
                  plot = FALSE)
pdf("plots/stationarity/acf_dublin.pdf")
plot(acf_dublin, main = "")
dev.off()

pacf_dublin <- pacf(COVID_county$weeklyCasesDiff, 
                    plot = FALSE)
pdf("plots/stationarity/pacf_dublin.pdf")
plot(pacf_dublin, main = "")
dev.off()


# Augmented Dickey-Fuller test -------------------------------------------
# Dublin
COVID_county <- COVID_data_week_smooth %>% 
  filter(CountyName == "Dublin")
adf.test(COVID_county$weeklyCasesDiff, k = 0)
# p-value significant 


# Other counties 
COVID_county <- COVID_data_week_smooth %>% 
  filter(CountyName == "Leitrim")
adf.test(COVID_county$weeklyCasesDiff, k = 0) 
# p-value significant 
COVID_county <- COVID_data_week_smooth %>% 
  filter(CountyName == "Cork")
adf.test(COVID_county$weeklyCasesDiff, k = 0)
# p-value significant 
COVID_county <- COVID_data_week_smooth %>% 
  filter(CountyName == "Mayo")
adf.test(COVID_county$weeklyCasesDiff, k = 0)
# p-value significant 
COVID_county <- COVID_data_week_smooth %>% 
  filter(CountyName == "Donegal")
adf.test(COVID_county$weeklyCasesDiff, k = 0)
# p-value significant 

# Heatmap -----------------------------------------------------------------
# plot heatmaps for 1-lag COVID-19 ID for years 2020, 2021 and 2022
COVID_heatmap <- COVID_final %>% 
  mutate(year = year(yw), 
         week = epiweek(yw))

# year 2020
ggplot(COVID_heatmap %>% filter(year == "2020"), 
       aes(x = week, 
           y = CountyName, 
           fill = weeklyCases)) +
  geom_tile() +
  scale_fill_gradient(low = "red", 
                      high = "green", 
                      limits = c(-692, 944)) +
  xlab("Time (in epid. calendar week)") +
  ylab("County") +
  labs(fill = TeX("$\\Delta$ weekly cases"))
ggsave("plots/dataVisualisation/heatmap_2020.pdf", 
       width = 23, height = 15, unit = "cm")


# year 2021
ggplot(COVID_heatmap %>% filter(year == "2021"), 
       aes(x = week, 
           y = CountyName, 
           fill = weeklyCases)) +
  geom_tile() +
  scale_fill_gradient(low = "red", 
                      high = "green", 
                      limits = c(-692, 944)) +
  xlab("Time (in epid. calendar week)") +
  ylab("County") +
  labs(fill = TeX("$\\Delta$ weekly cases"))
ggsave("plots/dataVisualisation/heatmap_2021.pdf", 
       width = 23, height = 15, unit = "cm")


# year 2022
ggplot(COVID_heatmap %>% filter(year == "2022"), 
       aes(x = week, 
           y = CountyName, 
           fill = weeklyCases)) +
  geom_tile() +
  scale_fill_gradient(low = "red", 
                      high = "green", 
                      limits = c(-692, 944)) +
  xlab("Time (in epid. calendar week)") +
  ylab("County") +
  labs(fill = TeX("$\\Delta$ weekly cases"))
ggsave("plots/dataVisualisation/heatmap_2022.pdf", 
       width = 23, height = 15, unit = "cm")


# Save pre-processed data -------------------------------------------------
# save data frame with 1-lag COVID-19 ID and smoothed winter 21/22
COVID_final <- COVID_data_week_smooth %>% 
  dplyr::select("CountyName", 
                "yw", 
                "weeklyCasesSum", 
                "weeklyCasesDiff", 
                "PopulationCensus16") %>% 
  rename(weeklyCases = weeklyCasesDiff)

# check if all necessary columns include 
COVID_final %>% colnames()

# save data set
write_csv(COVID_final,
          file = "data/COVID/ireland_covid_weekly.csv")

