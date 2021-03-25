

set.seed(42075)

# Load libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(zoo)
library(magrittr)
library(tidyr)
library(reshape2)
library(gridExtra)
library(pracma)
library(zoo)
library(mvtnorm)
library(GGally)
library(xtable)


##### Prepare COVID-19 data #####

# Read in data
covid19.df = read.csv(file="time_series_covid19_deaths_global.csv")
covid19_feb21.df = read.csv(file="time_series_covid19_deaths_global_090221.csv")
covid19_feb21_reduced.df = covid19_feb21.df[,c(1:2,191:387)]
# Left join subsequent data 
covid19_append.df = covid19.df %>% left_join(covid19_feb21_reduced.df, by=c("Province.State"="Province.State", "Country.Region"="Country.Region"))
covid19_append.df[is.na(covid19_append.df)] = 0

# Calculate daily cumulative sum by country
covid19_bycountry.df = aggregate(covid19_append.df[, 5:length(covid19_append.df)], by=list(covid19.df$Country.Region), FUN=sum)
covid19_bycountry.df = cbind(names(covid19_bycountry.df), t(covid19_bycountry.df))
covid19_bycountry.df[1,] = ifelse(covid19_bycountry.df[1,]=="United Kingdom", "UK", covid19_bycountry.df[1,])
col_names <- covid19_bycountry.df[1,]
colnames(covid19_bycountry.df) = col_names
covid19_bycountry.df = as.data.frame(covid19_bycountry.df[2:nrow(covid19_bycountry.df) ,])

# Manipulate date variables
covid19_bycountry.df = separate(as.data.frame(covid19_bycountry.df), Group.1, c("Month", "Day", "Year"))
covid19_bycountry.df$Month = sub(".", "", covid19_bycountry.df$Month)
covid19_bycountry.df$Year = paste0("20", covid19_bycountry.df$Year)
covid19_bycountry.df %<>% mutate(Date = as.Date(ISOdate(year = Year, month = Month, day = Day))) 
covid19_bycountry.df %<>% mutate_if(is.character, function(x) as.numeric(x)) %>%
  select(-c(Day, Month, Year)) %>%
  select(Date, everything()) 
covid19_bycountry.df %<>% mutate(Global = rowSums(covid19_bycountry.df[2:length(covid19_bycountry.df)]))

# Create dataset of number of new cases daily
covid19_bycountry_dailycounts.df = covid19_bycountry.df %>%
  mutate_if(is.numeric,  function(x) x- lag(x, default = NA)) 
covid19_bycountry_dailycounts.df = covid19_bycountry_dailycounts.df[-1 ,]
covid19_bycountry_dailycounts.df = melt(covid19_bycountry_dailycounts.df, id.vars = 'Date', variable.name = 'Country', value.name = "n")
# Correct negative numbers of events
covid19_bycountry_dailycounts.df = covid19_bycountry_dailycounts.df %>%
  mutate(n = ifelse(n < 0, 0, n)) %>%
  mutate(t = as.numeric(covid19_bycountry_dailycounts.df$Date - min(covid19_bycountry_dailycounts.df$Date)+1)) %>%
  mutate(Country = as.character(Country))

# Compute 5 day rolling average
scatter_events = as.data.frame(cbind(rep(covid19_bycountry_dailycounts.df$t, covid19_bycountry_dailycounts.df$n), 
                                     rep(as.character(covid19_bycountry_dailycounts.df$Country), covid19_bycountry_dailycounts.df$n)))
colnames(scatter_events) = c("t", "Country")
scatter_events %<>% mutate(t = as.numeric(as.character(t))) %>%
  mutate(adj = floor(runif(nrow(scatter_events))*7 - 3)) %>%
  mutate(t_adj = t + adj)
adj_daily_counts = scatter_events %>% 
  group_by(t_adj, Country) %>%
  tally() %>%
  arrange(Country, t_adj) %>%
  ungroup() %>%
  mutate(Date = as.Date(t_adj, origin = min(covid19_bycountry_dailycounts.df$Date)-1)) %>%
  mutate(Country = as.character(Country)) %>%
  select(Date, Country, n)
# Update daily counts
covid19_bycountry_dailycounts.df %<>% 
  select("Date", "Country") %>%
  left_join(adj_daily_counts, by=c("Country", "Date")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  filter(Date <= (max(covid19_bycountry.df$Date) - 3))
