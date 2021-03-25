##### Data functions #####


# Function to filter data manually specifying change points
filter_data_cp_peak = function(country, cps, i_star=NULL){

  # Filter data to country of interest
  modelling_data.df = covid19_bycountry_dailycounts.df %>%
    filter(Country == country) 
  # Start observation period once there is 10 events 
  t1_ind = as.numeric(which(cumsum(modelling_data.df$n) >=10)[1])
  
  # Calculate number of days in sample
  start_date = modelling_data.df$Date[t1_ind]
  if (!(country %in% c("UK", "Brazil", "US"))){
    end_date = max(modelling_data.df$Date)
  } else if (country == "UK") {
    end_date = as.Date("2021-01-24")
  } else if (country == "Brazil") {
    end_date = as.Date("2021-01-07")
  } else if (country == "US") {
    end_date = as.Date("2021-01-12")
  }
  max_T = as.numeric(difftime(end_date,  start_date, units="days")) + 1
  

  # Modelling data
  modelling_data.df %<>%
    mutate(time = as.numeric(difftime(Date, start_date, units="days") + 1))  %>%
    select(time, n) 

  # Calculate individual events and their previous events
  event_count = as.matrix(sapply(modelling_data.df, as.numeric))
  decay_times = lapply(1:max_T, function(t) as.numeric(t- event_count[event_count[,1] <t & event_count[,2] >0 ,1]))
  decay_times_counts = lapply(1:max_T, function(t) as.numeric(event_count[event_count [,1] <t & event_count[,2] >0 ,2]))
  
  # Create global variables
  assign("start_date", start_date, envir = .GlobalEnv)
  assign("event_count", event_count, envir = .GlobalEnv)
  
  
  breaks = list()
  if (country=="Brazil"){
    breaks$new_end = 74
    breaks$cp_ind = 1
  } else if (country=="China"){
    breaks$new_end = 81
    breaks$cp_ind = 2
  } else if (country=="India"){
    breaks$new_end = 81
    breaks$cp_ind = 1
  } else if (country=="Spain"){
    breaks$new_end = 101
    breaks$cp_ind = 2
  } else if (country=="US"){
    breaks$new_end = 109
    breaks$cp_ind = 2
  } else if (country=="France"){
    breaks$new_end = 141
    breaks$cp_ind = 2
  } else if (country=="Germany"){
    breaks$new_end = 135
    breaks$cp_ind = 2
  } else if (country=="Italy"){
    breaks$new_end = 151
    breaks$cp_ind = 2
  } else if (country=="Sweden"){
    breaks$new_end = 130
    breaks$cp_ind = 2
  } else if (country=="UK"){
    breaks$new_end = 136
    breaks$cp_ind = 2
  } 
  
  if (!is.null(i_star)){
    cps_len = which.max(cps>i_star)
  } else {
    cps_len = length(cps)
  }

  dim = cps_len
  data = list()
  
  data$event_count = list() 
  data$decay_times = list()
  data$decay_times_counts = list()
  
  for (i in 1:dim){
    if (i==1){
      if (cps_len == i & !is.null(i_star)){
        data$event_count[[i]] = event_count[(t1_ind):(i_star+t1_ind-1),]
        data$decay_times[[i]] = decay_times[1:i_star]
        data$decay_times_counts[[i]] = decay_times_counts[1:i_star]
      } else {
        if (breaks$cp_ind == i){
          data$event_count[[i]] = event_count[(t1_ind):(breaks$new_end+t1_ind-1),]
          data$decay_times[[i]] = decay_times[1:breaks$new_end]
          data$decay_times_counts[[i]] = decay_times_counts[1:breaks$new_end]
        } else {
          data$event_count[[i]] = event_count[(t1_ind):(cps[i]+t1_ind-1),]
          data$decay_times[[i]] = decay_times[1:cps[i]]
          data$decay_times_counts[[i]] = decay_times_counts[1:cps[i]]
        }
      }
    } else {
      if (cps_len == i & !is.null(i_star)){
        data$event_count[[i]] = event_count[(t1_ind+cps[(i-1)]):(i_star+t1_ind-1),]
        data$decay_times[[i]] = decay_times[(cps[(i-1)]+1):i_star]
        data$decay_times_counts[[i]] = decay_times_counts[(cps[(i-1)]+1):i_star]
      } else {
        if (breaks$cp_ind == i){
          data$event_count[[i]] = event_count[(t1_ind+cps[(i-1)]):(breaks$new_end+t1_ind-1),]
          data$decay_times[[i]] = decay_times[(cps[(i-1)]+1):breaks$new_end]
          data$decay_times_counts[[i]] = decay_times_counts[(cps[(i-1)]+1):breaks$new_end]
        } else {
          data$event_count[[i]] = event_count[(t1_ind+cps[(i-1)]):(cps[i]+t1_ind-1),]
          data$decay_times[[i]] = decay_times[(cps[(i-1)]+1):cps[i]]
          data$decay_times_counts[[i]] = decay_times_counts[(cps[(i-1)]+1):cps[i]]
        }
      }
    }
  }

  assign("data", data, envir = .GlobalEnv)

  return()
}


