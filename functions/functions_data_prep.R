##### Data functions #####

# Function to filter data 
filter_data_cp_peak = function(country, cps, i_star=NULL){

  # Filter data to country of interest
  modelling_data.df = covid19_bycountry_dailycounts.df %>%
    filter(Country == country) 
  
  # Start observation period once there is 10 events 
  t1_ind = as.numeric(which(cumsum(modelling_data.df$n) >=10)[1])
  
  # Calculate number of days in sample
  start_date = modelling_data.df$Date[t1_ind]

  if (!(country %in% c("China", "US", "Spain", "India"))){
    end_date = max(modelling_data.df$Date)
  } else if (country == "China"){
    end_date = as.Date("2020-04-13")
  } else if (country == "Spain"){
    end_date = as.Date("2020-06-15")
  } else if (country == "India"){
    end_date = as.Date("2020-06-12")
  } else if (country == "US"){
    end_date = as.Date("2020-06-21")
  } 
  max_T = as.numeric(difftime(end_date,  start_date, units="days")) + 1
  
  # Modelling data
  modelling_data.df %<>%
    filter(Date <= end_date) %>%
    mutate(time = as.numeric(difftime(Date, start_date, units="days") + 1))  %>%
    select(time, n) 

  # Calculate individual events and their previous events
  event_count = as.matrix(sapply(modelling_data.df, as.numeric))
  decay_times = lapply(1:max_T, function(t) as.numeric(t- event_count[event_count[,1] <t & event_count[,2] >0 ,1]))
  decay_times_counts = lapply(1:max_T, function(t) as.numeric(event_count[event_count [,1] <t & event_count[,2] >0 ,2]))
  
  # Create global variables
  assign("start_date", start_date, envir = .GlobalEnv)
  assign("event_count", event_count, envir = .GlobalEnv)
  
  if (length(cps)==1){
    
    # Upward trajectory    
    max_T_upward = cps[1]
    event_count_upward = event_count[t1_ind:(max_T_upward+t1_ind-1),]
    decay_times_upward = decay_times[1:max_T_upward]
    decay_times_counts_upward = decay_times_counts[1:max_T_upward]

    # Create global variables
    assign("event_count_upward", event_count_upward, envir = .GlobalEnv)
    assign("decay_times_upward", decay_times_upward, envir = .GlobalEnv)
    assign("decay_times_counts_upward", decay_times_counts_upward, envir = .GlobalEnv)
    
    # Downward trajectory   
    if (!(country %in% c("Brazil", "India"))){
      event_count_downward = event_count[(max_T_upward+t1_ind):(max_T+t1_ind-1),]
      decay_times_downward = decay_times[(max_T_upward+1):max_T]
      decay_times_counts_downward = decay_times_counts[(max_T_upward+1):max_T]
      assign("decay_times_downward", decay_times_downward, envir = .GlobalEnv)
      assign("decay_times_counts_downward", decay_times_counts_downward, envir = .GlobalEnv)
      assign("event_count_downward", event_count_downward, envir = .GlobalEnv)
      
    }
    
    # If doing cross validation
    if (!is.null(i_star)){
      if (i_star <= cps[1]){
        # Upward trajectory
        event_count_upward_val = event_count_upward[1:i_star,]
        decay_times_upward_val = decay_times_upward[1:i_star]
        decay_times_counts_upward_val = decay_times_counts_upward[1:i_star]
      } else {
        event_count_upward_val = event_count_upward
        decay_times_upward_val = decay_times_upward
        decay_times_counts_upward_val = decay_times_counts_upward
        
        if (!(country %in% c("Brazil", "India"))){
          # Downward trajectory
          if (i_star == (cps[1]+1)){
            event_count_downward_val = t(event_count_downward[1:(i_star-cps[1]),])
          } else {
            event_count_downward_val = event_count_downward[1:(i_star-cps[1]),]
          }
          decay_times_downward_val = decay_times_downward[1:(i_star-cps[1])]
          decay_times_counts_downward_val = decay_times_counts_downward[1:(i_star-cps[1])]
          # Create global variables
          assign("decay_times_downward_val", decay_times_downward_val, envir = .GlobalEnv)
          assign("decay_times_counts_downward_val", decay_times_counts_downward_val, envir = .GlobalEnv)
          assign("event_count_downward_val", event_count_downward_val, envir = .GlobalEnv)
        }
      }
      
      # Create global variables
      assign("decay_times_upward_val", decay_times_upward_val, envir = .GlobalEnv)
      assign("decay_times_counts_upward_val", decay_times_counts_upward_val, envir = .GlobalEnv)
      assign("event_count_upward_val", event_count_upward_val, envir = .GlobalEnv)
      
    }
  } 

  return()
}


