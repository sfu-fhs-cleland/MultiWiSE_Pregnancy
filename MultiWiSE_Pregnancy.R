# MultiWiSE_Pregnancy.R
# This code provides an example for how to use the MultiWiSE approach to estimate WFS exposure during pregnancy.
# The script follows the approach outlined in the manuscript 'Multiyear Wildfire Smoke Exposure (MultiWiSE) metrics: 
# A data-driven approach to characterizing episodic PM2.5 exposures for epidemiologic research on chronic health effects',
# and is adapted from the code provided in the 'MultiWiSE.R' script (https://github.com/sfu-fhs-cleland/MultiWiSE_Metrics_Manuscript). 
# The script was modified to use information on key maternal dates (date of conception and date of birth) and individual-level
# residential history to assemble a PM2.5 exposure profile and calculate the MultiWiSE metrics during pregnancy and each trimester.
# This approach uses data on PM2.5 in the 3 years prior to the date of birth to calculate the counterfactual value and modified z-scores,
# which are then used to estimate WFS and non-WFS PM2.5 during pregnancy and calculate the metrics. 

####################################################################################################################################
############################################# PACKAGES & WORKING DIRECTORY #########################################################
####################################################################################################################################

#### Install/Load packages ####
packages <- c(
  "dplyr", "data.table",
  "lubridate", "ggplot2", "tidyr",
  "readr", "patchwork"
)

installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install)) {
  install.packages(to_install)
}

lapply(packages, library, character.only = TRUE)
rm(packages, installed, to_install)

#### Set a working directory ####
# setwd("<working_directory>")
setwd('/Users/stephaniecleland/Library/CloudStorage/OneDrive-SimonFraserUniversity(1sfu)/Research/Characterizing Prolonged WFS/MultiWiSE Pregnancy')

####################################################################################################################################
############################################ READ IN & CHECK MATERNAL DATA #########################################################
####################################################################################################################################
## To calculate the MultiWiSE metrics during pregnancy, the user must provide two files:
##  1) A .csv file with the date of conception and date of birth for each individual:
##     File must have 3 columns: IndividualID, DateOfConception, and DateOfBirth. There must be only one row per IndividualID 
##     and DateOfConception and DateOfBirth must be provided as YYYY-MM-DD.
##  2) A .csv file with the maternal residential history:
##     File must have 4 columns: IndividualID, GeoID, StartDate, EndDate. The StartDate and EndDate represent the days the 
##     individual moved in and out of the GeoID and must be provided as YYYY-MM-DD. It must provide residential data from 
##     3 years prior to the DateOfBirth through the DateOfBirth with no gaps or overlap in the Dates.
## The files used in this script provide an example for what these should look like.

# Set length of exposure window to use for calculating counterfactual and modified z-scores
exp_wind_length = 3

### Read in the files
maternal_dates <- read.csv('./maternal_dates.csv')
maternal_dates <- maternal_dates %>%
  mutate(
    DateOfConception = as.Date(DateOfConception),
    DateOfBirth   = as.Date(DateOfBirth)
  ) %>% unique()

residential_history <- read.csv('./residential_history.csv')
residential_history <- residential_history %>%
  mutate(
    StartDate = as.Date(StartDate),
    EndDate   = as.Date(EndDate)
  ) %>% unique()

### Run some basic sanity checks on both files
# Check for required columns
if (!all(c("IndividualID", "GeoID", "StartDate", "EndDate") %in% names(residential_history))) {
  stop(sprintf("Error: residential_history must contain the following columns: IndividualID, GeoID, StartDate, EndDate"))
}
if (!all(c("IndividualID", "DateOfConception", "DateOfBirth") %in% names(maternal_dates))) {
  stop(sprintf("Error: maternal_dates must contain the following columns: IndividualID, DateOfConception, DateOfBirths"))
}
# Check for missing data
if (any(is.na(residential_history)) || any(is.na(maternal_dates))) {
  stop(sprintf(
    "Error: There is missing data in either or both of the files."
  ))
}

# Check that IDs are present in each file
res_ids <- unique(residential_history$IndividualID[!is.na(residential_history$GeoID)])
dat_ids <- unique(maternal_dates$IndividualID[!is.na(maternal_dates$DateOfConception) & !is.na(maternal_dates$DateOfBirth)])

missing_in_dat <- setdiff(res_ids, dat_ids)
missing_in_res <- setdiff(dat_ids, res_ids)
if (length(missing_in_dat) > 0 || length(missing_in_res) > 0) {
  stop(sprintf(
    "Error: There is a mismatch between the IndividualIDs in the files."
  ))
}

# Check that residential history covers required period (3 years prior to date of birth through date of birth)
res_hist_range <- residential_history %>% 
  group_by(IndividualID) %>% 
  summarise(Start = min(StartDate), End = max(EndDate))
check_res_hist <- merge(maternal_dates, res_hist_range, by = 'IndividualID')
if (any(check_res_hist$Start > check_res_hist$DateOfBirth - years(exp_wind_length)) || 
    any(check_res_hist$End < check_res_hist$DateOfBirth)) {
  stop(sprintf(
    "Error: The residential history does not provide data for the required period."
  ))
}

# Check that dates are in chronological order
if (any(residential_history$EndDate <= residential_history$StartDate)) {
  stop(sprintf(
    "Error: One or more of the EndDates in the residential history file occur before the StartDates"
  ))
}
if (any(maternal_dates$DateOfBirth <= maternal_dates$DateOfConception)) {
  stop(sprintf(
    "Error: One or more of the DateOfBirths in the residential history file occur before the DateOfConceptions"
  ))
}

# Check that there are no gaps or overlaps in residential history
span_check <- residential_history %>% group_by(IndividualID) %>%
  summarize(total_span = max(EndDate) - min(StartDate),
            sum_span = sum(EndDate - StartDate) + (n()-1))
if (any(span_check$sum_span != span_check$sum_span)) {
  stop(sprintf(
    "Error: There are gaps and/or overlaps present in the residential history file."
  ))
}
####################################################################################################################################
############################################# READ IN & PREPARE DAILY PM2.5 DATA ###################################################
####################################################################################################################################
### Read in daily PM2.5 data, combine into one dataset, and aggregate to epiweeks
# NOTE: Structured to work with the publicly available CanOSSEM data (https://www.sfu.ca/fhs/canossem.html)

all_geoIDs <- unique(residential_history$GeoID) # List of all postal codes in residential history file
all_years <- year(min(residential_history$StartDate)):year(max(residential_history$EndDate)) # List of all years in residential history file

pm25_filepath <- './Daily PM25 Data/'

pm25_daily <- data.table()

for (year in all_years) {
  # Read in data
  load(paste0(pm25_filepath,'PC_CanOSSEM_BC_',year,'_WFS.RData'))
  
  # Subset to postal codes in residential history
  PC_CanOSSEM_WFS <- PC_CanOSSEM_WFS[PC_CanOSSEM_WFS$PC %in% all_geoIDs, c('DATE','PC','PC_PM2.5')]
  
  pm25_daily <- rbind(pm25_daily, PC_CanOSSEM_WFS)
}

### Aggregate daily PM2.5 data to weekly level
# Create a complete sequence of dates from start to end of PM2.5 data
dates <- seq(
  from = min(pm25_daily$DATE),    
  to   = max(pm25_daily$DATE),      
  by   = "days"                     
)

# Build a helper data.frame mapping each date to its epiweek number
weeks_df <- data.frame(
  date = dates,
  epiweek_num = epiweek(as.POSIXct(dates))
)

# Compute a continuous epiweek variable
weeks_df$cont_epiweek_num <- rleid(weeks_df$epiweek_num) - 1

# Remove incomplete weeks at beginning and end of period
weeks_df <- weeks_df %>% 
  group_by(cont_epiweek_num) %>% 
  mutate(n = n(),
         start_date = min(date),
         end_date = max(date),
         epiweek = as.integer(cont_epiweek_num)) %>% 
  filter(n == 7)

# Add the epiweek data to the daily PM2.5 data
pm25_daily <- merge(pm25_daily, weeks_df[, c('date', 'epiweek', 'start_date', 'end_date')], by.x = 'DATE', by.y = 'date') 

names(pm25_daily) <- c('date', 'geoid', 'pm25', 'epiweek', 'start_date', 'end_date')

# Aggregate to weekly level
pm25_weekly <- pm25_daily %>%
  group_by(geoid, epiweek) %>%
  summarize(
    weekly_n_days    = n(),                 # How many days had data for this week
    weekly_avg_pm25  = mean(pm25),          # Average of daily PM2.5 values
    weekly_sum_pm25  = sum(pm25),           # Total sum of daily PM2.5 values
    year_start       = min(year(date)),     # Year of the first day in the epiweek
    start_date       = unique(start_date),  # First calendar date in the epiweek
    end_date         = unique(end_date),    # Last calendar date in the epiweek
  )

####################################################################################################################################
############################ DEFINE FUNCTIONS FOR PROCESSING DATA & GENERATING METRICS #############################################
####################################################################################################################################
# Function to build exposure profiles for an individual - exposure profile spans from 3 years prior to the date of birth through the date of birth
build_ind_pm25_profile <- function(ind_residential_history, # data.frame with individual's residential history
                                   ind_maternal_dates, # data.frame with individual's maternal dates
                                   ind_pm25_data, # data.frame with the weekly PM2.5 data for the individual's residential history
                                   exp_wind_length = 3 # number of years before birth date to use to calculate counterfactual & z-scores
) {
  exp_start_date <- ind_maternal_dates$DateOfBirth - years(exp_wind_length) # Set exposure window start date to 3 years before birth date
  exp_end_date <- ind_maternal_dates$DateOfBirth # Set exposure window end date to birth date
  
  res_dt <- as.data.table(ind_residential_history)[order(StartDate)]
  
  pm25_list <- list(); ctr <- 1L
  
  # Use residential history to generate weekly time series of PM2.5
  for (i in seq_len(nrow(res_dt))) {
    geoid    <- res_dt$GeoID[i]
    geo_start <- res_dt$StartDate[i]
    geo_end   <- res_dt$EndDate[i]
    
    overlap_start <- max(geo_start, exp_start_date)
    overlap_end   <- min(geo_end, exp_end_date)
    if (overlap_end < overlap_start) next
    
    pm25 <- data.table(ind_pm25_data[ind_pm25_data$geoid == geoid &
                                       ind_pm25_data$end_date >= overlap_start & 
                                       ind_pm25_data$start_date <= overlap_end,])
    
    if (nrow(pm25) == 0) next
    
    pm25[, `:=`(
      partial_interval_start = pmax(start_date, overlap_start),
      partial_interval_end   = pmin(end_date,   overlap_end)
    )]
    pm25[, days_overlap     := as.numeric(partial_interval_end - partial_interval_start) + 1L]
    pm25[, days_overlap     := ifelse(days_overlap < 7, days_overlap, weekly_n_days)]
    pm25[, fraction_of_week := days_overlap / weekly_n_days]
    pm25[, partial_pm25_weekly_sum := weekly_sum_pm25 * fraction_of_week]
    pm25[, partial_n_days          := weekly_n_days * fraction_of_week]
    
    pm25_list[[ctr]] <- pm25[, .(
      epiweek, start_date, end_date, year_start,
      partial_pm25_weekly_sum,
      partial_n_days, geoid
    )]
    ctr <- ctr + 1L
  }
  
  if (ctr == 1L) {
    stop("Error: No overlapping data found for this individual's history.", call. = FALSE)
  }
  
  combine_partial <- data.table::rbindlist(pm25_list, fill = TRUE)
  
  # Handle weeks split between two locations
  weekly_exposure <- combine_partial[, .(
    weekly_sum_pm25 = sum(partial_pm25_weekly_sum, na.rm = TRUE),
    n_days          = sum(partial_n_days,          na.rm = TRUE),
    geoid           = geoid[which.max(partial_n_days)]
  ), by = .(epiweek, start_date, end_date, year_start)][order(start_date)]
  
  weekly_exposure[, weekly_avg_pm25 := weekly_sum_pm25 / n_days]
  weekly_exposure <- weekly_exposure[n_days > 0]
  
  # Drop incomplete first and last epi-weeks
  if (nrow(weekly_exposure) > 0 && weekly_exposure$n_days[1] < 7){
    weekly_exposure <- weekly_exposure[-1]
  }
  if (nrow(weekly_exposure) > 0 && weekly_exposure$n_days[nrow(weekly_exposure)] < 7){
    weekly_exposure <- weekly_exposure[-nrow(weekly_exposure)]
  }
  
  weekly_exposure$epiweek <- rleid(weekly_exposure$epiweek)
  
  # Add flag for wildfire season and calculate modified Z-score
  weekly_exposure <- weekly_exposure %>%
    arrange(epiweek) %>%
    mutate(
      is_wildfire_season = (end_date >= as.Date(paste0(year_start, "-05-01")) &
                              start_date < as.Date(paste0(year_start, "-11-01"))),
      modified_z = calc_modified_z(weekly_sum_pm25)
    )
  
  weekly_exposure
}

# Function to calculate modified Z-Score
calc_modified_z <- function(values) {
  med <- median(values, na.rm = TRUE)
  mad_val <- mad(values, constant = 1.4826, na.rm = TRUE)
  if (mad_val == 0) return(rep(NA_real_, length(values)))
  (values - med) / mad_val
}

# Function to calculate counterfactual
calc_counterfactual <- function(
    ind_df    # data.frame with weekly individual-level PM2.5 data for the 3-year exposure profile
) {
  ind_df %>%
    filter(!is.na(modified_z),
           dplyr::between(modified_z, -2, 2),
           n_days == 7) %>%
    pull(weekly_sum_pm25) %>% median(na.rm = TRUE) %>% replace_na(0)
}

# Function to identify smoke episodes
identify_smoke_episodes <- function(
    is_smoke_impacted,                 # Vector identifying if an epiweek is WFS-impacted or not (1 or 0) 
    cont_week,                         # Vector of continuous epiweeks
    weekly_sum_wfs_pm25_values,        # Vector of the sum of WFS PM2.5 for each epiweek
    episode_threshold,                 # Value that must be exceeded to count as an episode (0 for WFS episode, 250 for severe episode)
    num_weeks_in_episode,              # Number of WFS-impacted weeks that need to be present to count as an episode (2 for WFS episode, 1 for severe episode)
    max_no_smoke_gap,                  # Maximum number of continuous non-WFS-impacted weeks that can be present and still count as an episode (3 weeks)
    weeks_needed_above_threshold       # Number of weeks above the threshold needed to count as an episode (1 week)
) {
  # Total number of epiweeks
  n <- length(is_smoke_impacted)
  
  # Initialize output vector (episode IDs), current episode counter, and flag
  epi_id <- integer(n)
  cur_id <- 0
  in_epi <- FALSE
  
  # Trackers within each candidate episode
  smoke_wks <- 0      # total smoke-impacted weeks
  wks_above <- 0      # weeks where PM2.5 > threshold
  gap       <- 0      # current non-smoke gap length
  bridge    <- integer(0)  # indices of non-smoke weeks within the gap
  
  # Helper to abandon a failed episode (reset its epi_id entries to 0)
  fail <- function(id) epi_id[epi_id == id] <<- 0
  
  # Iterate week by week
  for (i in seq_len(n)) {
    
    # Break episode if non-consecutive week index
    if (i > 1 && cont_week[i] - cont_week[i - 1] > 1) {
      # If in an episode that doesn’t meet criteria, drop it
      if (in_epi &&
          (smoke_wks < num_weeks_in_episode ||
           wks_above < weeks_needed_above_threshold)) {
        fail(cur_id)
      }
      # Reset all episode trackers
      in_epi    <- FALSE
      smoke_wks <- 0
      wks_above <- 0
      gap       <- 0
      bridge    <- integer(0)
    }
    
    # Handle a smoke-impacted week
    if (is_smoke_impacted[i] == 1) {
      
      # Start a new episode if not already in one
      if (!in_epi) {
        cur_id <- cur_id + 1
        in_epi <- TRUE
        smoke_wks <- 0
        wks_above <- 0
        gap <- 0
      }
      
      # Any bridged non-smoke weeks now belong to this episode
      if (length(bridge) > 0) {
        epi_id[bridge] <- cur_id
        bridge <- integer(0)
      }
      
      # Update counters
      smoke_wks <- smoke_wks + 1
      if (weekly_sum_wfs_pm25_values[i] > episode_threshold) {
        wks_above <- wks_above + 1
      }
      
      epi_id[i] <- cur_id
      gap <- 0  # reset gap
      
      # Handle a non-smoke week within an ongoing episode
    } else if (in_epi) {
      gap <- gap + 1
      bridge <- c(bridge, i)  # tentatively include this week
      
      # If gap is too long, close out the episode
      if (gap > max_no_smoke_gap) {
        
        # Abandon if criteria not met
        if (smoke_wks < num_weeks_in_episode ||
            wks_above < weeks_needed_above_threshold) {
          fail(cur_id)
        }
        in_epi <- FALSE
        bridge <- integer(0)
      }
    }
  }
  
  # Final check at end of series
  if (in_epi &&
      (smoke_wks < num_weeks_in_episode ||
       wks_above < weeks_needed_above_threshold)) {
    fail(cur_id)
  }
  
  # Trim trailing non-smoke weeks after last smoke week in each episode
  for (id in setdiff(unique(epi_id), 0)) {
    idx <- which(epi_id == id)
    smoke_idx <- idx[is_smoke_impacted[idx] == 1]
    epi_id[idx[idx > max(smoke_idx)]] <- 0
  }
  
  # Reset episode IDs to 1,2,3, etc. in order of discovery
  uniq <- setdiff(unique(epi_id), 0)
  for (k in seq_along(uniq)) {
    epi_id[epi_id == uniq[k]] <- k
  }
  epi_id
}

# Function to process data for a given individual to generate WFS and non-WFS PM2.5 and identify WFS episodes during pregnancy
process_ind_pm25_profile <- function(
    ind_df,                                     # data.frame with weekly individual-level PM2.5 data for the duration of pregnancy
    ind_counterfactual,                         # counterfactual value for individual based on 3-year exposure profile
    threshold_episode               = 0,        # Value that must be exceeded to count as a WFS episode (0 ug/m3)
    num_wfs_weeks_in_episode        = 2,        # Number of WFS-impacted weeks that need to be present to count as a WFS episode (2 weeks)
    threshold_severe_episode        = 250,      # Value that must be exceeded to count as a severe episode (250 ug/m3)
    num_wfs_weeks_in_severe_episode = 1,        # Number of WFS-impacted weeks that need to be present to count as a severe episode (1 weeks)
    max_no_smoke_gap                = 3,        # Maximum number of continuous non-WFS-impacted weeks that can be present and still count as an episode (3 weeks)
    weeks_needed_above_threshold    = 1,        # Maximum number of continuous non-WFS-impacted weeks that can be present and still count as an episode (3 weeks)
    min_smoke_thresh                = 0.005     # Minimum value for WFS PM2.5 - used to remove near-zero estimates
) {
  
  # Calculate cumulative total PM2.5
  ind_df <- ind_df %>%
    arrange(epiweek) %>%
    mutate(
      total_cumulative = cumsum(weekly_sum_pm25)
    )
  
  # Separate WFS and non-WFS components
  ind_df <- ind_df %>%
    mutate(
      non_wfs_weekly_sum_pm25 = ifelse(modified_z > 2 & is_wildfire_season,
                                       ind_counterfactual, weekly_sum_pm25),
      non_wfs_weekly_avg_pm25 = ifelse(modified_z > 2 & is_wildfire_season,
                                       ind_counterfactual / 7, weekly_avg_pm25),
      wfs_weekly_sum_pm25 = weekly_sum_pm25 - non_wfs_weekly_sum_pm25,
      wfs_weekly_avg_pm25 = weekly_avg_pm25 - non_wfs_weekly_avg_pm25,
      non_wfs_cumulative = cumsum(non_wfs_weekly_sum_pm25),
      wfs_cumulative = cumsum(wfs_weekly_sum_pm25)
    ) %>%
    mutate(
      wfs_weekly_avg_pm25 = if_else(wfs_weekly_avg_pm25 < min_smoke_thresh, 0, wfs_weekly_avg_pm25),
      wfs_weekly_sum_pm25 = if_else(wfs_weekly_sum_pm25 < min_smoke_thresh, 0, wfs_weekly_sum_pm25),
      is_smoke_impacted = as.integer(wfs_weekly_avg_pm25 > 0)
    )
  
  # Identify WFS episodes & severe episodes
  norm_ids <- identify_smoke_episodes(
    ind_df$is_smoke_impacted, ind_df$epiweek, ind_df$wfs_weekly_sum_pm25,
    episode_threshold = threshold_episode,
    num_weeks_in_episode = num_wfs_weeks_in_episode,
    max_no_smoke_gap = max_no_smoke_gap,
    weeks_needed_above_threshold = weeks_needed_above_threshold
  )
  
  sev_ids <- identify_smoke_episodes(
    ind_df$is_smoke_impacted, ind_df$epiweek, ind_df$wfs_weekly_sum_pm25,
    episode_threshold = threshold_severe_episode,
    num_weeks_in_episode = num_wfs_weeks_in_severe_episode,
    max_no_smoke_gap = max_no_smoke_gap,
    weeks_needed_above_threshold = weeks_needed_above_threshold
  )
  
  epi_id <- integer(nrow(ind_df))
  cur    <- 0; in_epi <- FALSE
  for (i in seq_len(nrow(ind_df))) {
    if (norm_ids[i] > 0 || sev_ids[i] > 0) {
      if (!in_epi) { cur <- cur + 1; in_epi <- TRUE }
      epi_id[i] <- cur
    } else in_epi <- FALSE
  }
  
  ind_df$episode_id <- epi_id
  ind_df$severe_episode_id <- sev_ids
  
  ind_df$epiweek <- rleid(ind_df$epiweek)
  
  # Return processed data
  ind_df
}

# Function to number metrics from #1-12
metric_numbering <- c(
  "Cumulative WFS PM2.5 (mg/m³)"           = "1. Cumulative WFS PM₂.₅ (mg/m³)",
  "WFS Fraction (%)"                       = "2. WFS Fraction (%)",
  "Average WFS PM2.5 (µg/m³)"              = "3. Average WFS PM₂.₅ (µg/m³)",
  "Any WFS (# weeks)"                      = "4. Any WFS (# weeks)",
  "WFS PM2.5 > 5 µg/m³ (# weeks)"          = "5. WFS PM₂.₅ > 5 µg/m³ (# weeks)",
  "Total PM2.5 > 25 µg/m³ (# weeks)"       = "6. Total PM₂.₅ > 25 µg/m³ (# weeks)",
  "WFS Episodes (# episodes)"              = "7. WFS Episodes (# episodes)",
  "Severe Episodes (# episodes)"           = "8. Severe Episodes (# episodes)",
  "Longest Episode (# weeks)"              = "9. Longest Episode (# weeks)",
  "Worst Episode (µg/m³)"                  = "10. Worst Episode (µg/m³)",
  "WFS from Severe Episodes (%)"           = "11. WFS from Severe Episodes (%)",
  "Average Recovery (# weeks)"             = "12. Average Recovery (# weeks)"
)
apply_metric_numbers <- function(df) {
  rename_with(df, ~ metric_numbering[.x] %||% .x)
}

# Function to compute the 12 MultiWiSE metrics for a given individual
compute_metrics_ind <- function(
    ind_df,                   # data.frame with weekly individual-level exposure data for period of interest (e.g., pregnancy, specific trimester, etc.)
    threshold_low    = 5,     # Threshold used for metric #5
    threshold_high   = 25,    # Threshold used for metric #6
    min_smoke_thresh = 0.005
) {
  
  # Conversion for ug to mg
  MICRO_TO_MILLI <- 1/1000
  
  ind_df <- ind_df %>%
    mutate(
      wfs_weekly_avg_pm25 = if_else(wfs_weekly_avg_pm25 < min_smoke_thresh, 0, wfs_weekly_avg_pm25),
      wfs_weekly_sum_pm25 = if_else(wfs_weekly_sum_pm25 < min_smoke_thresh, 0, wfs_weekly_sum_pm25)
    )
  
  # 1. Cumulative WFS PM2.5
  wfs_cumulative <- sum(ind_df$wfs_weekly_sum_pm25, na.rm = TRUE) * MICRO_TO_MILLI
  
  # 2. WFS Fraction
  wfs_fraction <- 100 * wfs_cumulative / (sum(ind_df$weekly_sum_pm25, na.rm = TRUE) * MICRO_TO_MILLI)
  
  # 3. Average WFS PM2.5
  mean_smoke_week <- ifelse(sum(ind_df$wfs_weekly_avg_pm25 > 0, na.rm = TRUE) > 0,
    mean(ind_df$wfs_weekly_avg_pm25[ind_df$wfs_weekly_avg_pm25 > 0]), 0)
  
  # 4. Any WFS
  total_smoke_weeks <- sum(ind_df$wfs_weekly_avg_pm25 > 0, na.rm = TRUE)
  
  # 5. WFS PM2.5 > 5
  weeks_wfs_over_5 <- sum(ind_df$wfs_weekly_avg_pm25 > threshold_low, na.rm = TRUE)
  
  # 6. Total PM2.5 > 25
  weeks_over_25 <- sum(ind_df$wfs_weekly_avg_pm25 > 0 & ind_df$weekly_avg_pm25 > threshold_high, na.rm = TRUE)
  
  # 7. WFS Episodes
  n_normal <- max(ind_df$episode_id, na.rm = TRUE)
  
  # 8. Severe Episodes
  n_severe <- max(ind_df$severe_episode_id, na.rm = TRUE)
  
  # 9. Longest Episode
  ep_summary_all <- ind_df %>%
    filter(episode_id > 0) %>%
    group_by(episode_id) %>%
    summarise(
      active_weeks = n(),
      .groups = "drop"
    )
  longest_episode_len <- ifelse(nrow(ep_summary_all) > 0, max(ep_summary_all$active_weeks), 0)
  
  # 10. Worst Episode
  ep_summary_smk <- ind_df %>%
    filter(episode_id > 0, is_smoke_impacted == 1) %>%
    group_by(episode_id) %>%
    summarise(
      avg_pm25 = mean(wfs_weekly_avg_pm25),
      .groups = "drop"
    )
  worst_exposure <- ifelse(nrow(ep_summary_smk) > 0, max(ep_summary_smk$avg_pm25), 0)
  
  # 11. WFS from Severe Episodes
  severe_pm25_cum <- ind_df %>%
    filter(severe_episode_id > 0, is_smoke_impacted == 1) %>%
    summarise(sum(wfs_weekly_sum_pm25, na.rm = TRUE)) %>% pull()
  severe_pm25_cum <- severe_pm25_cum * MICRO_TO_MILLI
  pct_wfs_from_severe <- ifelse(wfs_cumulative > 0, 100 * severe_pm25_cum / wfs_cumulative, 0)
  
  # 12. Average Recovery 
  period_start <- min(ind_df$start_date, na.rm = TRUE)
  period_end   <- max(ind_df$end_date,   na.rm = TRUE)
  if (n_normal == 0) {    # No episodes -> recovery is the full window
    avg_time_between <- round(as.numeric(period_end - period_start + 1) / 7, 2)
  } else {
    episode_info <- ind_df %>%
      filter(episode_id > 0, is_smoke_impacted == 1) %>%
      group_by(episode_id) %>%
      summarise(
        start_date = min(start_date),
        end_date = max(end_date),
        .groups = "drop"
      ) %>%
      arrange(start_date)
    gap_starts <- c(period_start, episode_info$end_date + 1) # first day after each episode
    gap_ends <- c(episode_info$start_date - 1, period_end) # last day before next episode
    gaps_days <- pmax(0, as.numeric(gap_ends - gap_starts + 1)) 
    avg_time_between <- round(mean(gaps_days) / 7, 2)
  }
  
  # Create a table (tibble) that assigns the metrics to the official metric name
  numbered <- tibble(
    `Cumulative WFS PM2.5 (mg/m³)` = round(wfs_cumulative, 4),
    `WFS Fraction (%)` = round(wfs_fraction, 2),
    `Average WFS PM2.5 (µg/m³)`= round(mean_smoke_week, 2),
    `Any WFS (# weeks)` = total_smoke_weeks,
    `WFS PM2.5 > 5 µg/m³ (# weeks)` = weeks_wfs_over_5,
    `Total PM2.5 > 25 µg/m³ (# weeks)` = weeks_over_25,
    `WFS Episodes (# episodes)` = n_normal,
    `Severe Episodes (# episodes)` = n_severe,
    `Longest Episode (# weeks)` = longest_episode_len,
    `Worst Episode (µg/m³)` = round(worst_exposure, 2),
    `WFS from Severe Episodes (%)` = round(pct_wfs_from_severe, 1),
    `Average Recovery (# weeks)` = avg_time_between
  ) |> apply_metric_numbers()
  
  # Additional variables (total, non-WFS, and counterfactual) to be included in results
  extras <- tibble::tibble(
    `Cumulative_Total_PM25`     = round(sum(ind_df$non_wfs_weekly_sum_pm25,  na.rm = TRUE) * MICRO_TO_MILLI, 2),
    `Average_Total_PM25`        = round(mean(ind_df$weekly_avg_pm25, na.rm = TRUE), 2),
    `Cumulative_NonWFS_PM25`    = round(sum(ind_df$weekly_sum_pm25, na.rm = TRUE) * MICRO_TO_MILLI, 2),
    `Average_NonWFS_PM25`       = round(mean(ind_df$non_wfs_weekly_avg_pm25, na.rm = TRUE), 2),
    `Counterfactual_Value`      = round(ind_counterfactual / 7, 2)
  )
  
  dplyr::bind_cols(numbered, extras)
}

####################################################################################################################################
############################ PROCESS INDIVIDUAL DATA & GENERATE METRICS ############################################################
####################################################################################################################################

all_ids <- unique(maternal_dates$IndividualID)

all_exp_profile_preg <- data.table()
all_metrics_preg <- data.table()

# For each individual: 1) assemble exposure profile, 2) estimate WFS and non-WFS PM2.5 and identify episodes, and 3) calculate metrics
for (ID in all_ids) {
  ## Subset datasets to individual
  ind_residential_history <- residential_history[residential_history$IndividualID == ID,]
  ind_maternal_dates <- maternal_dates[maternal_dates$IndividualID == ID,] 
  ind_pm25_data <- pm25_weekly[pm25_weekly$geoid %in% ind_residential_history$GeoID,] 
  
  ## Assemble 3-year exposure profile and calculate modified Z-score
  ind_exp_profile <- build_ind_pm25_profile(ind_residential_history, ind_maternal_dates, ind_pm25_data)
  
  ## Calculate counterfactual using 3-year exposure profile
  ind_counterfactual <- calc_counterfactual(ind_exp_profile)
  
  # Subset processed data to pregnancy period
  ind_exp_profile_preg <- ind_exp_profile[start_date >= ind_maternal_dates$DateOfConception &
                                            end_date <= ind_maternal_dates$DateOfBirth]
  
  ## Process exposure during pregnancy - calculate WFS and non-WFS PM2.5 and identify WFS episodes
  ind_processed_preg <- process_ind_pm25_profile(ind_exp_profile_preg,
                                                 ind_counterfactual)
  ind_processed_preg$IndividualID <- rep(ID, nrow(ind_processed_preg))
  
  ## Add flags for trimesters
  tri_1_enddate <- ind_maternal_dates$DateOfConception + weeks(13) # 13 weeks from date of conception
  tri_2_enddate <- ind_maternal_dates$DateOfConception + weeks(27) # 27 weeks from date of conception
  ind_processed_preg$trimester <- ifelse(ind_processed_preg$end_date <= tri_1_enddate,
                                         'Tri1',
                                         ifelse(ind_processed_preg$start_date <= tri_1_enddate,
                                                'Tri1/Tri2',
                                                ifelse(ind_processed_preg$end_date <= tri_2_enddate,
                                                       'Tri2',
                                                       ifelse(ind_processed_preg$start_date <= tri_2_enddate,
                                                              'Tri2/Tri3',
                                                              'Tri3'))))
  # For weeks with two trimesters, assign to trimester with most overlap
  ind_processed_preg[trimester == 'Tri1/Tri2',]$trimester <- ifelse(abs(ind_processed_preg[trimester == 'Tri1/Tri2',]$start_date - tri_1_enddate) > 3,
                                                                    'Tri1', 'Tri2')
  ind_processed_preg[trimester == 'Tri2/Tri3',]$trimester <- ifelse(abs(ind_processed_preg[trimester == 'Tri2/Tri3',]$start_date - tri_2_enddate) > 3,
                                                                    'Tri2', 'Tri3')
  
  ## Calculate metrics for periods of interest
  # Metrics for entire pregnancy
  metrics_preg <- compute_metrics_ind(ind_processed_preg) 
  metrics_preg$IndividualID <- ID; metrics_preg$window <- 'Pregnancy'
  
  # # Metrics for trimester 1
  metrics_tri1 <- compute_metrics_ind(ind_processed_preg[trimester == 'Tri1'])
  metrics_tri1$IndividualID <- ID; metrics_tri1$window <- 'Trimester1'

  # Metrics for trimester 2
  metrics_tri2 <- compute_metrics_ind(ind_processed_preg[trimester == 'Tri2'])
  metrics_tri2$IndividualID <- ID; metrics_tri2$window <- 'Trimester2'

  # Metrics for trimester 3
  metrics_tri3 <- compute_metrics_ind(ind_processed_preg[trimester == 'Tri3'])
  metrics_tri3$IndividualID <- ID; metrics_tri3$window <- 'Trimester3'

  ## Add data to master datasets
  all_exp_profile_preg <- bind_rows(all_exp_profile_preg, ind_processed_preg)
  all_metrics_preg <- bind_rows(all_metrics_preg, metrics_preg, metrics_tri1, metrics_tri2, metrics_tri3)
  
}

# Save processed individual-level data to a .csv file
write_csv(all_exp_profile_preg[,c('IndividualID',
                           'epiweek',
                           'n_days',
                           'start_date',
                           'end_date',
                           'trimester',
                           'weekly_avg_pm25',
                           'wfs_weekly_avg_pm25',
                           'non_wfs_weekly_avg_pm25',
                           'weekly_sum_pm25',
                           'wfs_weekly_sum_pm25',
                           'non_wfs_weekly_sum_pm25',
                           'modified_z',
                           'is_wildfire_season',
                           'is_smoke_impacted',
                           'episode_id',
                           'severe_episode_id')],
          './Ind_Pregnancy_Weekly_PM25_Episode_Estimates.csv')

# Save individual-level metrics for pregnancy and each trimester to .csv file
all_metrics_preg <- all_metrics_preg[,c(18:19,1:17)]
names(all_metrics_preg) <- c('IndividualID', 'Window', '1_Cumulative_WFS_PM25', 
                            '2_WFS_Fraction', 
                            '3_Average_WFS_PM25',
                            '4_Any_WFS' ,
                            '5_WFS_PM25_exceeds_5',
                            '6_Total_PM25_exceeds_25',
                            '7_WFS_Episodes',
                            '8_Severe_Episodes',
                            '9_Longest_Episode',
                            '10_Worst_Episode',
                            '11_WFS_from_Severe_Episodes',
                            '12_Average_Recovery',
                            'Cumulative_Total_PM25',
                            'Average_Total_PM25',
                            'Cumulative_NonWFS_PM25',
                            'Average_NonWFS_PM25',
                            'Counterfactual_Value'
)
write_csv(all_metrics_preg, "./Ind_Pregnancy_MultiWiSE_Metrics.csv")


####################################################################################################################################
######################################### GENERATE FIGURES #########################################################################
####################################################################################################################################

# Function to generate cumulative PM2.5 plot
plot_cumulative_pm25_by_month <- function(ind_df, y_limits = NULL, line_thin = 0.25) {
  
  ind_df <- ind_df %>% 
    mutate(total_cumulative = total_cumulative / 1000,
           non_wfs_cumulative = non_wfs_cumulative / 1000, 
           wfs_cumulative = wfs_cumulative / 1000)
  
  if (is.null(y_limits)) { 
    y_limits <- c(0, max(c(ind_df$total_cumulative, ind_df$non_wfs_cumulative, ind_df$wfs_cumulative), na.rm = TRUE) * 1.05) 
    }
  
  ggplot(ind_df, aes(x = start_date)) +
    geom_line(aes(y = total_cumulative,       color = "Total"),    linewidth = 1) +
    geom_line(aes(y = non_wfs_cumulative, color = "Non-WFS"), linewidth = 1) +
    geom_line(aes(y = wfs_cumulative,  color = "WFS"),     linewidth = 1) +
    scale_color_manual(
      values = c("Total" = "black", "Non-WFS" = "#B0C6DF", "WFS" = "#FF8572"),
      name   = NULL
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b-%y",
      expand      = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = y_limits
    ) +
    labs( title = "A",
          x = "Month",
          y = expression("Cumulative PM"[2.5]*" (mg/m"^3*")")
    ) +
    theme_classic(base_size = 22) +
    theme(
      legend.position = "top",
      axis.line = element_line(linewidth = line_thin),
      plot.title = element_text(hjust = 0),
      axis.ticks  = element_line(linewidth = line_thin)
    )
}

# Function to generate PM2.5 time-series plot
plot_weekly_avg_pm25_by_month <- function(ind_df, y_limits = NULL, line_thin = 0.25) {
  
  if (is.null(y_limits)) { 
    y_limits <- c(0, max(c(ind_df$weekly_avg_pm25, ind_df$non_wfs_weekly_avg_pm25, ind_df$wfs_weekly_avg_pm25), na.rm = TRUE) * 1.05) 
    }
  
  ggplot(ind_df, aes(x = start_date)) +
    geom_line(aes(y = weekly_avg_pm25, color = "Total"),    linewidth = 1) +
    geom_line(aes(y = non_wfs_weekly_avg_pm25, color = "Non-WFS"), linewidth = 1) +
    geom_line(aes(y = wfs_weekly_avg_pm25, color = "WFS"),     linewidth = 1) +
    scale_color_manual(
      values = c("Total" = "black", "Non-WFS" = "#B0C6DF", "WFS" = "#FF8572"),
      name   = NULL
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b-%y",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = y_limits
    ) +
    labs(title = "B",
         x = "Month",
         y = expression("Mean PM"[2.5]*" ("*mu*"g/m"^3*")")
    ) +
    theme_classic(base_size = 22) +
    theme(
      legend.position = "top",
      axis.line = element_line(linewidth = line_thin),
      plot.title = element_text(hjust = 0),
      axis.ticks  = element_line(linewidth = line_thin)
    )
}

# Function to generate histogram of weekly sums (used to identify counterfactual value)
plot_slope_histogram <- function(ind_df, y_limits = NULL, x_limits = NULL, binwidth = 8, line_thin = 0.25) {
  
  if (is.null(y_limits)) { 
    tmp <- ggplot_build(ggplot(ind_df, aes(x = weekly_sum_pm25)) + geom_histogram(binwidth = binwidth))
    y_max = max(tmp$data[[1]]$count, na.rm = TRUE)
    y_limits <- c(0, y_max * 1.1) 
  }
  if (is.null(x_limits)) {
    x_limits <- c(0, max(ind_df$weekly_sum_pm25))
  }
  
  # Counterfactual value
  median_slope <- median(ind_df[modified_z >= -2 & modified_z <= 2 & n_days == 7,]$weekly_sum_pm25, na.rm = TRUE)
  
  z_vals <- ind_df$modified_z
  z_vals_norm <- which(!is.na(z_vals) & z_vals >= -2 & z_vals <= 2)
  if (length(z_vals_norm) > 0) {
    min_norm_slope <- min(ind_df[z_vals_norm,]$weekly_sum_pm25, na.rm = TRUE)
    max_norm_slope <- max(ind_df[z_vals_norm,]$weekly_sum_pm25, na.rm = TRUE)
  } else {
    min_norm_slope <- max_norm_slope <- NA_real_
  }
  
  norm_slopes_rect <- data.frame(
    xmin = min_norm_slope, xmax = max_norm_slope,
    ymin = 0, ymax = max(y_limits)
  )
  
  ggplot() +
    # Z‐range rectangle
    geom_rect(
      data = norm_slopes_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "-2 < Z-Score < 2"),
      alpha = 0.3, inherit.aes = FALSE
    ) +
    # Histogram bars
    geom_histogram(
      data        = ind_df,
      aes(x = weekly_sum_pm25),
      binwidth    = binwidth,
      boundary    = 0,
      fill        = "grey90",
      color       = "black",
      linewidth   = line_thin,
      show.legend = FALSE
    ) +
    # Cutoff lines (no legend)
    { if (!is.na(min_norm_slope)) geom_vline(xintercept = min_norm_slope, linetype = "dashed", show.legend = FALSE) } +
    { if (!is.na(max_norm_slope)) geom_vline(xintercept = max_norm_slope, linetype = "dashed", show.legend = FALSE) } +
    # Median line
    geom_vline(
      aes(xintercept = median_slope, color = "Median"),
      linewidth = 1.2
    ) +
    # Manual scales
    scale_fill_manual(
      name   = NULL,
      values = c("-2 < Z-Score < 2" = "yellow")
    ) +
    scale_color_manual(
      name   = NULL,
      values = c("Median" = "#0085E4")
    ) +
    guides(
      fill  = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0), oob = scales::oob_keep) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0), oob = scales::oob_keep) +
    labs(
      title = "C",
      x     = expression("Weekly Sum of PM"[2.5]*" ("*mu*"g/m"^3*" per epiweek)"),
      y     = "Count"
    ) +
    theme_classic(base_size = 26) +
    theme(
      plot.title      = element_text(hjust = 0),
      legend.position = "top",
      axis.line       = element_line(linewidth = line_thin),
      panel.grid      = element_blank(),
      axis.ticks      = element_line(linewidth = line_thin)
    )
}

## Use functions to generate figures
ID <- 'Individual23'
cumu_plot <- plot_cumulative_pm25_by_month(all_exp_profile_preg[IndividualID == ID])
ts_plot <- plot_weekly_avg_pm25_by_month(all_exp_profile_preg[IndividualID == ID])
hist_plot <- plot_slope_histogram(all_exp_profile_preg[IndividualID == ID])
title <- ggplot() +
  geom_text(aes(x = 0.5, y = 0.5, label = ID), 
            size = 12, fontface = "bold") + xlim(0, 1) + ylim(0, 1) + theme_void()

title / cumu_plot / ts_plot / hist_plot
