#getting weather data using daymet
#libraries 
library(easypackages) # Assists on loading packages and installing them if they are not
packages('tidyverse', 'dplyr', 'tidyr') # Data wrangling
packages('lubridate') # Dates operations
packages('kableExtra') # Table formatting
packages('daymetr','chirps','nasapower') # Weather databases
packages('vegan') # Shannon Diversity Index
packages('skimr') # For checking weather data after sourcing
packages('fuzzyjoin')
packages('ggplot2')

data <- read.csv("daymet.csv")
head(data)
str(data)

# Convert Start and End columns to Date format with the specified format
data$Start <- as.Date(data$Start, format = "%m/%d/%Y")
data$End <- as.Date(data$End, format = "%m/%d/%Y")
str(data)
#using daymet to get the weather data
#view data
kable(data) %>% 
  kable_styling(latex_options = c("striped"), position = "center", font_size = 10)

# Constants for ET0 (Cobaner et al., 2017)
# Solar constant:
Gsc <- 0.0820 # (MJ m-2 min-1)
# Radiation adjustment coefficient (Samani, 2004)
kRs <- 0.17

## a. DAYMET
weather.daymet <- function(data){ 
  # Downloads the daily weather data from the DAYMET database and process it
  # Args:
  #  input = input file containing the locations and the start & end dates for the time series
  #  dpp = days prior to the Start
  # Returns:
  #  a tibble of DAYMET weather variables for the requested time period
  data %>%
    dplyr::mutate(
      Weather = purrr::pmap(list(ID = location,
                                 lat = latitude,
                                 lon = longitude,
                                 sta = Start,
                                 end = End),
                            
                            # Retrieving daymet data:
                            function(ID, lat, lon, sta, end) {
                              daymetr::download_daymet(site = ID,
                                                       lat = lat, 
                                                       lon = lon,
                                                       # Extracting year:
                                                       start = as.numeric(substr(sta, 1, 4)),
                                                       end = as.numeric(substr(end, 1, 4)),
                                                       internal = TRUE, 
                                                       simplify = TRUE)})) %>% 
    # Organizing dataframe (Re-arranging rows and columns)
    dplyr::mutate(Weather = Weather %>% 
                    # Adjusting dates format:
                    purrr::map(~ 
                                 dplyr::mutate(., 
                                               Date = as.Date(as.numeric(yday) - 1, # Day of the year
                                                              origin = paste0(year, '-01-01')),
                                               Year = year(Date),
                                               Month = month(Date),
                                               Day = mday(Date))) %>% 
                    purrr::map(~ 
                                 dplyr::select(., yday, Year, Month, Day,
                                               Date, measurement, value)) %>%
                    # NOTE: spread has been superseded. Use pivot_wider:
                    # purrr::map(~ spread(., 'measurement', 'value'))  %>% 
                    purrr::map(~ 
                                 tidyr::pivot_wider(.,
                                                    names_from = measurement, values_from = value)) %>%
                    # Renaming variables:
                    # NOTE: rename_all() has been superseded by rename_with()
                    purrr::map(~ rename_with(., ~c(
                      "DOY",   # Date as Day of the year
                      "Year",  # Year
                      "Month", # Month 
                      "Day",   # Day of the month
                      "Date",  # Date as normal format
                      "DL",    # Day length (sec)
                      "PP",    # Precipitation (mm)
                      "Rad",   # Radiation (W/m2)
                      "SWE",   # Snow water (kg/m2)
                      "Tmax",  # Max. temp. (degC)
                      "Tmin",  # Min. temp. (degC)
                      "VPD")))) %>%   # Vap Pres Def (Pa)
    # Processing data given start and ending dates:
    dplyr::mutate(Weather = purrr::pmap(list(sta = Start, 
                                             end = End, 
                                             data = Weather), # Requested period
                                        function(sta, end, data) {
                                          dplyr::filter(data, Date >= sta & Date <= end) 
                                        })) %>% 
    tidyr::unnest(cols = c(Weather)) %>% 
    
    # Converting units or adding variables:
    dplyr::mutate(Rad = Rad*0.000001*DL, # Radiation (W/m2 to MJ/m2)
                  Tmean = (Tmax+Tmin)/2, # Mean temperature (degC),
                  VPD = VPD / 1000, # VPD (Pa to kPa),
                  # Creating variables for ET0 estimation:
                  lat_rad = latitude*0.0174533,
                  dr = 1 + 0.033*cos((2*pi/365)*DOY),
                  Sd = 0.409*sin((2*pi/365)*DOY - 1.39),
                  ws = acos(-tan(lat_rad)*tan(Sd)),
                  Ra = (24*60)/(pi) * Gsc * dr * (ws*sin(lat_rad)*sin(Sd) + cos(lat_rad)*sin(ws)),
                  ET0_HS = 0.0135 * kRs * (Ra / 2.45) * (sqrt(Tmax-Tmin)) * (Tmean + 17.8),
                  # Extreme PP events
                  EPE_i = case_when((PP > 25) ~ 1, TRUE ~ 0),
                  # Extreme Temp events
                  ETE_i = case_when((Tmax >= 30) ~ 1, TRUE ~ 0),
                  # Day length (hours)
                  DL = (DL/60)/60, 
                  #diurnal temp
                  DT = Tmax-Tmin,
                  # Thermal time (Growing Degree Days)
                  TT = ifelse(((Tmax + Tmin) / 2) - 5 < 0, 0, ((Tmax + Tmin) / 2) - 5)) %>% 
    dplyr::select(-lat_rad, -dr, -Sd, -ws, -Ra)
}

df.weather.daymet <- weather.daymet(data = data)
# Overview of the variables (useful checking for missing values):
skimr::skim(df.weather.daymet)

# Exporting data as a .csv file
write.csv(df.weather.daymet, row.names = FALSE, na = '', file = paste0('output_daymet_acr.csv'))

str(df.weather.daymet)
summary(df.weather.daymet)

#calculating mean of weather data for locations
library(dplyr)
library(fuzzyjoin)

# Load weather data and ensure 'Date' is in the correct format
weather <- read.csv("output_daymet_acr.csv")
weather$Date <- as.Date(weather$Date, format = "%Y-%m-%d")  # Adjust the format as per your data
head(weather)

# Load agronomic data and convert 'sd' and 'hd' into Date format
data <- read.csv("acr_ag.csv")
data <- data %>%
  mutate(
    sd = as.Date(sd, format = "%m/%d/%Y"),  # Adjust the format as per your data
    hd = as.Date(hd, format = "%m/%d/%Y")   # Adjust the format as per your data
  )
head(data)

# Perform a fuzzy join based on location and date range conditions
merged_data <- fuzzy_left_join(
  data, weather,
  by = c("location" = "location"),
  match_fun = list(`==`)  # Matches exactly on location
) %>%
  filter(Date >= sd & Date <= hd) %>%  # Ensures weather data falls between seeding and harvesting dates
  select(-location.y)  # Removes the redundant 'location.y' column

# Rename location.x to just 'location' for clarity
merged_data <- rename(merged_data, location = location.x)

# View the structure and summary of the cleaned merged data
summary(merged_data)
str(merged_data)

#######
weather_summary <- merged_data %>%
  # Filter rows to include only dates between sd and hd for each location
  filter(Date >= sd & Date <= hd) %>%
  group_by(location, sd, hd, genotype, id, block) %>%
  summarize(
    mean_DL = mean(DL, na.rm = TRUE),
    mean_PP = mean(PP, na.rm = TRUE),
    mean_Rad = mean(Rad, na.rm = TRUE),
    mean_SWE = mean(SWE, na.rm = TRUE),
    mean_Tmax = mean(Tmax, na.rm = TRUE),
    mean_Tmin = mean(Tmin, na.rm = TRUE),
    mean_VPD = mean(VPD, na.rm = TRUE),    mean_Tmean = mean(Tmean, na.rm = TRUE),
    mean_ET0_HS = mean(ET0_HS, na.rm = TRUE),
    mean_DT =mean (DT, na.rm = TRUE),
    mean_TT = mean(TT,na.rm= TRUE),
    # Sum of GDD (TT) only within the sd to hd period
    total_GDD = sum(TT, na.rm = TRUE)
  ) %>%
  ungroup()


str(weather_summary)
write.csv(weather_summary, row.names = FALSE, na = '', file = paste0('weather_acr.csv'))



# Ensure merged_data is loaded and date types are correct
# Compute summary statistics for the grain filling period
weather_summary_grain_filling <- merged_data %>%
  # Filter to include only dates between heading_jd (start of grain filling) and maturity (end of grain filling)
  filter(DOY >= heading_jd & DOY <= maturity) %>%
  # Further filter to include only days with Tmean >= 10°C
  filter(Tmean >= 10) %>%
  group_by(location, genotype, id, block) %>%
  summarize(
    mean_DL = mean(DL, na.rm = TRUE),
    mean_PP = mean(PP, na.rm = TRUE),
    mean_Rad = mean(Rad, na.rm = TRUE),
    mean_SWE = mean(SWE, na.rm = TRUE),
    mean_Tmax = mean(Tmax, na.rm = TRUE),
    mean_Tmin = mean(Tmin, na.rm = TRUE),
    mean_VPD = mean(VPD, na.rm = TRUE),
    mean_Tmean = mean(Tmean, na.rm = TRUE),
    mean_ET0_HS = mean(ET0_HS, na.rm = TRUE),
    mean_DT = mean(DT, na.rm = TRUE),
    mean_TT = mean(TT, na.rm = TRUE),
    total_GDD = sum(TT, na.rm = TRUE),  # Sum the daily Growing Degree Days (TT)
    total_PP = sum(PP, na.rm = TRUE),   # Total precipitation (Σ x)
    total_Tmean = sum(Tmean, na.rm = TRUE),  # Total temperature (Σ t)
    HTC = (total_PP / total_Tmean) * 10,  # Calculate Hydrothermal Coefficient
    .groups = "drop"  # Remove automatic grouping to prevent unwanted summary printing
  )

# View the resulting dataset
print(weather_summary_grain_filling)


# Optionally, check the structure of the output to verify everything is as expected
str(weather_summary_grain_filling)

write.csv(weather_summary_grain_filling, row.names = FALSE, na = '', file = paste0('weather_summary_grain_filling.csv'))

#merging data
# Load the agronomic data from 'acr_all.csv'
agr_data <- read.csv("acr_ag.csv")

# Load the weather summary data for the grain filling period
grain_filling_weather <- read.csv("weather_summary_grain_filling.csv")# Convert dates in agr_data if necessary, assuming 'sd' and 'hd' are date fields

agr_data <- agr_data %>%
  mutate(
    sd = as.Date(sd, format = "%m/%d/%Y"),
    hd = as.Date(hd, format = "%m/%d/%Y")
  )

#check the column names for consistency
print(names(agr_data))
print(names(grain_filling_weather))

# Merge the data based on common columns
merged_data <- left_join(agr_data, grain_filling_weather, by = c("location", "genotype", "id", "block"))

# Check the structure of the merged data to ensure the join was successful
str(merged_data)
# Check for missing values in key columns
summary(merged_data)

# Write the merged data to a new CSV file
write.csv(merged_data, "acr_all.csv", row.names = FALSE)



