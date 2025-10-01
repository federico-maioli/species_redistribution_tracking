nc_time_to_date <- function(time_vals, time_units, tz = "UTC") {
  # Extract the units (days, hours, minutes, seconds) and origin date
  matches <- regmatches(time_units, regexec("^(\\w+) since (.+)$", time_units))[[1]]
  
  if(length(matches) != 3) stop("Time units not recognized")
  
  unit <- matches[2]           # "days", "hours", "minutes", or "seconds"
  origin <- matches[3]         # reference date-time
  
  # Convert origin to POSIXct
  origin_dt <- as.POSIXct(origin, format = "%Y-%m-%d %H:%M:%S", tz = tz)
  if(is.na(origin_dt)) {        # fallback if no time component
    origin_dt <- as.POSIXct(origin, format = "%Y-%m-%d", tz = tz)
  }
  
  # Convert time values to seconds
  multiplier <- switch(
    tolower(unit),
    "days"    = 86400,
    "hours"   = 3600,
    "minutes" = 60,
    "seconds" = 1,
    stop("Unknown time unit: ", unit)
  )
  
  # Add time values to origin
  origin_dt + time_vals * multiplier
}
