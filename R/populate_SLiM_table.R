populate.SLiM.table <- function(pruned_inf_record, inf_id, RM_sample_time) {
  # Filter to infections that resulted from inf_id
  dat <- pruned_inf_record[which(pruned_inf_record$origin_inf==inf_id),]
  # Determine infector ID
  infector <- unique(dat$infector)
  # Determine who was infected over the course of this infection
  infected <- dat$infected
  # Determine the infector type
  infector_type <- substring(infector,1,1)
  # Determine when these infections began in Ross-Macdonald sim time
  RM_start_times <- dat$start_t
  # If not a seed infection...
  if (length(grep("seed", inf_id))==0) {
    # Are any of the resulting infections sampled?
    is_sample <- any(dat[dat$inf_id, "end_t"] == "sample")
    if (!is.na(is_sample) & is_sample==T) {
      infected <- c(infected, paste0(infector_type, "_sample"))
      RM_start_times <- c(RM_start_times, RM_sample_time)
    }
    # Determine origin of inf_id
    origin <- pruned_inf_record[pruned_inf_record$inf_id==inf_id,"origin_inf"]
    # Determine infector id
    infector_id <- substring(infector, 2)
    # Determine when the infector was infected in RM time
    infector_start_time <- pruned_inf_record[pruned_inf_record$inf_id==inf_id,"start_t"]
    # If a seed infection...
  } else {
    origin <- NA
    infector_id <- "seed"
    infector_start_time <- NA
  }
  # Determine when the infections began wrt to infector's status
  start_times <- suppressWarnings(RM_start_times - infector_start_time)
  # Return output
  output <- c(inf_id,
              origin,
              infector_type,
              infector_id,
              paste(infected, collapse=";"),
              paste(RM_start_times, collapse=";"),
              paste(start_times, collapse=";"))
  return(output)
}
