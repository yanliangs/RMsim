#' @export
prune.infection.record <- function(inf_record, focal_infections) {
  `%ni%` <- Negate(`%in%`)
  # Determine seed infections
  seed_infs <- inf_record[grep("seed", inf_record$infector),"inf_id"]
  # Keep all focal transmissions
  keep_infs <- c()
  keep_infs[(length(keep_infs)+1):(length(keep_infs)+length(focal_infections))] <- focal_infections
  # Loop over each focal infection and trace back until a seed infection is reached
  for (i in focal_infections) {
    j <- i
    while (j %ni% seed_infs) {
      # update j
      j <- inf_record[inf_record$inf_id==j, "origin_inf"]
      # keep this infection in pruned record
      keep_infs[length(keep_infs)+1] <- j
    }
  }
  # Sort and remove duplicates
  keep_infs <- sort(unique(as.numeric(keep_infs)))
  pruned_inf_record <- inf_record[inf_record$inf_id %in% keep_infs,]
  pruned_inf_record[pruned_inf_record$inf_id %in% focal_infections, "end_t"] <- "sample"
  return(pruned_inf_record)
}
