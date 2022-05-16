populate.hypno.reservoir <- function(hID, t,
                                     mean_hypno,
                                     inf_record,
                                     hypno_reservoir) {

  # Identify active infection
  origin_inf <- which(inf_record$infected==hID & is.na(inf_record$end_t))
  # Choose number of hypnozoites
  n_hypnozoites <- rgeom(1, prob=1/(1+mean_hypno))
  # Add new hypnozoites to existing reservoir
  current_hyps <- hypno_reservoir[[hID]]
  if (length(n_hypnozoites) > 0) {
    output <- c(current_hyps, rep(origin_inf, n_hypnozoites))
  } else {
    output <- current_hyps 
  }
  return(output)
}
