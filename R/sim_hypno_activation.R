sim.hypno.activation <- function(h, t,
                                 hypno_reservoir,
                                 hyp_act_rate,
                                 indiv_status) {

  infections <- NULL
  # Decide which hypnozoites activate
  hypnozoites <- hypno_reservoir[[h]]
  activated_hypnozoite <- hypnozoites[rbinom(length(hypnozoites), 1, hyp_act_rate) == 1]
  # If more than one activates, choose just one
  if (length(activated_hypnozoite) > 1) {
    message(paste0("WARNING: ", length(activated_hypnozoite)," hypnozoites were activated in time step ", t, " - only 1 was kept"))
    activated_hypnozoite <- sample(activated_hypnozoite, size=1)
  }

  # Activate hypnozoites
  if (length(activated_hypnozoite) == 1) {
    # Remove hypnozoite from reservoir
    hypnozoites <- hypnozoites[-match(activated_hypnozoite, hypnozoites)]
    # Update host infection status
    indiv_status[h, paste0("T", t)] <- 1
    # Record transmission event
    infections <- c(activated_hypnozoite, h, h, t, NA)
  }
  return(list(indiv_status=indiv_status, new.infections=infections, hypnozoites=hypnozoites))
}
