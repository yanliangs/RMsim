sim.hypno.death <- function(h,
                            hypno_reservoir,
                            hyp_death_rate) {

  # Determine while hypnozoites die
  hypnozoites <- hypno_reservoir[[h]]
  hyp_deaths <- rbinom(n=length(hypnozoites), size=1, prob=hyp_death_rate);sum(hyp_deaths)
  # Remove these hypnozoites from the reservoir
  if (sum(hyp_deaths) > 0) {
    hypnozoites <- hypnozoites[-which(hyp_deaths==1)]
  }
  return(hypnozoites)
}
