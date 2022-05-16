sim.transmission <- function(h, v, t,
                             v_inf_status,
                             v_lag,
                             h_lag,
                             hv_trans_rate,
                             vh_trans_rate,
                             indiv_status,
                             inf_record) {

  infections <- NULL
  infection_occurs <- 0
  # Check host infection status
  h_inf_status <- indiv_status[h,paste0("T",t-1)]

  # If vector is infected and host is susceptible...
  if (v_inf_status == 1 & h_inf_status == 0) {
    # Does the host become infected?
    inf_age <- t - as.numeric(inf_record[inf_record$infected==v & is.na(inf_record$end_t), "start_t"])+1
    if (inf_age > v_lag) {
      infection_occurs <- rbinom(n=1, size=1, prob=vh_trans_rate)
    }
    if (infection_occurs == 1) {
      # Update individual infection status
      indiv_status[h, paste0("T",t)] <- 1
      # Determine index/source of the vector's current infection
      origin_inf <- rownames(inf_record[inf_record$infected==v & is.na(inf_record$end_t),])
      # Record transmission event
      infections <- c(origin_inf, v, h, t, NA)
    }
  }

  # If host is infected and vector is susceptible...
  if (v_inf_status == 0 & h_inf_status == 1) {
    # Does the vector become infected?
    inf_age <- t - as.numeric(inf_record[inf_record$infected==h & is.na(inf_record$end_t), "start_t"])+1
    if (inf_age > h_lag) {
      infection_occurs <- rbinom(n=1, size=1, prob=hv_trans_rate)
    }
    if (infection_occurs == 1) {
      # Update individual status
      indiv_status[v, paste0("T",t)] <- 1
      # Determine index of the host's current infection
      origin_inf <- rownames(inf_record[inf_record$infected==h & is.na(inf_record$end_t),])
      # Record transmission event
      infections <- c(origin_inf, h, v, t, NA)
    }
  }

  return(list(indiv_status=indiv_status, new.infections=infections))
}
