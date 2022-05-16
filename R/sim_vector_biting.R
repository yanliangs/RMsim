sim.vector.biting <- function(v, hIDs, t,
                              bite_rate,
                              v_lag,
                              h_lag,
                              hv_trans_rate,
                              vh_trans_rate,
                              indiv_status,
                              inf_record) {

  # Determine vector's infection status
  v_inf_status <- indiv_status[v, paste0("T",t-1)]
  # Determine number of bites
  n_bites <- rpois(n=1, lambda=bite_rate)
  if (n_bites > 0) {
    # Choose hosts to bite
    hosts_bitten <- sample(x=hIDs, size=n_bites, replace=T)
    # Bite hosts and transmit or become infected
    dat_v <- lapply(X=hosts_bitten, FUN=sim.transmission,
                    v=v, t=t,
                    v_inf_status=v_inf_status,
                    v_lag=v_lag,
                    h_lag=h_lag,
                    hv_trans_rate=hv_trans_rate,
                    vh_trans_rate=vh_trans_rate,
                    indiv_status=indiv_status,
                    inf_record=inf_record)
    names(dat_v) <- hosts_bitten
    # Update individual status
    indiv_status[,paste0("T",t)] <- summarize.status(dat_v, t)
    # Record new infections
    new_infections <- as.data.frame(do.call(rbind,lapply(dat_v, function(l) l[[2]])))
    if (ncol(new_infections)==0) {
      new_infections <- as.data.frame(matrix(ncol=5, nrow=0))
    }
  } else {
    new_infections <- as.data.frame(matrix(ncol=5, nrow=0))
  }
  colnames(new_infections) <- c("origin_inf","infector", "infected","start_t","end_t")
  return(list(indiv_status=indiv_status, new_infections=new_infections))
}
