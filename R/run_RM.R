#' @export
run.RM <- function(N_h,
                   N_h_t0,
                   N_v,
                   N_v_t0,
                   runtime,
                   bite_rate,
                   raw_hv_trans_rate=NULL,
                   raw_vh_trans_rate=NULL,
                   eff_hv_trans_rate=NULL,
                   eff_vh_trans_rate=NULL,
                   h_rec_rate,
                   v_rec_rate,
                   h_lag,
                   v_lag,
                   h_max_duration,
                   v_max_duration,
                   mean_hypno=0,
                   hyp_act_rate=NULL,
                   hyp_death_rate=NULL,
                   prev_sim_output=NULL) {

  # Check number of seed individuals
  if (N_h < N_h_t0 || N_v < N_v_t0) {
    stop("Number of seed individuals exceeds total number of individuals.")
  }

  # Check hypnozoite parameters
  if (mean_hypno < 0) {
    stop("Mean number of hypnozoites cannot be less than 0")
  } else if (mean_hypno > 0) {
    if (is.null(hyp_act_rate) | is.null(hyp_death_rate)) {
      stop("If implementing hypnozoites, both 'hyp_act_rate' and 'hyp_death_rate' must be defined.")
    }
  }

  # Check that all rates are between 0 and 1
  if (h_rec_rate < 0 | h_rec_rate > 1 | v_rec_rate < 0 | v_rec_rate > 1 | hyp_act_rate < 0 | hyp_act_rate > 1 | hyp_act_rate < 0 | hyp_act_rate > 1) {
    stop("Recovery and hypnozoite rate parameters must be between 0 and 1.")
  }

  # Calculate raw and effective transmission probabilities
  if (is.null(eff_hv_trans_rate) & !is.null(raw_hv_trans_rate)) {
    eff_hv_trans_rate <- calc.trans.rates(lag=h_lag,
                                          max_duration=h_max_duration,
                                          recovery_rate=h_rec_rate,
                                          raw_trans_rate=raw_hv_trans_rate)
  } else if (!is.null(eff_hv_trans_rate) & is.null(raw_hv_trans_rate)) {
    raw_hv_trans_rate <- calc.trans.rates(lag=h_lag,
                                          max_duration=h_max_duration,
                                          recovery_rate=h_rec_rate,
                                          eff_trans_rate=eff_hv_trans_rate)
  } else if (!is.null(eff_hv_trans_rate) & !is.null(raw_hv_trans_rate)) {
    stop("'raw_hv_trans_rate' and 'eff_hv_trans_rate' cannot both be specified.")
  }  else {
    stop("Either 'raw_hv_trans_rate' or 'eff_hv_trans_rate' must be specified.")
  }
  if (is.null(eff_vh_trans_rate) & !is.null(raw_vh_trans_rate)) {
    eff_vh_trans_rate <- calc.trans.rates(lag=v_lag,
                                          max_duration=v_max_duration,
                                          recovery_rate=v_rec_rate,
                                          raw_trans_rate=raw_vh_trans_rate)
  } else if (!is.null(eff_vh_trans_rate) & is.null(raw_vh_trans_rate)) {
    raw_vh_trans_rate <- calc.trans.rates(lag=v_lag,
                                          max_duration=v_max_duration,
                                          recovery_rate=v_rec_rate,
                                          eff_trans_rate=eff_vh_trans_rate)
  } else if (!is.null(eff_vh_trans_rate) & !is.null(raw_vh_trans_rate)) {
    stop("'raw_vh_trans_rate' and 'eff_vh_trans_rate' cannot both be specified.")
  }  else {
    stop("Either 'raw_vh_trans_rate' or 'eff_vh_trans_rate' must be specified.")
  }

  # Calculate equation parameters
  RM_params <- calc.RM.params(N_h=N_h,
                              N_v=N_v,
                              bite_rate=bite_rate,
                              eff_hv_trans_rate=eff_hv_trans_rate,
                              eff_vh_trans_rate=eff_vh_trans_rate,
                              h_rec_rate=h_rec_rate,
                              v_rec_rate=v_rec_rate)
  for (param in names(RM_params)) {
    assign(param, RM_params[[param]])
  }

  # Store parameters
  call <- list(N_h=N_h,
               N_h_t0=N_h_t0,
               N_v=N_v,
               N_v_t0=N_v_t0,
               runtime=runtime,
               bite_rate=bite_rate,
               raw_hv_trans_rate=raw_hv_trans_rate,
               raw_vh_trans_rate=raw_vh_trans_rate,
               eff_hv_trans_rate=eff_hv_trans_rate,
               eff_vh_trans_rate=eff_vh_trans_rate,
               h_rec_rate=h_rec_rate,
               v_rec_rate=v_rec_rate,
               h_max_duration=h_max_duration,
               v_max_duration=v_max_duration,
               h_lag=h_lag,
               v_lag=v_lag,
               prev_sim_output=prev_sim_output)

  model_params <- list(lambda=lambda,
                       phi=phi,
                       gamma=gamma,
                       epsilon=epsilon,
                       H_eq=H_eq,
                       V_eq=V_eq)

  # Create vectors of host and vector names for convenience
  hIDs <- paste0(rep("H", N_h), 1:N_h)
  vIDs <- paste0(rep("V", N_v), 1:N_v)

  # Create matrix to store individual infection states over time
  indiv_status <- as.data.frame(matrix(0, nrow=N_h+N_v, ncol=0))
  rownames(indiv_status) <- c(hIDs, vIDs)
  
  # If specified, create objects to store information about hypnozoite reservoirs
  if (mean_hypno > 0) {
    hypno_reservoir <- list();length(hypno_reservoir) <- N_h;names(hypno_reservoir) <- hIDs
    n_hypno <- as.data.frame(matrix(0, nrow=N_h, ncol=runtime));rownames(n_hypno) <- hIDs;colnames(n_hypno) <- paste0("T",1:runtime)
  }

  # Create dataframe to store infection record
  inf_record <- as.data.frame(matrix(ncol=5, nrow=0))
  colnames(inf_record) <- c("origin_inf","infector", "infected","start_t","end_t")

  ## Run simulation sequentially over each time step
  ## -----------------------------------------------
  for (t in 1:runtime) {
    if (t == 1 & is.null(prev_sim_output)) {

      # -- Step 0 - Initialize w/ infected hosts/vectors
      indiv_status[,t] <- 0 # new time step (adds column to indiv_status)
      colnames(indiv_status)[t] <- paste0("T",t)
      # Choose infected individuals
      new_inf_indivs <- c(sample(hIDs, size=N_h_t0), sample(vIDs, size=N_v_t0))
      indiv_status[new_inf_indivs,t] <- 1 # infect these
      for (type in c("H", "V")) {
        if (length(grep(type, new_inf_indivs)) > 0) {
          entry_rows <- (nrow(inf_record)+1):(nrow(inf_record)+length(grep(type, new_inf_indivs)))
          inf_record[entry_rows,] <- cbind(paste0(c("H","V")[-match(type, c("H","V"))], "-seed"), paste0(c("H","V")[-match(type, c("H","V"))], "-seed"), new_inf_indivs[grep(type, new_inf_indivs)], t, NA)
        }
      }
      # Create hypnozoite reservoir for newly infected hosts
      if (mean_hypno > 0) {
        new_inf_hosts <- new_inf_indivs[new_inf_indivs %in% hIDs]
        if (length(new_inf_hosts) > 0) {
          hypno_reservoir[new_inf_hosts] <- lapply(X=new_inf_hosts, populate.hypno.reservoir, t=t, mean_hypno=mean_hypno, inf_record=inf_record, hypno_reservoir=hypno_reservoir)
        }
      }
    } else {

      # -- Step 1 - All vectors bite hosts, transmitting and/or becoming infected
      indiv_status[,t] <- indiv_status[,t-1] # new time step (adds column to indiv_status)
      colnames(indiv_status)[t] <- paste0("T",t)
      # Simulate vector biting
      dat <- lapply(X=vIDs, FUN=sim.vector.biting,
                    hIDs=hIDs, t=t,
                    bite_rate=bite_rate,
                    vh_trans_rate=raw_vh_trans_rate,
                    hv_trans_rate=raw_hv_trans_rate,
                    h_lag=h_lag,
                    v_lag=v_lag,
                    indiv_status=indiv_status[,paste0("T",(t-1):t)],
                    inf_record=inf_record)
      names(dat) <- vIDs
      # Update individual statuses
      indiv_status[,t] <- summarize.status(dat, t)
      # Record infections
      new_infections <- as.data.frame(do.call(rbind, lapply(dat, function(l) l[[2]])))
      if (nrow(new_infections) > 0) {
        # If an individual is infected more than once in this time step, only keep one
        co_infected <- which(table(new_infections$infected) > 1)
        if (length(co_infected) > 0) {
          co_infected <- names(table(new_infections$infected))[co_infected]
          if (length(co_infected)>1) {msg_plural="s";msg_were="were"} else {msg_plural="";msg_were="was"}
          message(paste0("WARNING: ", length(co_infected), " co-infection", msg_plural, " in time step ", t, " ", msg_were, " simplified"))
          new_infections <- rbind(new_infections[-which(new_infections$infected %in% co_infected),], new_infections[match(co_infected, new_infections$infected),])
        }
        rownames(new_infections) <- (nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections))
        inf_record[(nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections)),] <- new_infections
      }
      # Create hypnozoite reservoir for hosts newly infected by a vector (not re-activations)
      if (mean_hypno > 0) {
        new_inf_hosts <- hIDs[indiv_status[hIDs,t] == 1 & indiv_status[hIDs,t-1] == 0]
        new_inf_hosts <- inf_record[intersect(which(inf_record$infected %in% new_inf_hosts & is.na(inf_record$end_t)), grep("V", inf_record$infector)), "infected"]
        if (length(new_inf_hosts) > 0) {
          hypno_reservoir[new_inf_hosts] <- lapply(X=new_inf_hosts, populate.hypno.reservoir, t=t, mean_hypno=mean_hypno, inf_record=inf_record, hypno_reservoir=hypno_reservoir)
        }
        # Record total number of hypnozoites per host
        n_hypno[,paste0("T",t)] <- sapply(hypno_reservoir, length)
      }
      # -- Step 2 - Dormant hypnozoites activate/die
      if (mean_hypno > 0) {
        # Identify hosts with hypnozoites but no active infection
        dormant_hosts <- rownames(indiv_status[which(indiv_status[hIDs,t]==0 & sapply(hypno_reservoir, FUN=length)>0),])
        if (length(dormant_hosts) > 0) {
          # Simulate hypnozoite death
          hypno_reservoir[hIDs] <- lapply(X=hIDs, sim.hypno.death, hypno_reservoir=hypno_reservoir, hyp_death_rate=hyp_death_rate)
          # Simulate hypnozoite reactivation
          dat <- lapply(X=dormant_hosts, t=t,
                        FUN=sim.hypno.activation,
                        hypno_reservoir=hypno_reservoir,
                        hyp_act_rate=hyp_act_rate,
                        indiv_status=indiv_status[,paste0("T",(t-1):t)])
          names(dat) <- dormant_hosts
          # Update individual statuses
          indiv_status[,t] <- summarize.status(dat, t)
          # Record infections
          new_infections <- as.data.frame(do.call(rbind, lapply(dat, function(l) l[[2]])))
          if (nrow(new_infections) > 0) {
            rownames(new_infections) <- (nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections))
            inf_record[(nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections)),] <- new_infections
          }
          # Update hypnozoite reservoir
          hypno_reservoir[dormant_hosts] <- lapply(dat, function(l) l[[3]])
        }
      }

      # -- Step 3 - Infected individuals clear infections
      for (type in c("h", "v")) {
        # Determine who is infected
        ids <- get(paste0(type, "IDs"))
        inf_indivs <- ids[which(indiv_status[ids,t-1]==1)]
        if (length(inf_indivs) > 0) {
          # Clear infections with 'recovery rate' probability
          indiv_status[inf_indivs,t] <- indiv_status[inf_indivs,t-1] - rbinom(n=length(inf_indivs), size=1, prob=get(paste0(type,"_rec_rate")))
          # Clear infections from individuals that have exceeded their maximum duration
          active_infs <- inf_record[inf_record$infected %in% inf_indivs & is.na(inf_record$end_t),]
          too_long <- active_infs[t - as.numeric(active_infs$start_t) == get(paste0(type,"_max_duration"))-1, "infected"]
          if (length(too_long) > 0) {
            type.full <- c("host", "vector")[match(type, c("h", "v"))]
            if (length(too_long)>1) {msg_plural="s";msg_were="were"} else {msg_plural="";msg_were="was"}
            message(paste0("WARNING: ", length(too_long), " ", type.full, " infection", msg_plural, " ", msg_were, " forcibly terminated in time step ", t))
            indiv_status[too_long,t] <- 0
          }
          # Update infection record
          cleared <- inf_indivs[indiv_status[inf_indivs,t]==0]
          inf_record[inf_record$infected %in% cleared & is.na(inf_record$end_t),"end_t"] <- t
        }
      }
    }
  }

  # Fix, store, and return output
  inf_record$start_t <- as.numeric(inf_record$start_t)
  inf_record$end_t <- as.numeric(inf_record$end_t)
  inf_record$duration <- inf_record$end_t - inf_record$start_t
  inf_record$inf_id <- rownames(inf_record)
  inf_record <- inf_record[,c("inf_id", "origin_inf","infector", "infected","start_t","end_t")]
  if (mean_hypno > 0) {
    output <- list(call=call, model_parameters=model_params, indiv_status=indiv_status, infection_record=inf_record, hypno_reservoir=hypno_reservoir, n_hypno=n_hypno)
  } else {
    output <- list(call=call, model_parameters=model_params, indiv_status=indiv_status, infection_record=inf_record)
  }
  return(output)
}
