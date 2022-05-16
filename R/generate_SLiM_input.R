#' @export
generate.SLiM.input <- function(pruned_inf_record, RM_sample_time, sim_params) {
  output <- as.data.frame(t(sapply(X=unique(pruned_inf_record$origin_inf),
                                   FUN=populate.SLiM.table,
                                   pruned_inf_record=pruned_inf_record,
                                   RM_sample_time=RM_sample_time)))
  colnames(output) <- c("infection_idx", "infection_source", "infector_type", "infector_id", "infected_ids", "RM_time_start", "inf_time_start")
  rownames(output) <- 1:nrow(output)
  # fill in infection start times
  for (infector_type in c("h", "v")) {
    seed_infs <- intersect(grep("NA", output$inf_time_start), which(output$infector_type==toupper(infector_type)))
    infected_type <- c("v", "h")[-match(infector_type, c("v", "h"))]
    trans_prob_distr <- (1-sim_params[[paste0(infector_type, "_rec_rate")]])^(1:(sim_params[[paste0(infector_type, "_max_duration")]])) * sim_params[["bite_rate"]] * sim_params[[paste0("raw_", infector_type, infected_type, "_trans_rate")]]
    trans_prob_distr[1:(sim_params[[paste0(infector_type, "_lag")]]-1)] <- 0
    trans_prob_distr <- as.data.frame(cbind(1:sim_params[[paste0(infector_type, "_max_duration")]], trans_prob_distr));colnames(trans_prob_distr) <- c("x", "prob")
    trans_prob_distr <- pdqr::new_r(trans_prob_distr, type="discrete")
    for (i in seed_infs) {
      start_times <- unlist(strsplit(output[i, "inf_time_start"], split=";"))
      output[i, "inf_time_start"] <- paste(trans_prob_distr(length(start_times)), collapse=";")
    }
  }
  return(output)
}
