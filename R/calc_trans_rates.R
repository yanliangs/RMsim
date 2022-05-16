#' @export
calc.trans.rates <- function(lag, 
                             max_duration,
                             recovery_rate,
                             raw_trans_rate=NULL,
                             eff_trans_rate=NULL) {
  if (is.null(raw_trans_rate) & is.null(eff_trans_rate)) {
    stop("Either 'raw_trans_rate' or 'eff_trans_rate' must be specified.")
  }
  if (!is.null(raw_trans_rate) & is.null(!eff_trans_rate)) {
    stop("'raw_trans_rate' and 'eff_trans_rate' cannot both be specified.")
  }

  trans_prob_vector <- c(rep(0, lag), rep(1, max_duration-lag))
  eff_trans_distr <- trans_prob_vector * dgeom(0:(max_duration-1), recovery_rate)
  if (!is.null(raw_trans_rate)) {
    eff_trans_rate <- raw_trans_rate * MESS::auc(1:max_duration,
                                                 eff_trans_distr,
                                                 type="spline")
    if (eff_trans_rate > 1) {
      stop(paste("Calculated effective transmission probability is", signif(eff_trans_rate, digits=3),"- cannot be greater than 1."))
    } else {
      return(eff_trans_rate)
    }
  }
  if (!is.null(eff_trans_rate)) {
    raw_trans_rate <- eff_trans_rate/MESS::auc(1:max_duration,
                                               eff_trans_distr,
                                               type="spline")
    if (raw_trans_rate > 1) {
      stop(paste("Calculated  transmission probability is", signif(raw_trans_rate, digits=3),"- cannot be greater than 1."))
    } else {
      return(raw_trans_rate)
    }
  }
}
