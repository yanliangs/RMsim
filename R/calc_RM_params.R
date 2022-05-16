#' @export
calc.RM.params <- function(N_h,
                           N_v,
                           bite_rate,
                           eff_hv_trans_rate,
                           eff_vh_trans_rate,
                           h_rec_rate,
                           v_rec_rate) {
  phi = bite_rate * eff_hv_trans_rate
  lambda = bite_rate * (N_v/N_h) * eff_vh_trans_rate
  gamma = h_rec_rate
  epsilon = v_rec_rate
  H_eq = (lambda*phi - gamma*epsilon)/(lambda*phi + phi*gamma)
  V_eq = (lambda*phi - gamma*epsilon)/(lambda*phi + lambda*epsilon)
  if (H_eq < 0) {
    stop("Calculated host equilibrium is < 0. Please choose different starting parameters.")
  } else if (V_eq < 0) {
    stop("Calculated vector equilibrium is < 0. Please choose different starting parameters.")
  } else {
    return(list(phi=phi,
                lambda=lambda,
                gamma=gamma,
                epsilon=epsilon,
                H_eq=H_eq,
                V_eq=V_eq))
  }
}