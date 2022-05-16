summarize.status <- function(dat, t) {
  dat0 <- lapply(dat, function(l) l[[1]][,paste0("T",t)])
  new_status <- apply(X=do.call(cbind,dat0), FUN=sum, MARGIN=1)
  if (any(new_status>1)) {
    new_status[new_status>1] <- 1
  }
  return(new_status)
}