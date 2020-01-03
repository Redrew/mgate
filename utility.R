# utility functions
find_and_empty <- function(v, target=NA) {
  if (is.na(target)) v[is.na(v)] <- rep("", sum(is.na(v))) 
  else v[v==target] <- rep("", sum(v==target)) 
  return(v)
}

cappend <- function(v, x) {
  if (is.null(v)) return(x)
  if (is.na(v)) return(x)
  return(c(v, x))
}

make_df <- function(column_names) {
  df <- data.frame(matrix(ncol=length(column_names), nrow=0))
  colnames(df) <- column_names
  return(df)
}

exclude <- function(x, ...) {
  return(x[!x %in% c(...)])
}

in_range <- function(x, range, sd, sig=.975) {
  top <- max(range) + qnorm(sig, sd=sd)
  bottom <- min(range) - qnorm(sig, sd=sd)
  return(x >= bottom & x <= top)
}

zscore <- function(sig, sd=1) {
  return(qnorm(0.5 + sig/2, sd=sd))
}
