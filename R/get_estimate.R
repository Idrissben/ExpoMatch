get_estimate <- function(Tr = Tr, X = X, Y = Y, matches = matches,Expo = TRUE, Expo_weights = NULL) {
  
  sum_treated = sum(Y[matches[,1]]*matches[,3])
  sum_control = sum(Y[matches[,2]]*matches[,3])
  
  return((sum_treated - sum_control) / sum(matches[,3]))
}
