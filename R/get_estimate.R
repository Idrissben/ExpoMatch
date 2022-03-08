#This function is aimed at getting the causal estimate of the matches provided (currently supports ATT only)

get_estimate <- function(Tr = Tr, X = X, Y = Y, matches = matches,Expo = TRUE, Expo_weights = NULL) {

  if (Expo){

