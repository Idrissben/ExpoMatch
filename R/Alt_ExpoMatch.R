#' Expands Genetic Matching
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export

# Balance matrix should NOT contain the outcome or treatment vector

# fit function #1 is Chris's biggest linear combination discrepancy metric


Alt_ExpoMatch <- function(Tr,
                          X,
                          pop.size = 10,
                          max.generations = 5,
                          domains = c(0.50, 1.99),
                          start.gm = TRUE,
                          start.weights = NULL,
                          print.level = 1) {

  start.time <- Sys.time()

  # series of QoL checks
  if (length(Tr) != nrow(X)) {
    return('ERROR - INCOMPATIBLE LENGTHS')
  }
  if (min(X) < 0) {
    return('CAREFUL - NEGATIVE INPUT VALUES')
  }
  if (is.na(max(X)**domains[2])) {
    return('HIGHEST POSSIBLE EXPONENT OUT OF BOUNDS')
  }
  
  if(!"Matching" %in% (.packages())){
        print("Missing Matching library, can not proceed")

  }
  if(!"rgenoud" %in% (.packages())){
    print("Missing rgenoud library, can not proceed")
  }
  # check if user inputted starting weights but didn't disable starting genmatch
  if (!is.null(start.weights) & start.gm) {
    print('Deactivating starting genmatch, because you inputted starting weights.')
    start.gm <- FALSE
  }

  # store number of variables and number of observations
  n.obs <- length(Tr)
  n.var <- ncol(X)

  #Creating the balance Matrix
  BM <- scale(X)

  
  #Creating the custom leximin imbalance function
  myfit.imb <- function(matches, BM) {
  # ASSUME
  # no outcome and no treatment column in BM
  # BM is scaled to mean=0, sd=1

  # prepare inputs for Chris's metric
  trt.X <- BM[matches[, 1], ]
  ctrl.X <- BM[matches[, 2], ]
  
  # find the difference between treated and control units
  # this returns a matrix with element-wise differences
  t <- trt.X - ctrl.X
  # sum each column, and square the value & then sum each element
  return(sum(colSums(t)^2))
  
  }

  if (start.gm) {
    # run an initial GenMatch to get starting parameters
    genout.start <- GenMatch(Tr = Tr,
                             X = X,
                             pop.size = pop.size*50,
                             max.generations = max.generations*2,
                             BalanceMatrix = BM,
                             hard.generation.limit = TRUE,
                             fit.func = myfit.imb)
    print("The initial GenMatch run determined the following variable weights:")
    print(genout.start$par)
    start.weights <- genout.start$par
    n.p.values <- length(genout.start$values)
  }

  if (is.null(start.weights)) {
    start.weights <- rep(1, n.var)
  }

  # prepare domains
  dom <- cbind(rep(domains[1], n.var), rep(domains[2], n.var))
  
  

  
  
  GenMatchWrapper <- function(exponents) {

    if (print.level == 2) {
      print(exponents)
    }

    XN <- X
    
    
    for (i in c(1:n.var)) {
      XN[, i] <- XN[, i]^exponents[i]
    }

    genout <- GenMatch(Tr = Tr,
                       X = XN,
                       print.level = 1,
                       project.path = paste(tempdir(), "/genoud.txt", sep = ""),
                       pop.size = pop.size,
                       max.generations = max.generations,
                       BalanceMatrix = BM,
                       starting.values = start.weights,
                       fit.func = myfit.imb)

    return(genout$value[1]) # = highest lowest p-value
  }

  genoudout <- genoud(GenMatchWrapper,
                      nvars = n.var,
                      max = FALSE,
                      pop.size = pop.size,
                      max.generations = max.generations,
                      Domains = dom,
                      boundary.enforcement = 2)

  # parse the file to find the best result's weights
  file_data <- read.delim(paste(tempdir(), "/genoud.txt", sep = ""), skip = 1, header = FALSE, nrows = 1)
  best.weights <- file_data[1, (2+n.p.values):(1+n.p.values+n.var)]
  best.weights <- as.numeric(best.weights[1,])

  XM <- X
  for (i in c(1:n.var)) {
    XM[, i] <- XM[, i]^genoudout$par[i]
  }

  if (print.level == 2) {
    print(head(XM))
  }

  genout.fin <- GenMatch(Tr = Tr,
                         X = XM,
                         pop.size = pop.size*50,
                         max.generations = max.generations*2,
                         BalanceMatrix = BM,
                         starting.values = best.weights,
                         fit.func = myfit.imb)
  mout.fin <- Match(Tr = Tr, X = XM, Weight.matrix = genout.fin)

  # prepare output
  end.time <- Sys.time()
  outputlist <- list(mout = mout.fin,
                     genout = genout.fin,
                     matches = genout.fin$matches,
                     pvalues = genout.fin$value,
                     weights = genout.fin$par,
                     exponents = genoudout$par,
                     time = end.time - start.time)
  if (start.gm) {
    outputlist[['start.weights']] = genout.start$par
    outputlist[['start.pvalues']] = genout.start$value
  }

  print('')
  print("##########################")
  print('')

  if (print.level == 2) {
    print(genout.fin$matches)
  }
  print(sprintf("p-values: %s", genout.fin$value))
  print(sprintf("weights: %s", genout.fin$par))
  print(sprintf("exponents: %s", genoudout$par))
  print(sprintf("time taken: %s", end.time - start.time))
  return(outputlist)
}
