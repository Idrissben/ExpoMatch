############SET-UP
df = data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("name", "time taken", "Supremum value",
                                                             "Best p-value","Metric used",'Ground truth ATE',
                                                             'Estimated ATE','Standard Error','pop.size','seed'))))

drive_deauth()
names = drive_ls("https://drive.google.com/drive/u/0/folders/1j60vqxqes-vs_I094sIZgQIllv_arMFT")
library(googledrive)
### MONTE CARLO SIMULATION
################################### EXPOMATCH

library(Matching)
library(rgenoud)

Alt_ExpoMatch <- function(Tr,
                          X,
                          Y,
                          pop.size = 10,
                          max.generations = 5,
                          wait.generations=5,
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
    
    t <- t*matches[,3]
    # sum each column, and square the value & then sum each element
    return(sum(colSums(t)^2))
    
  }
  
  if (start.gm) {
    # run an initial GenMatch to get starting parameters
    genout.start <- GenMatch(Tr = Tr,
                             X = X,
                             pop.size = pop.size*10,
                             max.generations = max.generations*2,
                             BalanceMatrix = BM,
                             wait.generations = 5,
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
                      pop.size = pop.size/5,
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
  
  mout.fin <- Match(Tr = Tr, X = XM,Y=Y, Weight.matrix = genout.fin)
  
  # prepare output
  end.time <- Sys.time()
  outputlist <- list("mout" = mout.fin,
                     "genout" = genout.fin,
                     "matches" = genout.fin$matches,
                     "fitness" = genout.fin$value,
                     "weights" = genout.fin$par,
                     "exponents" = genoudout$par,
                     "time" = end.time - start.time)
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
  
  t <- t*matches[,3]
  # sum each column, and square the value & then sum each element
  return(sum(colSums(t)^2))
}

#Extracting TE
TE_finder <- function(s) as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", s))

#main function

maccabee <- function(dtaname,id,seed=123,pop.size=200,fit.func = myfit.imb,wait.generations){
  dta = try(read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", toString(id)),skipNul = TRUE,check.names=FALSE))
  if(inherits(dta, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    return(numeric(length = 10))
  }
  print('data imported successfully')
  dta$X0 = dta$X0 +1 
  dta$X1 = dta$X1 +1
  dta$X2 = dta$X2+1
  dta = dta[,2:6]
  print('data rescaled successfully')
  
  d.X <- cbind(dta$X0,dta$X1,dta$X2,dta$X0^2,dta$X1^2,dta$X2^2)
  colnames(d.X) <- c('X0','X1','X2','X0^2','X1^2','X2^2')
  d.Tr <- dta$T
  d.Y <- dta$Y
  d.form <- T ~ X0 + X1 + X2 + I(X0^2) + I(X1^2) + I(X2^2)
  new.form <- ~ X0 + X1 + X2 + I(X0^2) + I(X1^2) + I(X2^2)
  
  
  # start timer
  start.time <- Sys.time()
  
  # check google sheets authorization upload if needed
  #if (!is.null(out.url)) {
  # gs4_auth()
  #}
  
  # create balance matrix (mean=0, sd=1)
  BM = scale(d.X)
  n.vars = ncol(d.X)
  new.BM <- scale(d.X)
  
  
  
  # generate weights (always include all 1s)
  set.seed(seed)
  
  metric = "Chris's metric"
    
  gm_out = Alt_ExpoMatch(Tr=d.Tr,X=d.X,Y=d.Y,pop.size=pop.size,max.generations = 10,wait.generations=5)

  
  # try matching to find p.value balance

  match_object = gm_out$mout
  
  mb <- MatchBalance(d.form, data = dta, match.out = match_object, print.level = 0)
  
  
  # find total time taken
  total.time <- Sys.time() - start.time
  print(sprintf("Time taken is %.3f",total.time))
  

  line <- c(toString(dtaname), total.time, gm_out$fitness, mb$AMsmallest.p.value,
                              metric,TE_finder(dtaname),match_object$est,match_object$se,pop.size, seed)

  return(line)}


for (i in 57:180){df[nrow(df)+1,] = maccabee(dtaname = names[i,1],id = names[i,2],seed=120,pop.size=10,wait.generations=5,fit.func = myfit.imb)}


gm_df = df
