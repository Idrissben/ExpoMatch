############SET-UP
df = data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("name", "time taken", "Supremum value",
                                                             "Best p-value","Metric used",'Ground truth ATE',
                                                             'Estimated ATE','Standard Error','pop.size','seed'))))

names = drive_ls("https://drive.google.com/drive/u/0/folders/1j60vqxqes-vs_I094sIZgQIllv_arMFT")
### MONTE CARLO SIMULATION
################################### EXPOMATCH  
         
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


for (i in 1:180){df[nrow(df)+1,] = maccabee(dtaname = names[i,1],id = names[i,2],seed=120,pop.size=10,wait.generations=5,fit.func = myfit.imb)}
