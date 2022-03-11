### FUNCTION TO GET ESTIMATE
get_estimate <- function(matches, Y) {
  
  sum_treated = sum(Y[matches[,1]]*matches[,3])
  sum_control = sum(Y[matches[,2]]*matches[,3])
  
  return((sum_treated - sum_control) / sum(matches[,3]))
}

### FUNCTION TO COMPARE

compare <- function(Tr = Tr,X = X, Y = Y, matches1 = null,matches2 = null,plot=TRUE) {
  
  count = 0

  for (i in 1:sum(Tr)){
    if (matches1[which(matches1[,1] == i),2] == matches2[which(matches2[,1] == i),2]) count = count +1
  }
  
  diffs.1 = colSums(BM[matches1[,1],]-BM[matches1[,2],])^2
  
  diffs.2 = colSums(BM[matches2[,1],]-BM[matches2[,2],])^2
  
  
  table = rbind(diffs.1,diffs.2)
  
  table
            
  colnames(table) = c("age","re74","re75","education","nodegree","married","black","hispanic")
            
  tes = c(get_estimate(matches1,Y=Y),get_estimate(matches2,Y=Y))
  
  sprintf('Matches 1 resulted in an ATT of: %f',get_estimate(matches1,Y=Y))
  sprintf('Matches 2 resulted in an ATT of: %f',get_estimate(matches2,Y=Y))
  
  if(plot){
    data1 = as.data.frame(t(diffs.1))
    data2 = as.data.frame(t(diffs.2))
    
    colnames(data1) = c("age","re74","re75","education","nodegree","married","black","hispanic")
    colnames(data2) = c("age","re74","re75","education","nodegree","married","black","hispanic")
    
    radarchart(rbind(rep(max(data1),8), rep(0,8),data1),axistype = 2)
    radarchart(rbind(rep(max(data2),8), rep(0,8),data2),axistype = 2)
    
  }
  return(list(table,count,tes))
}

### FUNCTION TO GET PLOT
get_plot = function(Tr = Tr,X = X, Y = Y, matches1 = null){
  
  diffs.1 = colSums(BM[matches1[,1],]-BM[matches1[,2],])^2
  data1 = as.data.frame(t(diffs.1))
  colnames(data1) = c("age","re74","re75","education","nodegree","married","black","hispanic")
  radarchart(rbind(rep(max(data1),8), rep(0,8),data1),axistype = 2)
  
}
