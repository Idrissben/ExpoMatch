install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)

etwd("parent_directory")

setwd("./ExpoMatch")
document()

devtools::document()

setwd("..")
install("ExpoMatch")
?ExpoMatch_function
