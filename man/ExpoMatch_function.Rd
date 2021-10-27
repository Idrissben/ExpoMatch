\name{ExpoMatch_function}
\alias{ExpoMatch_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ExpoMatch function
}
\description{
Find the best control units for Matching exploring a large set of weights
}
\usage{
ExpoMatch_function(Tr, X, pop.size = 10, max.generations = 5, domains = c(0.5, 1.99), start.gm = TRUE, start.weights = NULL, print.level = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tr}{
The treatment assignment vector
}
  \item{X}{
The covariates
}
  \item{pop.size}{
The population size of a single generation of the genetic algorithm. Currently, it applies to both layers of the ExpoMatch Function.
}
  \item{max.generations}{
%%     ~~Describe \code{max.generations} here~~
}
  \item{domains}{
%%     ~~Describe \code{domains} here~~
}
  \item{start.gm}{
%%     ~~Describe \code{start.gm} here~~
}
  \item{start.weights}{
%%     ~~Describe \code{start.weights} here~~
}
  \item{print.level}{
%%     ~~Describe \code{print.level} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Tr, X, pop.size = 10, max.generations = 5, domains = c(0.5,
    1.99), start.gm = TRUE, start.weights = NULL, print.level = 1)
{
    start.time <- Sys.time()
    if (length(Tr) != nrow(X)) {
        return("ERROR - INCOMPATIBLE LENGTHS")
    }
    if (min(X) < 0) {
        return("CAREFUL - NEGATIVE INPUT VALUES")
    }
    if (is.na(max(X)^domains[2])) {
        return("HIGHEST POSSIBLE EXPONENT OUT OF BOUNDS")
    }
    if (!"Matching" \%in\% (.packages())) {
        library(Matching)
    }
    if (!"rgenoud" \%in\% (.packages())) {
        library(rgenoud)
    }
    if (!is.null(start.weights) & start.gm) {
        print("Deactivating starting genmatch, because you inputted starting weights.")
        start.gm <- FALSE
    }
    n.obs <- length(Tr)
    n.var <- ncol(X)
    if (start.gm) {
        genout.start <- GenMatch(Tr = Tr, X = X, pop.size = pop.size *
            50, max.generations = max.generations * 2, BalanceMatrix = X,
            hard.generation.limit = TRUE)
        print("The initial GenMatch run determined the following variable weights:")
        print(genout.start$par)
        start.weights <- genout.start$par
        n.p.values <- length(genout.start$values)
    }
    if (is.null(start.weights)) {
        start.weights <- rep(1, n.var)
    }
    dom <- cbind(rep(domains[1], n.var), rep(domains[2], n.var))
    GenMatchWrapper <- function(exponents) {
        if (print.level == 2) {
            print(exponents)
        }
        XN <- X
        for (i in c(1:n.var)) {
            XN[, i] <- XN[, i]^exponents[i]
        }
        genout <- GenMatch(Tr = Tr, X = XN, print.level = 1,
            project.path = paste(tempdir(), "/genoud.txt", sep = ""),
            pop.size = pop.size, max.generations = max.generations,
            BalanceMatrix = X, starting.values = start.weights)
        return(genout$value[1])
    }
    genoudout <- genoud(GenMatchWrapper, nvars = n.var, max = TRUE,
        pop.size = pop.size, max.generations = max.generations,
        Domains = dom, boundary.enforcement = 2, starting.values = rep(1,
            n.var))
    file_data <- read.delim(paste(tempdir(), "/genoud.txt", sep = ""),
        skip = 1, header = FALSE, nrows = 1)
    best.weights <- file_data[1, (2 + n.p.values):(1 + n.p.values +
        n.var)]
    best.weights <- as.numeric(best.weights[1, ])
    XM <- X
    for (i in c(1:n.var)) {
        XM[, i] <- XM[, i]^genoudout$par[i]
    }
    if (print.level == 2) {
        print(head(XM))
    }
    genout.fin <- GenMatch(Tr = Tr, X = XM, pop.size = pop.size *
        50, max.generations = max.generations * 2, BalanceMatrix = X,
        starting.values = best.weights)
    mout.fin <- Match(Tr = Tr, X = XM, Weight.matrix = genout.fin)
    end.time <- Sys.time()
    outputlist <- list(mout = mout.fin, genout = genout.fin,
        matches = genout.fin$matches, pvalues = genout.fin$value,
        weights = genout.fin$par, exponents = genoudout$par,
        time = end.time - start.time)
    if (start.gm) {
        outputlist[["start.weights"]] = genout.start$par
        outputlist[["start.pvalues"]] = genout.start$value
    }
    print("")
    print("##########################")
    print("")
    if (print.level == 2) {
        print(genout.fin$matches)
    }
    print(sprintf("p-values: \%s", genout.fin$value))
    print(sprintf("weights: \%s", genout.fin$par))
    print(sprintf("exponents: \%s", genoudout$par))
    print(sprintf("time taken: \%s", end.time - start.time))
    return(outputlist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
