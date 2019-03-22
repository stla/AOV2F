#' Run the two-way ANOVA child.
#'
#' @description Run knitr child for two-way ANOVA; to be used as inline code in
#' a knitr document.
#'
#' @param dat dataframe
#' @param factors the column names or indices of the factors
#' @param resp the column name or index of the response
#' @param factor.names the names of the factor for the report
#' (default value is the column names)
#' @param resp.name the name of the response for the report
#' (default value is the column name)
#' @param factor1.ref the reference level of factor1
#' @param factor2.ref the reference level of factor2
#' @param title study title
#'
#' @export
#' @importFrom knitr knit_child
childANOVA2F <- function(dat, factors=c(1,2), resp=3,
                         factor.names=names(dat[,factors]),
                        resp.name=names(dat[,resp,drop=FALSE]),
                        factor1.ref=levels(dat[,factors[1]])[1],
                        factor2.ref=levels(dat[,factors[2]])[1],
                        title="title") {
  if(!is.factor(dat[, factors[1]])){
    dat[, factors[1]] <- factor(dat[, factors[1]])
  }
  if(!is.factor(dat[, factors[2]])){
    dat[, factors[2]] <- factor(dat[, factors[2]])
  }
  dat <- droplevels(dat)
  myenv <- new.env()
  args <- formals(childANOVA2F)
  for(arg in names(args)) assign(arg, get(arg), envir=myenv)
  knitr::knit_child(
    system.file(package = "AOV2F", "knitr", "AOV2F_child.Rmd"),
    envir = myenv) # faire une possibilité knit
}


#' Simulate a two-way ANOVA design with random effects
#' @description Simulate a two-way ANOVA dataset with random effects.
#' @param I number of levels of first factor
#' @param J number of leveles of second factor
#' @param Kmin minimum number of repeats
#' @param Kmax maximum number of repeats
#' @param p numeric vector giving the probabilities to sample between
#' \code{Kmin} and \code{Kmax}
#' @param sigmaP standard deviation of first factor
#' @param sigmaO standard deviation of second factor
#' @param sigmaPO standard deviation of interaction
#' @param sigmaE residual standard deviation
#' @param factor.names names of the two factors
#' @param resp.name name of the response
#' @param keep.intermediate keep intermediate calculations in the output
#' @return A dataframe.
#'
#' @examples
#' SimDataAV2(I=3, J=2, Kmin=0, Kmax=2, p=c(0.1,0.2))
#'
#' @export
#' @importFrom stats rnorm
SimDataAV2 <- function(I, J, Kmin, Kmax, p=NULL, mu=0, sigmaP=1, sigmaO=1,
                       sigmaPO=1, sigmaE=1, factor.names=c("Operator","Part"),
                       resp.name="y", keep.intermediate=FALSE){
  Operator <- rep(1:J, each=I)
  Oj <- rep(rnorm(J, 0, sigmaO), each=I)
  Part <- rep(1:I, times=J)
  Pi <- rep(rnorm(I, 0, sigmaP), times=J)
  POij <- rnorm(I*J, 0, sigmaPO)
  simdata0 <- data.frame(Part, Operator, Pi, Oj, POij)
  simdata0$Operator <- factor(simdata0$Operator)
  levels(simdata0$Operator) <- sprintf(paste0("%0", floor(log10(J))+1, "d"), 1:J)
  simdata0$Part <- factor(simdata0$Part)
  levels(simdata0$Part) <- sprintf(paste0("%0", floor(log10(I))+1, "d"), 1:I)
  II <- 0 ; JJ <- 0
  while(II<I | JJ <J){
    if(Kmin < Kmax){
      Kij <- sample(Kmin:Kmax, I*J, replace=TRUE, prob=c(p,1-sum(p)))
    }else{
      Kij <- rep(Kmin, I*J)
    }
    simdata <- droplevels(
      as.data.frame(
        sapply(simdata0, function(v) rep(v, times=Kij), simplify=FALSE)))
    JJ <- length(levels(simdata$Operator)); II <- length(levels(simdata$Part))
  }
  Eijk <- rnorm(sum(Kij), 0, sigmaE)
  simdata <- cbind(simdata, Eijk)
  simdata[[resp.name]] <- mu + with(simdata, Oj+Pi+POij+Eijk)
  levels(simdata[,1]) <- paste0("A", levels(simdata[,1]))
  levels(simdata[,2]) <- paste0("B", levels(simdata[,2]))
  names(simdata)[1:2] <- factor.names
  if(!keep.intermediate) simdata <- simdata[,c(factor.names,resp.name)]
  simdata
}

#' Format a table of type \code{ftable}
#' @description Format a table of type \code{ftable} for HTML printing.
#' @note This function is based on \code{R2HTML:::HTML.ftable}
#'
#' @param x a table of type \code{\link{ftable}}
#' @param digits number of digits to print
#' @return A table which can be used in \code{kable}.
#'
#' @export
format_ftable <- function(x, digits = getOption("digits"))
{
  if (!inherits(x, "ftable"))
    stop("x must be an `ftable'")
  ox <- x
  makeLabels <- function(lst) {
    lens <- sapply(lst, length)
    cplensU <- c(1, cumprod(lens))
    cplensD <- rev(c(1, cumprod(rev(lens))))
    y <- NULL
    for (i in rev(seq(along = lst))) {
      ind <- 1 + seq(from = 0, to = lens[i] - 1) * cplensD[i+1]
      tmp <- character(length = cplensD[i])
      tmp[ind] <- lst[[i]]
      y <- cbind(rep(tmp, times = cplensU[i]), y)
    }
    y
  }
  makeNames <- function(x) {
    nmx <- names(x)
    if (is.null(nmx))
      nmx <- rep("", length = length(x))
    nmx
  }
  xrv <- attr(x, "row.vars")
  xcv <- attr(x, "col.vars")
  LABS <- cbind(rbind(matrix("", nrow = length(xcv), ncol = length(xrv)),
                      makeNames(xrv), makeLabels(xrv)),
                c(makeNames(xcv), rep("", times = nrow(x) + 1)))
  DATA <- rbind(t(makeLabels(xcv)), rep("", times = ncol(x)),
                format(unclass(x), digits = digits))
  cbind(apply(LABS, 2, format, justify = "left"),
        apply(DATA, 2, format, justify = "right"))
}
