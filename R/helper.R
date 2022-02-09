#' expit
#'
#' Helper function to 1/(1+exp(-x))
#'
#' @param x real number
#'
#' @importFrom stats plogis
#' @export

expit = function(x) {plogis(x);}


#' logit
#'
#' Helper function to log(x/(1-x))
#'
#' @param x number between 0 and 1
#'
#' @importFrom stats qlogis
#' @export

logit = function(x) {qlogis(x);}



#' Error-handling function
#'
#'
#' Borrowed from the R package simsalapar
#'
#' @param expr R expression
#'
#' @export

tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- c(W,w)
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}



