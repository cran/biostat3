eform <- function (object, ...) 
  UseMethod("eform")
eform.default <- function(object, parm, level = 0.95, method=c("Delta","Profile"), name="exp(beta)", ...) {
  method <- match.arg(method)
  if (missing(parm))
      parm <- TRUE
  if (method == "Profile") class(object) <- c(class(object),"glm")
  estfun <- switch(method, Profile = confint, Delta = stats::confint.default)
  val <- exp(cbind(coef = coef(object), estfun(object, level = level)))
  colnames(val) <- c(name,colnames(val)[-1])
  val[parm, ]
}
irr <- function(..., name = "IRR") eform(..., name = name)
or <- function(..., name = "OR") eform(..., name = name)
methods::setGeneric("eform")
