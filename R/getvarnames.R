#' getvarnames
#'
#'
#' 'getvarnames' is used to get all variable names from a regression model
#'
#'@param formula the formula from which to extract variable names
#'@param data data frame with variables from the formula
#'@return 'getvarnames' returns a list with 'varnames' (referring to all variable names), 'xvar' (referring to the predictors), and 'yvar' (referring to the outcome)

getvarnames = function(formula, data = NULL)
{
  if (is.character(formula))
    return(list(varnames=formula, xvar=formula, yvar=NULL))
  if (is.null(formula)) return(list(varnames=NULL, xvar=NULL, yvar=NULL))

  formula <- formula(formula)
  lyv <- NULL
  lxv <- lvnm <- all.vars(formula[1:2])
  if (length(formula)==3) {
    lyv <- lxv
    lxv <- all.vars(formula[-2])
    if ("." %in% lxv) {
      if (length(data)==0)
        stop("!getvarnames! '.' in formula and no 'data'")
      lform <- formula(terms(formula, data=data))
      lxv <- all.vars(lform[-2])
    }
    lvnm <- c(lxv, lvnm)
  }
  return(list(varnames=lvnm, xvar=lxv, yvar=lyv))
}
