#' summary.mahi
#'
#' "summary.mahi" is used to display the results of the mediation analyzes done with "multimediate".
#'
#'
#' @param object element of the class "mahi".
#' @param opt a character string indicating the details of the analysis "navg" for the average causal effects for t=0,1 and "avg"  for the average causal effects.
#' @param logit a character string indicating, when the outcome is binary, the scale of the average causal effects. "effects" for average causal effects, " OR" average causal effects on OR scale, "logOR" average causal effects on logOR scale and "all" for all scale.
#' @param ... additional arguments affecting the summary produced
#'
#' @return table summarizing the causal analysis
#' @export
#'

summary.mahi = function(object,opt='avg',logit='all',...){
  cat('Resulst for the mahi selection with n =', object$n, "and p =", object$K,"\n \n",
      'The Kmax = ',object$Kmax, 'first potentials mediators after stability selection are : \n', object$Kmaxrank,'\n')

  if (object$step2==TRUE){
    cat('The multiple mediation analysis with the previous potentials mediators.')
    print(summary(object$multimed,opt=opt,logit=logit))

    cat('The corrected p-values using the p.adjust.method specified are : \n', object$pvalscorr)

  }

}


