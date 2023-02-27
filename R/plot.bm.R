#' plot.bm
#'
#' "plot.bm" is used to display the results of the mediation analyzes done with "blockmediate".
#'
#'
#' @param x element of the class "bm".
#' @param treatment a character string indicating the baseline treatment value of the estimated causal mediation effect and direct effect to plot. Can be either "control", "treated", "average" or "three". If "NULL"(default), three sets of estimates are plotted.
#' @param logit a character string indicating, when the outcome is binary, the scale of the average causal effects. "effects" for average causal effects, " OR" average causal effects on OR scale, "logOR" average causal effects on logOR scale.
#' @param labels a vector of character strings indicating the labels for estimated effects. The default labels wiil be used if NULL.
#' @param effect.type a vector indicating which quantities of interest to plot. Default is to plot all three quantities (indirects, direct and total effects).
#' @param xlim range of the horizontal axis.
#' @param ylim range of the vertical axis.
#' @param xlab label of the horizontal axis.
#' @param ylab label of the vertical axis.
#' @param main main title.
#' @param is.legend a logical value indicating the presence of the legend on top right. If the legend obscures some results, is.legend can be switched to "FALSE" or the x-axis can be changed manually with the parameters xlim.
#' @param lwd width of the horizontal bars for confidence intervals.
#' @param cex size of the dots for point estimates.
#' @param col color of the dots and horizontal bars for the estimates.
#' @param ... additional parameters passed to 'plot'.
#'
#'
#' @return plot summarizing the causal analysis
#'
#' @export
#'
#' @import graphics
#'
#'

plot.bm = function(x, treatment = NULL,logit="logOR", labels = NULL, effect.type = c("indirect","direct", "total"), xlim = NULL, ylim = NULL, xlab = "", ylab = "", main = "Estimates and confidence intervals",is.legend=TRUE, lwd = 1.5, cex = 0.85, col = "black",...){

  effect.type <- match.arg(effect.type, several.ok = TRUE)
  IND <- "indirect" %in% effect.type
  DIR <- "direct" %in% effect.type
  TOT <- "total" %in% effect.type
  if (is.null(treatment)) {
    treatment <- c(0, 1, 2)
  }
  else {
    treatment <- switch(treatment, control = 0, treated = 1, average = 2,
                        three = c(0, 1, 2))
  }

  if (logit=="logOR"){
    param <- plot.process(x,logit="logOR")
    main=paste(main, "on log OR scale")
    vline <- 0
  }
  else if (logit=="OR"){
    param <- plot.process(x,logit="OR")
    main=paste(main, "on OR scale")
    vline <- 1
  }
  else{
    param <- plot.process(x,logit="effects")
    vline <- 0
  }

  NM=max(x$clust)
  y.axis <- (IND*(NM+1) + DIR + TOT):1
  if (is.null(xlim)) {
    plusx=c(-.05,0.05)
    if (length(treatment) > 1) {
      xlim <- range(param$range.1, param$range.0) + plusx
    }
    else if (treatment == 1) {
      xlim <- param$range.1 + plusx
    }
    else if (treatment == 0){
      xlim <- param$range.0 + plusx
    }
    else{
      xlim <- param$range.avg + plusx
    }
  }
  if (is.null(ylim)) {
    ylim <- c(min(y.axis) - 0.5, max(y.axis) + 0.5)
  }
  plot(rep(0, IND*(NM+1) + DIR + TOT), y.axis, type = "n", xlab = xlab,
       ylab = ylab, yaxt = "n", xlim = xlim, ylim = ylim, main = main,
       ...)
  if (length(treatment) == 1) {
    adj <- c(0,0,0)
  }
  else  {
    adj <- c(-0.2,0,0.2)
  }
  if (1 %in% treatment) {

    if (IND && DIR) {


      points(param$coef.vec.1.NM, y.axis[1:NM] + adj[1], type = "p",
             pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1.NM, y.axis[1:NM] + adj[1], param$upper.vec.1.NM,
               y.axis[1:NM] + adj[1], lwd = lwd, col = col)

      points(param$coef.vec.1, y.axis[(NM+1):(NM+2)] + adj[1], type = "p",
             pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1, y.axis[(NM+1):(NM+2)] + adj[1], param$upper.vec.1,
               y.axis[(NM+1):(NM+2)] + adj[1], lwd = lwd, col = col)


    }
    if (IND && !DIR) {

      points(param$coef.vec.1.NM, y.axis[1:NM] + adj[1], type = "p",
             pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1.NM, y.axis[1:NM] + adj[1], param$upper.vec.1.NM,
               y.axis[1:NM] + adj[1], lwd = lwd, col = col)

      points(param$coef.vec.1[1], y.axis[NM+1] + adj[1], type = "p",
             pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1[1], y.axis[NM+1] + adj[1], param$upper.vec.1[1],
               y.axis[NM+1] + adj[1], lwd = lwd, col = col)
    }
    if (!IND && DIR) {
      points(param$coef.vec.1[2], y.axis[1] + adj[1], type = "p",
             pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1[2], y.axis[1] + adj[1], param$upper.vec.1[2],
               y.axis[1] + adj[1], lwd = lwd, col = col)
    }
  }
  if (0 %in% treatment) {
    if (IND && DIR) {

      points(param$coef.vec.0.NM, y.axis[1:NM] + adj[3], type = "p",
             pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0.NM, y.axis[1:NM] + adj[3], param$upper.vec.0.NM,
               y.axis[1:NM] + adj[3], lwd = lwd, lty = 3, col = col)

      points(param$coef.vec.0, y.axis[(NM+1):(NM+2)] + adj[3], type = "p",
             pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0, y.axis[(NM+1):(NM+2)] + adj[3], param$upper.vec.0,
               y.axis[(NM+1):(NM+2)] + adj[3], lwd = lwd, lty = 3, col = col)
    }
    if (IND && !DIR) {

      points(param$coef.vec.0.NM, y.axis[1:NM] + adj[3], type = "p",
             pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0.NM, y.axis[1:NM] + adj[3], param$upper.vec.0.NM,
               y.axis[1:NM] + adj[3], lwd = lwd, lty = 3, col = col)

      points(param$coef.vec.0[1], y.axis[NM+1] + adj[3], type = "p",
             pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0[1], y.axis[NM+1] + adj[3], param$upper.vec.0[1],
               y.axis[NM+1] + adj[3], lwd = lwd, lty = 3, col = col)
    }
    if (!IND && DIR) {
      points(param$coef.vec.0[2], y.axis[1] + adj[3], type = "p",
             pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0[2], y.axis[1] + adj[3], param$upper.vec.0[2],
               y.axis[1] + adj[3], lwd = lwd, lty = 3, col = col)
    }
  }
  if (2 %in% treatment) {
    if (IND && DIR) {

      points(param$coef.vec.avg.NM, y.axis[1:NM] + adj[2], type = "p",
             pch = 10, cex = cex, col = col)
      segments(param$lower.vec.avg.NM, y.axis[1:NM] + adj[2], param$upper.vec.avg.NM,
               y.axis[1:NM] + adj[2], lwd = lwd, lty = 5, col = col)

      points(param$coef.vec.avg, y.axis[(NM+1):(NM+2)] + adj[2], type = "p",
             pch = 10, cex = cex, col = col)
      segments(param$lower.vec.avg, y.axis[(NM+1):(NM+2)] + adj[2], param$upper.vec.avg,
               y.axis[(NM+1):(NM+2)] + adj[2], lwd = lwd, lty = 5, col = col)
    }
    if (IND && !DIR) {

      points(param$coef.vec.avg.NM, y.axis[1:NM] + adj[2], type = "p",
             pch = 10, cex = cex, col = col)
      segments(param$lower.vec.avg.NM, y.axis[1:NM] + adj[2], param$upper.vec.avg.NM,
               y.axis[1:NM] + adj[2], lwd = lwd, lty = 5, col = col)

      points(param$coef.vec.avg[1], y.axis[NM+1] + adj[2], type = "p",
             pch = 10, cex = cex, col = col)
      segments(param$lower.vec.avg[1], y.axis[NM+1] + adj[2], param$upper.vec.avg[1],
               y.axis[NM+1] + adj[2], lwd = lwd, lty = 5, col = col)
    }
    if (!IND && DIR) {
      points(param$coef.vec.avg[2], y.axis[1] + adj[2], type = "p",
             pch = 10, cex = cex, col = col)
      segments(param$lower.vec.avg[2], y.axis[1] + adj[2], param$upper.vec.avg[2],
               y.axis[1] + adj[2], lwd = lwd, lty = 5, col = col)
    }
  }
  if (TOT) {
    points(param$tau.vec[1], 1, type = "p", pch = 19, cex = cex,
           col = col)
    segments(param$tau.vec[2], 1, param$tau.vec[3], 1, lwd = lwd,
             col = col)
  }
  if (is.null(labels)) {
    labels <- c(paste("ACME\nBlock ",1:max(x$clust),sep=""),"ACME\nJoint", "ADE", "Total\nEffect")[c(rep(IND,NM+1), DIR,TOT)]
  }
  axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE,
       ...)
  abline(v = vline, lty = 2)
  if (is.legend){
    legend("topright",legend=c("control value","average","treated value"),pch=c(1,10,19),lty=c(3,5,1),title="Effects calculated with")

  }
}
