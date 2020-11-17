#' summary.bm
#'
#' "summary.bm" is used to display the results of the mediation analyzes done with "blockmediate".
#'
#'
#' @param object element of the class "blockmediate".
#' @param opt a character string indicating the details of the analysis "navg" for the average causal effects for t=0,1 and "avg"  for the average causal effects.
#' @param logit a character string indicating, when the outcome is binary, the scale of the average causal effects. "effects" for average causal effects, " OR" average causal effects on OR scale, "logOR" average causal effects on logOR scale and "all" for all scale.
#' @param ... additional arguments affecting the summary produced
#'
#' @return table summarizing the causal analysis
#' @export
#'

summary.bm=function(object,opt="navg",logit="all",...){


    nom.navg=c("ACME.joint.treat","PM(treat)","ACME.joint.control","PM(control)",paste(c("ACME.treat.","PM(treat).","ACME.control.","PM(control)."),"Block.",rep(1:max(object$clust),each=4),sep=""),"ADE.treat","ADE.control","Total Effect")
    nom.avg=c("ACME.joint","PM.joint",paste(c("ACME.","PM."),"Block.",rep(1:max(object$clust),each=2),sep=""),"ADE","Total Effect")


  if (max(object$clust)>1){
    navg=data.frame("."        =nom.navg,
                    Estimation=round(c(object$d1,object$n1,object$d0,object$n0,triout.NM(object$d1.NM,object$n1.NM,object$d0.NM,object$n0.NM),object$z1,object$z0,object$tau.coef),4),
                    IC.inf    =round(c(object$d1.ci[1],object$n1.ci[1],object$d0.ci[1],object$n0.ci[1],triout.ci.NM(object$d1.ci.NM,object$n1.ci.NM,object$d0.ci.NM,object$n0.ci.NM)[,1],object$z1.ci[1],object$z0.ci[1],object$tau.ci[1]),4),
                    IC.sup    =round(c(object$d1.ci[2],object$n1.ci[2],object$d0.ci[2],object$n0.ci[2],triout.ci.NM(object$d1.ci.NM,object$n1.ci.NM,object$d0.ci.NM,object$n0.ci.NM)[,2],object$z1.ci[2],object$z0.ci[2],object$tau.ci[2]),4),
                    P.val     =round(c(object$d1.p,object$n1.p,object$d0.p,object$n0.p,triout.NM(object$d1.p.NM,object$n1.p.NM,object$d0.p.NM,object$n0.p.NM),object$z1.p,object$z0.p,object$tau.p),4)
    )
    avg=data.frame("."        =nom.avg,
                   Estimation=round(c(object$d.avg,object$n.avg,triout.avg.NM(object$d.avg.NM,object$n.avg.NM),object$z.avg,object$tau.coef),4),
                   IC.inf    =round(c(object$d.avg.ci[1],object$n.avg.ci[1],triout.ci.avg.NM(object$d.avg.ci.NM,object$n.avg.ci.NM)[,1],object$z.avg.ci[1],object$tau.ci[1]),4),
                   IC.sup    =round(c(object$d.avg.ci[2],object$n.avg.ci[2],triout.ci.avg.NM(object$d.avg.ci.NM,object$n.avg.ci.NM)[,2],object$z.avg.ci[2],object$tau.ci[2]),4),
                   P.val     =round(c(object$d.avg.p,object$n.avg.p,triout.avg.NM(object$d.avg.p.NM,object$n.avg.p.NM),object$z.avg.p,object$tau.p),4)
    )

    if (!is.null(object$model.y$family)){
      if (object$model.y$family$link=="logit"){
        ORnavg=data.frame("."        =paste("OR",nom.navg),
                          Estimation=round(c(object$ORd1,object$ORn1,object$ORd0,object$ORn0,triout.NM(object$ORd1.NM,object$ORn1.NM,object$ORd0.NM,object$ORn0.NM),object$ORz1,object$ORz0,object$ORtau.coef),4),
                          IC.inf    =round(c(object$ORd1.ci[1],object$ORn1.ci[1],object$ORd0.ci[1],object$ORn0.ci[1],triout.ci.NM(object$ORd1.ci.NM,object$ORn1.ci.NM,object$ORd0.ci.NM,object$ORn0.ci.NM)[,1],object$ORz1.ci[1],object$ORz0.ci[1],object$ORtau.ci[1]),4),
                          IC.sup    =round(c(object$ORd1.ci[2],object$ORn1.ci[2],object$ORd0.ci[2],object$ORn0.ci[2],triout.ci.NM(object$ORd1.ci.NM,object$ORn1.ci.NM,object$ORd0.ci.NM,object$ORn0.ci.NM)[,2],object$ORz1.ci[2],object$ORz0.ci[2],object$ORtau.ci[2]),4),
                          P.val     =round(c(object$ORd1.p,object$ORn1.p,object$ORd0.p,object$ORn0.p,triout.NM(object$ORd1.p.NM,object$ORn1.p.NM,object$ORd0.p.NM,object$ORn0.p.NM),object$ORz1.p,object$ORz0.p,object$ORtau.p),4)
        )
        ORavg=data.frame("."        =paste("OR",nom.avg),
                         Estimation=round(c(object$ORd.avg,object$ORn.avg,triout.avg.NM(object$ORd.avg.NM,object$ORn.avg.NM),object$ORz.avg,object$ORtau.coef),4),
                         IC.inf    =round(c(object$ORd.avg.ci[1],object$ORn.avg.ci[1],triout.ci.avg.NM(object$ORd.avg.ci.NM,object$ORn.avg.ci.NM)[,1],object$ORz.avg.ci[1],object$ORtau.ci[1]),4),
                         IC.sup    =round(c(object$ORd.avg.ci[2],object$ORn.avg.ci[2],triout.ci.avg.NM(object$ORd.avg.ci.NM,object$ORn.avg.ci.NM)[,2],object$ORz.avg.ci[2],object$ORtau.ci[2]),4),
                         P.val     =round(c(object$ORd.avg.p,object$ORn.avg.p,triout.avg.NM(object$ORd.avg.p.NM,object$ORn.avg.p.NM),object$ORz.avg.p,object$ORtau.p),4)
        )

        logORnavg=data.frame("."        =paste("logOR",nom.navg),
                             Estimation=round(c(object$logORd1,object$logORn1,object$logORd0,object$logORn0,triout.NM(object$logORd1.NM,object$logORn1.NM,object$logORd0.NM,object$logORn0.NM),object$logORz1,object$logORz0,object$logORtau.coef),4),
                             IC.inf    =round(c(object$logORd1.ci[1],object$logORn1.ci[1],object$logORd0.ci[1],object$logORn0.ci[1],triout.ci.NM(object$logORd1.ci.NM,object$logORn1.ci.NM,object$logORd0.ci.NM,object$logORn0.ci.NM)[,1],object$logORz1.ci[1],object$logORz0.ci[1],object$logORtau.ci[1]),4),
                             IC.sup    =round(c(object$logORd1.ci[2],object$logORn1.ci[2],object$logORd0.ci[2],object$logORn0.ci[2],triout.ci.NM(object$logORd1.ci.NM,object$logORn1.ci.NM,object$logORd0.ci.NM,object$logORn0.ci.NM)[,2],object$logORz1.ci[2],object$logORz0.ci[2],object$logORtau.ci[2]),4),
                             P.val     =round(c(object$logORd1.p,object$logORn1.p,object$logORd0.p,object$logORn0.p,triout.NM(object$logORd1.p.NM,object$logORn1.p.NM,object$logORd0.p.NM,object$logORn0.p.NM),object$logORz1.p,object$logORz0.p,object$logORtau.p),4)
        )
        logORavg=data.frame("."        =paste("logOR",nom.avg),
                            Estimation=round(c(object$logORd.avg,object$logORn.avg,triout.avg.NM(object$logORd.avg.NM,object$logORn.avg.NM),object$logORz.avg,object$logORtau.coef),4),
                            IC.inf    =round(c(object$logORd.avg.ci[1],object$logORn.avg.ci[1],triout.ci.avg.NM(object$logORd.avg.ci.NM,object$logORn.avg.ci.NM)[,1],object$logORz.avg.ci[1],object$logORtau.ci[1]),4),
                            IC.sup    =round(c(object$logORd.avg.ci[2],object$logORn.avg.ci[2],triout.ci.avg.NM(object$logORd.avg.ci.NM,object$logORn.avg.ci.NM)[,2],object$logORz.avg.ci[2],object$logORtau.ci[2]),4),
                            P.val     =round(c(object$logORd.avg.p,object$logORn.avg.p,triout.avg.NM(object$logORd.avg.p.NM,object$logORn.avg.p.NM),object$logORz.avg.p,object$logORtau.p),4)
        )
      }}
  }
  else {
    navg=data.frame("."        =c("ACME.treat","PM(treat)","ACME.control","PM(control)","ADE.treat","ADE.control","Total Effect"),
                    Estimation=round(c(object$d1,object$n1,object$d0,object$n0,object$z1,object$z0,object$tau.coef),4),
                    IC.inf    =round(c(object$d1.ci[1],object$n1.ci[1],object$d0.ci[1],object$n0.ci[1],object$z1.ci[1],object$z0.ci[1],object$tau.ci[1]),4),
                    IC.sup    =round(c(object$d1.ci[2],object$n1.ci[2],object$d0.ci[2],object$n0.ci[2],object$z1.ci[2],object$z0.ci[2],object$tau.ci[2]),4),
                    P.val     =round(c(object$d1.p,object$n1.p,object$d0.p,object$n0.p,object$z1.p,object$z0.p,object$tau.p),4)
    )

    avg=data.frame("."        =c("ACME","PM","ADE","Total Effect"),
                   Estimation=round(c(object$d.avg,object$n.avg,object$z.avg,object$tau.coef),4),
                   IC.inf    =round(c(object$d.avg.ci[1],object$n.avg.ci[1],object$z.avg.ci[1],object$tau.ci[1]),4),
                   IC.sup    =round(c(object$d.avg.ci[2],object$n.avg.ci[2],object$z.avg.ci[2],object$tau.ci[2]),4),
                   P.val     =round(c(object$d.avg.p,object$n.avg.p,object$z.avg.p,object$tau.p),4)
    )
    if (!is.null(object$model.y$family)){
      if (object$model.y$family$link=="logit"){
        ORnavg=data.frame("."        =c("OR.ACME.treat","OR.PM(treat)","OR.ACME.control","OR.PM(control)","OR.ADE.treat","OR.ADE.control","OR.Total Effect"),
                          Estimation=round(c(object$ORd1,object$ORn1,object$ORd0,object$ORn0,object$ORz1,object$ORz0,object$ORtau.coef),4),
                          IC.inf    =round(c(object$ORd1.ci[1],object$ORn1.ci[1],object$ORd0.ci[1],object$ORn0.ci[1],object$ORz1.ci[1],object$ORz0.ci[1],object$ORtau.ci[1]),4),
                          IC.sup    =round(c(object$ORd1.ci[2],object$ORn1.ci[2],object$ORd0.ci[2],object$ORn0.ci[2],object$ORz1.ci[2],object$ORz0.ci[2],object$ORtau.ci[2]),4),
                          P.val     =round(c(object$ORd1.p,object$ORn1.p,object$ORd0.p,object$ORn0.p,object$ORz1.p,object$ORz0.p,object$ORtau.p),4)
        )

        ORavg=data.frame("."        =c("OR.ACME","OR.PM","OR.ADE","OR.Total Effect"),
                         Estimation=round(c(object$ORd.avg,object$ORn.avg,object$ORz.avg,object$ORtau.coef),4),
                         IC.inf    =round(c(object$ORd.avg.ci[1],object$ORn.avg.ci[1],object$ORz.avg.ci[1],object$ORtau.ci[1]),4),
                         IC.sup    =round(c(object$ORd.avg.ci[2],object$ORn.avg.ci[2],object$ORz.avg.ci[2],object$ORtau.ci[2]),4),
                         P.val     =round(c(object$ORd.avg.p,object$ORn.avg.p,object$ORz.avg.p,object$ORtau.p),4)
        )

        logORnavg=data.frame("."        =c("logOR.ACME.treat","logOR.PM(treat)","logOR.ACME.control","logOR.PM(control)","logOR.ADE.treat","logOR.ADE.control","logOR.Total Effect"),
                             Estimation=round(c(object$logORd1,object$logORn1,object$logORd0,object$logORn0,object$logORz1,object$logORz0,object$logORtau.coef),4),
                             IC.inf    =round(c(object$logORd1.ci[1],object$logORn1.ci[1],object$logORd0.ci[1],object$logORn0.ci[1],object$logORz1.ci[1],object$logORz0.ci[1],object$logORtau.ci[1]),4),
                             IC.sup    =round(c(object$logORd1.ci[2],object$logORn1.ci[2],object$logORd0.ci[2],object$logORn0.ci[2],object$logORz1.ci[2],object$logORz0.ci[2],object$logORtau.ci[2]),4),
                             P.val     =round(c(object$logORd1.p,object$logORn1.p,object$logORd0.p,object$logORn0.p,object$logORz1.p,object$logORz0.p,object$logORtau.p),4)
        )
        logORavg=data.frame("."        =c("logOR.ACME","logOR.PM","logOR.ADE","logOR.Total Effect"),
                            Estimation=round(c(object$logORd.avg,object$logORn.avg,object$logORz.avg,object$logORtau.coef),4),
                            IC.inf    =round(c(object$logORd.avg.ci[1],object$logORn.avg.ci[1],object$logORz.avg.ci[1],object$logORtau.ci[1]),4),
                            IC.sup    =round(c(object$logORd.avg.ci[2],object$logORn.avg.ci[2],object$logORz.avg.ci[2],object$logORtau.ci[2]),4),
                            P.val     =round(c(object$logORd.avg.p,object$logORn.avg.p,object$logORz.avg.p,object$logORtau.p),4)
        )
      }}
  }
  if (opt=="avg"){
    if (is.null(object$model.y$family) || object$model.y$family$link=="probit")
    {res=avg}
    else{
      if (object$model.y$family$link=="logit" & logit=="all"){
        res=cbind(avg,ORavg,logORavg)
      }
      else if (object$model.y$family$link=="logit" & logit=="OR"){
        res=ORavg
      }
      else if (object$model.y$family$link=="logit" & logit=="logOR"){
        res=logORavg
      }
      else{
        res=avg
      }
    }

  }
  else{

    if (is.null(object$model.y$family) || object$model.y$family$link=="probit")
    {res=navg}
    else {
      if (object$model.y$family$link=="logit" & logit=="all"){
        res=cbind(navg,ORnavg,logORnavg)
      }
      else if (object$model.y$family$link=="logit" & logit=="OR"){
        res=ORnavg
      }
      else if (object$model.y$family$link=="logit" & logit=="logOR"){
        res=logORnavg
      }
      else{
        res=navg
      }
    }

  }
  return(res)
}
