plot.process = function(model,logit="effects"){

  if (logit=="OR"){
    coef.vec.1 <- c(model$ORd1, model$ORz1)
    lower.vec.1 <- c(model$ORd1.ci[1], model$ORz1.ci[1])
    upper.vec.1 <- c(model$ORd1.ci[2], model$ORz1.ci[2])
    coef.vec.1.NM <- model$ORd1.NM
    lower.vec.1.NM <- model$ORd1.ci.NM[,1]
    upper.vec.1.NM <- model$ORd1.ci.NM[,2]
    tau.vec <- c(model$ORtau.coef, model$ORtau.ci[1], model$ORtau.ci[2])
    range.1 <- range(model$ORd1.ci[1],model$ORd1.ci.NM[,1], model$ORz1.ci[1], model$ORtau.ci[1],
                     model$ORd1.ci[2],model$ORd1.ci.NM[,2], model$ORz1.ci[2], model$ORtau.ci[2])

    coef.vec.0 <- c(model$ORd0, model$ORz0)
    lower.vec.0 <- c(model$ORd0.ci[1], model$ORz0.ci[1])
    upper.vec.0 <- c(model$ORd0.ci[2], model$ORz0.ci[2])
    coef.vec.0.NM <- model$ORd0.NM
    lower.vec.0.NM <- model$ORd0.ci.NM[,1]
    upper.vec.0.NM <- model$ORd0.ci.NM[,2]
    range.0 <- range(model$ORd0.ci[1],model$ORd0.ci.NM[,1], model$ORz0.ci[1], model$ORtau.ci[1],
                     model$ORd0.ci[2],model$ORd0.ci.NM[,2], model$ORz0.ci[2], model$ORtau.ci[2])

    coef.vec.avg <- c(model$ORd.avg, model$ORz.avg)
    lower.vec.avg <- c(model$ORd.avg.ci[1], model$ORz.avg.ci[1])
    upper.vec.avg <- c(model$ORd.avg.ci[2], model$ORz.avg.ci[2])
    coef.vec.avg.NM <- model$ORd.avg.NM
    lower.vec.avg.NM <- model$ORd.avg.ci.NM[,1]
    upper.vec.avg.NM <- model$ORd.avg.ci.NM[,2]
    range.avg <- range(model$ORd.avg.ci[1], model$ORd.avg.ci.NM[,1],model$ORz.avg.ci[1], model$ORtau.ci[1],
                       model$ORd.avg.ci[2],model$ORd.avg.ci.NM[,2], model$ORz.avg.ci[2], model$ORtau.ci[2])
  }
  else if (logit=="logOR"){
    coef.vec.1 <- c(model$logORd1, model$logORz1)
    lower.vec.1 <- c(model$logORd1.ci[1], model$logORz1.ci[1])
    upper.vec.1 <- c(model$logORd1.ci[2], model$logORz1.ci[2])
    coef.vec.1.NM <- model$logORd1.NM
    lower.vec.1.NM <- model$logORd1.ci.NM[,1]
    upper.vec.1.NM <- model$logORd1.ci.NM[,2]
    tau.vec <- c(model$logORtau.coef, model$logORtau.ci[1], model$logORtau.ci[2])
    range.1 <- range(model$logORd1.ci[1],model$logORd1.ci.NM[,1], model$logORz1.ci[1], model$logORtau.ci[1],
                     model$logORd1.ci[2],model$logORd1.ci.NM[,2], model$logORz1.ci[2], model$logORtau.ci[2])

    coef.vec.0 <- c(model$logORd0, model$logORz0)
    lower.vec.0 <- c(model$logORd0.ci[1], model$logORz0.ci[1])
    upper.vec.0 <- c(model$logORd0.ci[2], model$logORz0.ci[2])
    coef.vec.0.NM <- model$logORd0.NM
    lower.vec.0.NM <- model$logORd0.ci.NM[,1]
    upper.vec.0.NM <- model$logORd0.ci.NM[,2]
    range.0 <- range(model$logORd0.ci[1],model$logORd0.ci.NM[,1], model$logORz0.ci[1], model$logORtau.ci[1],
                     model$logORd0.ci[2],model$logORd0.ci.NM[,2], model$logORz0.ci[2], model$logORtau.ci[2])

    coef.vec.avg <- c(model$logORd.avg, model$logORz.avg)
    lower.vec.avg <- c(model$logORd.avg.ci[1], model$logORz.avg.ci[1])
    upper.vec.avg <- c(model$logORd.avg.ci[2], model$logORz.avg.ci[2])
    coef.vec.avg.NM <- model$logORd.avg.NM
    lower.vec.avg.NM <- model$logORd.avg.ci.NM[,1]
    upper.vec.avg.NM <- model$logORd.avg.ci.NM[,2]
    range.avg <- range(model$logORd.avg.ci[1], model$logORd.avg.ci.NM[,1],model$logORz.avg.ci[1], model$logORtau.ci[1],
                       model$logORd.avg.ci[2],model$logORd.avg.ci.NM[,2], model$logORz.avg.ci[2], model$logORtau.ci[2])
  }
  else{
    coef.vec.1 <- c(model$d1, model$z1)
    lower.vec.1 <- c(model$d1.ci[1], model$z1.ci[1])
    upper.vec.1 <- c(model$d1.ci[2], model$z1.ci[2])
    coef.vec.1.NM <- model$d1.NM
    lower.vec.1.NM <- model$d1.ci.NM[,1]
    upper.vec.1.NM <- model$d1.ci.NM[,2]
    tau.vec <- c(model$tau.coef, model$tau.ci[1], model$tau.ci[2])
    range.1 <- range(model$d1.ci[1],model$d1.ci.NM[,1], model$z1.ci[1], model$tau.ci[1],
                     model$d1.ci[2],model$d1.ci.NM[,2], model$z1.ci[2], model$tau.ci[2])

    coef.vec.0 <- c(model$d0, model$z0)
    lower.vec.0 <- c(model$d0.ci[1], model$z0.ci[1])
    upper.vec.0 <- c(model$d0.ci[2], model$z0.ci[2])
    coef.vec.0.NM <- model$d0.NM
    lower.vec.0.NM <- model$d0.ci.NM[,1]
    upper.vec.0.NM <- model$d0.ci.NM[,2]
    range.0 <- range(model$d0.ci[1],model$d0.ci.NM[,1], model$z0.ci[1], model$tau.ci[1],
                     model$d0.ci[2],model$d0.ci.NM[,2], model$z0.ci[2], model$tau.ci[2])

    coef.vec.avg <- c(model$d.avg, model$z.avg)
    lower.vec.avg <- c(model$d.avg.ci[1], model$z.avg.ci[1])
    upper.vec.avg <- c(model$d.avg.ci[2], model$z.avg.ci[2])
    coef.vec.avg.NM <- model$d.avg.NM
    lower.vec.avg.NM <- model$d.avg.ci.NM[,1]
    upper.vec.avg.NM <- model$d.avg.ci.NM[,2]
    range.avg <- range(model$d.avg.ci[1], model$d.avg.ci.NM[,1],model$z.avg.ci[1], model$tau.ci[1],
                       model$d.avg.ci[2],model$d.avg.ci.NM[,2], model$z.avg.ci[2], model$tau.ci[2])
  }






  return(list(coef.vec.1 = coef.vec.1, lower.vec.1 = lower.vec.1, upper.vec.1 = upper.vec.1,coef.vec.1.NM=coef.vec.1.NM,lower.vec.1.NM=lower.vec.1.NM, upper.vec.1.NM= upper.vec.1.NM,
              coef.vec.0 = coef.vec.0, lower.vec.0 = lower.vec.0, upper.vec.0 = upper.vec.0,coef.vec.0.NM=coef.vec.0.NM,lower.vec.0.NM=lower.vec.0.NM, upper.vec.0.NM= upper.vec.0.NM,
              coef.vec.avg = coef.vec.avg, lower.vec.avg = lower.vec.avg, upper.vec.avg = upper.vec.avg,coef.vec.avg.NM=coef.vec.avg.NM,lower.vec.avg.NM=lower.vec.avg.NM, upper.vec.avg.NM= upper.vec.avg.NM,
              tau.vec = tau.vec,
              range.1 = range.1, range.0 = range.0, range.avg = range.avg))
}
