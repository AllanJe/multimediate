plot.process=function (model)
{
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





  return(list(coef.vec.1 = coef.vec.1, lower.vec.1 = lower.vec.1, upper.vec.1 = upper.vec.1,coef.vec.1.NM=coef.vec.1.NM,lower.vec.1.NM=lower.vec.1.NM, upper.vec.1.NM= upper.vec.1.NM,
              coef.vec.0 = coef.vec.0, lower.vec.0 = lower.vec.0, upper.vec.0 = upper.vec.0,coef.vec.0.NM=coef.vec.0.NM,lower.vec.0.NM=lower.vec.0.NM, upper.vec.0.NM= upper.vec.0.NM,
              coef.vec.avg = coef.vec.avg, lower.vec.avg = lower.vec.avg, upper.vec.avg = upper.vec.avg,coef.vec.avg.NM=coef.vec.avg.NM,lower.vec.avg.NM=lower.vec.avg.NM, upper.vec.avg.NM= upper.vec.avg.NM,
              tau.vec = tau.vec,
              range.1 = range.1, range.0 = range.0, range.avg = range.avg))
}
