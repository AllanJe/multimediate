DichotomiePD = function(fPD,e,covMiMj,pMi,bornep,p,pMj,borneq,q){
  u=-1
  v=1
  K=0
  while (abs(u-v)>=e/2) {
    m=(u+v)/2
    fu=fPD(rho=u,covMiMj,pMi,bornep,p,pMj,borneq,q)
    fm=fPD(rho=m,covMiMj,pMi,bornep,p,pMj,borneq,q)
    if (fu*fm<0)
    {v=m}
    else
    {u=m}
    K=K+1
  }
  return(u)
}
