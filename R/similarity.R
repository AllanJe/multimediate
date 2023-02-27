similarity = function(Mi, Mj) {
  if(abs(cor(Mi,Mj))>=0.99){
    out=(-1/2)*log(1-0.99^2)}
  else{
    out=(-1/2)*log(1-cor(Mi,Mj)^2)
  }
  return(out)
}
