risk = function(h){
  sqrt(exp(h))
}

s = function(h,s0){
  ( sqrt(exp(h))/2) * s0
}

f = function(h, s0){
  return(c( sqrt(exp(h)), ( sqrt(exp(h))/2) * s0))
}