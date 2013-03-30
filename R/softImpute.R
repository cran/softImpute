
softImpute=function(x, rank.max = 2,lambda=0, type=c("als","svd"),thresh = 1e-05, maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE){
   this.call=match.call()
   type=match.arg(type)
   fit=softImpute.x(x, J=rank.max,lambda,type, thresh,maxit,trace.it,warm.start,final.svd)
  attr(fit,"call")=this.call
  fit
 }
softImpute.x.matrix=function(x, J,lambda,type, thresh,maxit,trace.it,warm.start,final.svd){
  switch(type,
         "als"=simpute.als(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd),
         "svd"=simpute.svd(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd)
         )
}
softImpute.x.Incomplete=function(x, J,lambda,type, thresh,maxit,trace.it,warm.start,final.svd){
  switch(type,
         "als"=Ssimpute.als(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd),
         "svd"=Ssimpute.svd(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd)
         )
}

setGeneric("softImpute.x",softImpute.x.matrix)
setMethod("softImpute.x","Incomplete",softImpute.x.Incomplete)
