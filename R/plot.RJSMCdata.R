#' @title function for plotting RJSMC data
#'
#' @param  d list containing two components,
#' \describe{
#' \item{time}{Vector of time-points of class 'Date'}
#' \item{labels}{Vector of labels, categorical (or integer)}
#' }
#' @return list with following objects
#' \describe{
#' \item{minimum_n}{minimum number of messages allowed in a segment}
#'}
#' @export
plot.RJSMCdata = function(d)
{
d$labels = droplevels(d$labels)
tab = table(x$label)
n = nrow(x)
m = length(tab)

x = list()
x$tid = as.POSIXct(x$tmsp,format="%Y.%m.%d %H:%M:%S")
x$dum = sample(1:m,n,replace=T)
par(mar=c(5,8,4,2)+0.1)
plot(x$tid,x$dum,type="n",xlab="Time",ylab="",yaxt="n")
ev = unique(x$label)
ntab = rep(NA,length(ev))
for(j in 1:m)
{
  sub.x = x[x$label==ev[j],]
  ntab[j] = nrow(sub.x)
}
d = data.frame(ev=ev,ntab=ntab)
d = d[order(d$ntab,decreasing=T),]
ev2 = d$ev
for(j in 1:m)
{
  sub.x = x[x$label==ev2[j],]
  points(sub.x$tid,rep(j,nrow(sub.x)),col=j,pch=16)
}
axis(2,labels=ev2,at=1:m,las=1,cex.axis=0.75)
}

