#' @title plot.RJSMC:  Plotting results from RJSMC::SMC routine
#' @description Plots an area plot over time with the areas corresponding to probabilities for the different discrete events. If 'truth' is available, also a piecewise horizontal line with levels equal to the state number is displayed within the same plot (with levels described on the right horizontal axis)
#' @param obj Object of class RJSMC, the output from RJSMC::SMC
#' @param origin Starting date for data
#' @param pl Equal to V, Z, Q or F depending on which variable to pl
#' @param truth A list containing 'V', a vector with true states within each segment and 'B', a vector containing the starting points for each segment
#' @param clnames A vector giving names on the different classes. If NULL (default), classes are named 0,1,...
#' @export

plot.RJSMC = function(obj,
                      origin=as.Date("2022-01-01"),
                      pl="V",
                      truth=list(B=Bvec,cl=Vvec),clnames=NULL)
{
 post_V = obj@posteriors_container_V
 post_Z = obj@posteriors_container_Z
 post_Q = obj@posteriors_container_Q
 post_F = obj@posteriors_container_F
 points = obj@points_container

 points = as.Date(points,origin=origin)


 if(pl=="V")
 {
   df = data.frame(date=points,V=post_V)
   if(!is.null(clnames))
       names(df)[-1] =paste(0:(length(clnames)-1),clnames)
   if(is.null(clnames))
       names(df)[-1] = 0:(ncol(post_V)-1)
 }

 else if (pl=="Z")
   df = data.frame(date=points,Z=post_Z)
 else if (pl=="Q")
   df = data.frame(date=points,Q=post_Q)
 else if (pl=="F")
   df = data.frame(date=points,F=post_F)

 dN = colMeans(df[,-1])
 ncl = max(which(dN>0))-1
 df = df[,1:(ncl+2)]

 # Got problems when Prob=1, scale it by factor 0.99999
 df[,-1] = df[,-1]*0.99999

 g.gath =  tidyr::gather(df,Class, Probabilities,names(df)[-1])
 g.gath <<- g.gath
 g = ggplot2::ggplot(data = g.gath,
                     ggplot2::aes(x = date,
                                  y = Probabilities,
                                  fill = Class)
                     ) + ggplot2::geom_area()

 #+ ggplot2::scale_fill_manual(values=c("gray","red","black", "pink","orange", "green", "blue", "yellow","brown" )  )
 if(!is.null(truth))
 {
    dfDate = as.POSIXct(df$date,origin=origin)
    d=as.POSIXct(truth$B,"%d",origin=origin)
     d = c(d,tail(d,1)+1)
     dat = sort(unique(g.gath$date))
     g2 = data.frame(date=dat)
     n = nrow(g2)
     g2$cl = rep(0,n)
     g2$Probabilities = rep(0.025,n)
     g2$Class = rep(0,n)
     for(i in length(truth$B):1)
     {
         g2$Class[dfDate<d[i+1]] = truth$cl[i]
     }
     g2$Class = as.character(g2$Class)
     g2$dat2 = g2$date+1
     g2 <<- g2
     d <<- d
     dfDate <<- dfDate
     g = g + geom_rect(data=g2,aes(xmin=date,xmax=dat2,ymin=1.01,ymax=1.1,fill=Class))
    }
 g
}
