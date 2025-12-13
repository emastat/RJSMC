#' @title plot.RJSMC:  Plotting results from RJSMC::SMC routine
#' @description Plots an area plot over time with the areas corresponding to probabilities for the different discrete events
#' @param obj Object of class RJSMC, the output from RJSMC::SMC
#' @param origin Starting date for data
#' @param pl Equal to V, Z, Q or F depending on which variable to pl
#' @export

plot.RJSMC = function(obj,
                      origin=as.Date("2022-01-01"),
                      pl="V",
                      truth=list(B=Bvec,cl=Vvec))
{
 post_V = obj@posteriors_container_V
 post_Z = obj@posteriors_container_Z
 post_Q = obj@posteriors_container_Q
 post_F = obj@posteriors_container_F
 points = obj@points_container

 points = as.Date(points*100,"%d",origin=origin)


 if(pl=="V")
   df = data.frame(date=points,V=post_V)
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

 #  mutate(date = as.Date(date)) %>%
 g.gath =  tidyr::gather(df,variable, Probabilities,names(df)[-1])
 if(!is.null(truth))
 {
    g.gath$cl = NA
    d=as.Date(truth$B*100,"%d",origin=origin)
    d = c(d,tail(d,1)+1)
    for(i in length(truth$B):1)
    {
       g.gath$cl[df$date<d[i+1]] = Vvec[i]
    }
 }
 g = ggplot2::ggplot(data = g.gath,
                     ggplot2::aes(x = date,
                                  y = Probabilities,
                                  fill = variable)
                     ) +
   ggplot2::geom_area() +
     ggplot2::scale_fill_manual(values=c("gray","red","black", "pink","orange", "green", "blue", "yellow","brown" )  )
   #+ ylim(0,1)+geom_ine(aes(y=y,x=d,))
 if(!is.null(truth))
 {
    #ncl = 5
     g = g + ggplot2::geom_line(
        ggplot2::aes(y=(cl)/ncl)) +
        ggplot2::scale_y_continuous(sec.axis=ggplot2::sec_axis(~.*ncl+1),
                                    name="True class")
 }
 g.gath.save <<- g.gath
 g
}
