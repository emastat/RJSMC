#' @title plot.RJSMC:  Plotting results from RJSMC::SMC routine
#' @description Plots an area plot over time with the areas corresponding to probabilities for the different discrete events. If 'truth' is available, also a piecewise horizontal line with levels equal to the state number is displayed within the same plot (with levels described on the right horizontal axis). If 'observations' is provided, dots for each observation are displayed above the truth rectangles.
#' @param obj Object of class RJSMC, the output from RJSMC::SMC
#' @param origin Starting date for data (used only when time_to_date=TRUE)
#' @param pl Equal to V, Z, Q or F depending on which variable to pl
#' @param truth A list containing 'V', a vector with true states within each segment and 'B', a vector containing the starting points for each segment
#' @param clnames A vector giving names on the different classes. If NULL (default), classes are named 0,1,...
#' @param observations A list containing 'Tvec' (time stamps of observations) and optionally 'Yvec' (message labels for coloring). If NULL (default), no observation dots are displayed.
#' @param time_to_date Logical. If TRUE (default), converts time values to Date format for better visualization. If FALSE, displays time in its original numeric format (useful for debugging).
#' @export

plot.RJSMC = function(obj,
                     origin=as.Date("2022-01-01"),
                     pl="V",
                     truth=NULL,clnames=NULL,
                     observations=NULL,
                     time_to_date=TRUE)
{
 post_V = obj@posteriors_container_V
 post_Z = obj@posteriors_container_Z
 post_Q = obj@posteriors_container_Q
 post_F = obj@posteriors_container_F
 points = obj@points_container

 # Convert to Date format only if time_to_date is TRUE
 if(time_to_date) {
   points = as.Date(points, origin=origin)
 } else {
   # Keep original numeric time values
   points = points
 }


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

 # Normalize each row to sum to 1 to ensure probabilities are valid
 # This handles any numerical issues where rows might not sum exactly to 1
 row_sums = rowSums(df[,-1])
 # Avoid division by zero (shouldn't happen, but safety check)
 row_sums[row_sums == 0] = 1
 df[,-1] = df[,-1] / row_sums

 g.gath =  tidyr::gather(df,Class, Probabilities,names(df)[-1])
 g.gath <<- g.gath
 g = ggplot2::ggplot(data = g.gath,
                     ggplot2::aes(x = date,
                                  y = Probabilities,
                                  fill = Class)
                     ) + ggplot2::geom_area()
 
 # Set x-axis label based on time format
 if(time_to_date) {
   g = g + ggplot2::xlab("Date")
 } else {
   g = g + ggplot2::xlab("Time")
 }

 # Get all unique classes from the main plot
 all_classes = sort(unique(as.character(g.gath$Class)))
 
 if(!is.null(truth))
 {
    # Convert truth breakpoints based on time_to_date setting
    if(time_to_date) {
      # Convert to Date (same scale as df$date)
      # truth$B contains numeric time values, convert to Date using origin
      d = as.Date(truth$B, origin=origin)
      # Add one day after the last breakpoint to define the final segment
      d = c(d, tail(d,1) + 1)
      # Rectangle width increment (one day)
      rect_width = 1
    } else {
      # Keep original numeric time values
      d = truth$B
      # Add small increment after the last breakpoint (use 1% of time range or 0.1, whichever is larger)
      time_range = max(points) - min(points)
      increment = max(0.1, time_range * 0.01)
      d = c(d, tail(d,1) + increment)
      # Rectangle width increment (small time increment)
      rect_width = increment
    }
    
    # Get unique sorted time points from the plot data
    dat = sort(unique(g.gath$date))
    g2 = data.frame(date=dat)
    n = nrow(g2)
    
    # Initialize Class column - will be assigned based on which segment each time point falls into
    g2$Class = rep(0, n)
    
    # Assign the true state class to each time point based on which segment it belongs to
    # Iterate backwards through segments to correctly assign classes
    for(i in length(truth$B):1)
    {
        # Assign truth$cl[i] to all time points that are before breakpoint d[i+1]
        # This means time points in segment i (from d[i] to d[i+1]) get class truth$cl[i]
        g2$Class[g2$date < d[i+1]] = truth$cl[i]
    }
    
    # Convert Class to character for ggplot2
    g2$Class = as.character(g2$Class)
    
    # Create end time points for rectangles
    g2$dat2 = g2$date + rect_width
    
    # Debug assignments (can be removed in production)
    g2 <<- g2
    d <<- d
    
    # Add rectangles to the plot showing true state sequence
    g = g + ggplot2::geom_rect(data=g2, 
                                ggplot2::aes(xmin=date, xmax=dat2, ymin=1.01, ymax=1.1, fill=Class),
                                inherit.aes = FALSE)
    
    # Update all_classes to include truth classes
    all_classes = sort(unique(c(all_classes, as.character(g2$Class))))
    }
 
 # Create color mapping: class "0" is always white, others get default colors
 # Generate colors for all classes except "0"
 n_classes = length(all_classes)
 color_values = grDevices::rainbow(n_classes)
 names(color_values) = all_classes
 
 # Set class "0" to white
 if("0" %in% all_classes) {
   color_values["0"] = "#ffa0a060"
 }
 
 # Apply manual color scale
 g = g + ggplot2::scale_fill_manual(values = color_values)
 
 # Add observation dots if observations are provided
 if(!is.null(observations))
 {
   # Convert observation time stamps based on time_to_date setting
   if(time_to_date) {
     obs_dates = as.Date(observations$Tvec, origin=origin)
   } else {
     obs_dates = observations$Tvec
   }
   
   # Create data frame for observations
   obs_df = data.frame(date = obs_dates)
   
   # If Yvec is provided, use it for coloring; otherwise use a single color
   if(!is.null(observations$Yvec))
   {
     obs_df$message = as.character(observations$Yvec)
     # Add points colored by message type, positioned above truth rectangles
     g = g + ggplot2::geom_point(data = obs_df,
                                 ggplot2::aes(x = date, y = 1.15, color = message),
                                 size = 0.5,
                                 inherit.aes = FALSE)
   } else
   {
     # Add points with single color, positioned above truth rectangles
     g = g + ggplot2::geom_point(data = obs_df,
                                   ggplot2::aes(x = date, y = 1.15),
                                   color = "black",
                                   size = 0.5,
                                   inherit.aes = FALSE)
   }
   
   # Extend y-axis to show observation dots (above truth rectangles at 1.1)
   g = g + ggplot2::coord_cartesian(ylim = c(0, 1.2))
 } else if(!is.null(truth))
 {
   # If only truth is provided, extend y-axis slightly to show truth rectangles
   g = g + ggplot2::coord_cartesian(ylim = c(0, 1.15))
 }
 
 g
}
