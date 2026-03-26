#' @title plot.RJSMC:  Plotting results from RJSMC::SMC routine
#' @description Plots an area plot over time with the areas corresponding to probabilities for the different discrete events. If 'truth' is available, also a piecewise horizontal line with levels equal to the state number is displayed within the same plot (with levels described on the right horizontal axis). If 'observations' is provided, dots for each observation are displayed above the truth rectangles.
#' @param obj Object of class RJSMC, the output from RJSMC::SMC
#' @param origin Starting date for data (used only when time_to_date=TRUE)
#' @param pl Equal to V, Z, Q, F, or B depending on which variable to plot (B for breakpoints)
#' @param truth A list containing 'V', a vector with true states within each segment and 'B', a vector containing the starting points for each segment
#' @param clnames A vector giving names on the different classes. If NULL (default), classes are named 0,1,...
#' @param observations A list containing 'Tvec' (time stamps of observations) and optionally 'Yvec' (message labels for coloring). If NULL (default), no observation dots are displayed.
#' @param time_to_date Logical. If TRUE (default), converts time values to Date format for better visualization. If FALSE, displays time in its original numeric format (useful for debugging).
#' @param title Optional character string; if non-NULL, displayed as the plot title (ggplot2 \code{labs(title = ...)}).
#' @export

plot.RJSMC = function(obj,
                      origin=as.Date("2022-01-01"),
                      pl="V",
                     truth=NULL,clnames=NULL,
                     observations=NULL,
                     time_to_date=TRUE,
                     title=NULL)
{
 post_V = obj@posteriors_container_V
 post_Z = obj@posteriors_container_Z
 post_Q = obj@posteriors_container_Q
 post_F = obj@posteriors_container_F
 post_B = obj@posteriors_container_B
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
 else if (pl=="B")
 {
   # For breakpoints, create a simple line plot with probabilities
   df = data.frame(date=points, Probability=post_B)
 }

 # Handle breakpoint plot differently from state variable plots
 if(pl=="B")
 {
   # For breakpoints, create a line plot
   g = ggplot2::ggplot(data = df,
                       ggplot2::aes(x = date, y = Probability)) + 
       ggplot2::geom_line(size = 1.5, color = "blue") +
       ggplot2::ylab("Probability of Breakpoint") +
       ggplot2::ylim(0, max(1, max(df$Probability, na.rm=TRUE) * 1.1))
 } else
 {
   # For state variables (V, Z, Q, F), use area plots
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
 
   # Y-axis label for state-variable probability plots
   g = g + ggplot2::ylab("Probability")
 }
 
 # Set x-axis label based on time format
 if(time_to_date) {
   g = g + ggplot2::xlab("Date")
 } else {
   g = g + ggplot2::xlab("Time")
 }

 # Get all unique classes from the main plot (only for state variable plots)
 if(pl != "B") {
   all_classes = sort(unique(as.character(g.gath$Class)))
 } else {
   all_classes = character(0)  # Empty for breakpoint plots
 }
 
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
    
    # For breakpoint plots, add vertical lines at true breakpoints
    if(pl == "B") {
      # Add vertical lines at true breakpoint locations
      truth_breakpoints = if(time_to_date) as.Date(truth$B, origin=origin) else truth$B
      g = g + ggplot2::geom_vline(xintercept = truth_breakpoints, 
                                   linetype = "dashed", 
                                   color = "red", 
                                   size = 1.0,
                                   alpha = 0.7)
    } else {
      # For state variable plots, use rectangles as before
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
      
      # Clear separator line between estimated states (y 0--1) and real V-state band
      sep_y = 1.02
      g = g + ggplot2::geom_hline(yintercept = sep_y, linewidth = 0.8, colour = "gray30", linetype = "solid")
      
      # Real V-state band: place rectangles in a dedicated band above the separator
      truth_ymin = 1.05
      truth_ymax = 1.18
      g = g + ggplot2::geom_rect(data=g2, 
                                  ggplot2::aes(xmin=date, xmax=dat2, ymin=truth_ymin, ymax=truth_ymax, fill=Class),
                                  inherit.aes = FALSE)
      
      # Label the top band as "Real V-state" (legend-style at top)
      g = g + ggplot2::annotate("text", x = Inf, y = (truth_ymin + truth_ymax) / 2,
                                 label = "Real V-state", hjust = -0.08, vjust = 0.5,
                                 size = 3.2, fontface = "bold")

      # Add left-side label for the top band: "True State" (positioned where true states are displayed)
      # Match style of the main y-axis title
      # Place a large, clearly visible \"True State\" label near the left side,
      # aligned vertically with the true-state band. We overshoot the size/offset
      # now so it is definitely visible and can be fine-tuned later.
      x_min <- min(g2$date)
      x_range <- max(g2$date) - x_min
      # Put label well inside the plotting area (10% of the x-range from the left)
      x_pos <- x_min - 0.03 * x_range
      g = g +
        ggplot2::annotate(
          "text",
          x = x_pos,
          y = (truth_ymin + truth_ymax) / 2,
          label = "True State",
          angle = 90,
          hjust = 0.35,
          vjust = 0.23,
          size = 3.2,        # deliberately large for visibility
        )
      
      # Update all_classes to include truth classes
      all_classes = sort(unique(c(all_classes, as.character(g2$Class))))
    }
    }
 
 # Create color mapping: class "0" is always white, others get default colors
 # Only apply for state variable plots (not breakpoint plots)
 if(pl != "B") {
   # Generate colors for all classes except "0"
   n_classes = length(all_classes)
   color_values = grDevices::rainbow(n_classes)
   names(color_values) = all_classes

   # Set class "0" to very light gray (slightly brighter)
   if("0" %in% all_classes) {
     color_values["0"] = "#F2F2F2"
   }

   # Apply manual color scale; legend title is "V-State"
   g = g + ggplot2::scale_fill_manual(values = color_values, name = "V-State")
 }
 
 # Add observation dots if observations are provided
 if(!is.null(observations))
 {
   # Convert observation time stamps based on time_to_date setting
   if(time_to_date) {
     obs_dates = as.Date(observations$Tvec, origin=origin)
   } else {
     obs_dates = observations$Tvec
   }
   
   # Create data frame for observations (include legend label for legend entry)
   obs_df = data.frame(date = obs_dates, legend_obs = "Observations")
   
   if(pl == "B") {
     # For breakpoint plots, show observations at the top of the probability scale
     # Position them at y = max probability + small offset
     max_prob = max(df$Probability, na.rm=TRUE)
     obs_y = max(1, max_prob * 1.1)
     
     # If Yvec is provided, use it for coloring; otherwise use a single color
     if(!is.null(observations$Yvec))
     {
       obs_df$message = as.character(observations$Yvec)
       g = g + ggplot2::geom_point(data = obs_df,
                                   ggplot2::aes(x = date, y = obs_y, color = message),
                                   size = 0.5,
                                   inherit.aes = FALSE)
     } else
     {
       g = g + ggplot2::geom_point(data = obs_df,
                                     ggplot2::aes(x = date, y = obs_y, color = legend_obs),
                                     size = 0.5,
                                     inherit.aes = FALSE) +
         ggplot2::scale_color_manual(values = c("Observations" = "black"), name = NULL)
     }
     # Adjust y-axis to show observations
     g = g + ggplot2::coord_cartesian(ylim = c(0, obs_y * 1.1))
   } else {
     # For state variable plots, use original positioning
     # If Yvec is provided, use it for coloring; otherwise use a single color
    obs_y_state = 1.25  # above Real V-state band (1.05--1.18)
    if(!is.null(observations$Yvec))
    {
       obs_df$message = as.character(observations$Yvec)
       # Add points colored by message type, positioned above Real V-state band
       g = g + ggplot2::geom_point(data = obs_df,
                                   ggplot2::aes(x = date, y = obs_y_state, color = message),
                                   size = 0.5,
                                   inherit.aes = FALSE)
     } else
     {
       # Add points with single color (black), positioned above Real V-state band; show in legend
       g = g + ggplot2::geom_point(data = obs_df,
                                     ggplot2::aes(x = date, y = obs_y_state, color = legend_obs),
                                     size = 0.5,
                                     inherit.aes = FALSE) +
         ggplot2::scale_color_manual(values = c("Observations" = "black"), name = NULL)
     }
     
     # Extend y-axis to show observation dots above the Real V-state band
     g = g + ggplot2::coord_cartesian(ylim = c(0, 1.32))
   }
 } else if(!is.null(truth) && pl != "B")
 {
   # If only truth is provided, extend y-axis to show separator and Real V-state band
   g = g + ggplot2::coord_cartesian(ylim = c(0, 1.22))
 }

 if (!is.null(title)) {
   g = g + ggplot2::labs(title = title)
 }

 g
}
