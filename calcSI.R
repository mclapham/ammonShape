#READ FILES, PERFORM DIGITIZATION AND CALCULATION FOR CEPHALOPOD SUTURES
#default is to perform Dijkstra's shortest path algorithm
#can specify dig_method="spline" for cubic spline fit
#Dijkstra method has moving window: window_edge sets starting point
#window_inc sets increment of moving window (in pixels); can be reduced if memory issues encountered
#spline method has smoothing factor, which can be varied with smooth= option
#parameter of 0.3-0.5 is reasonable for most, but you can use trial and error
#filled_line option (set to T) allows digitization of sutures where one half is filled with solid black
calc.si <- function(file_name, dig_method="dijkstra", window_edge=50, window_inc=50, smooth=0.4, filled_line=F) {
  
  #read file from disk and convert to data frame with x, y coords and pixel intensity
  suture_image <- readbitmap::read.bitmap(paste(file_name, ".png", sep=""))
  suture <- data.frame(x=sort(rep(seq(ncol(suture_image)),nrow(suture_image))), y=rep(seq(nrow(suture_image)), ncol(suture_image)), pixel=c(suture_image))
  
  if (filled_line==T) {
    blackpts <- subset(suture, suture$pixel==0) #extract black pixels
    whitepts <- subset(suture, suture$pixel==1) #extract white pixels
    
    #calculate distance from each black pixel to the nearest white pixel
    blackpts$dists <- apply(blackpts, 1, function(x) {
      sqrt(min((x[1]-whitepts[,1])^2 + (x[2]-whitepts[,2])^2))
    })
    
    #selects only black pixel that touch a white pixel
    suture <- subset(blackpts, blackpts$dists==1)
  }

  suture_sub <- subset(suture[1:2], suture$pixel == 0) #selects only darkest pixels
  
  plot(suture_sub$x, suture_sub$y, pch=15, cex=0.5, ylim=rev(range(suture_sub$y))) #plots raw suture pixels    
    
  if (dig_method=="spline") {

    suture.ma <- smooth.spline(suture_sub$x, suture_sub$y, spar = smooth) #cubic spline smoothing to generate line
    
    suture_pts <<- data.frame(x=suture.ma$x, y=suture.ma$y) #extracts x/y points from spline fit
    
    points(suture_pts, col=rainbow(nrow(suture_pts))) #plots points as quality control check
    
  } else if (dig_method=="dijkstra") {
    
    path_x <- numeric(0) #creates empty variable to store suture trace
    path_y <- numeric(0)
    
    startx <- min(suture_sub$x) #finds starting point for new segment (=end point for previous segment)
    starty <- min(suture_sub$y[which(suture_sub$x==startx)])
    
    #increments window to select parts of suture for digitization
    while (window_edge < max(suture_sub$x)+(window_inc-1)) {
      
      #sets endpoint of final window increment to actual end of suture
      #rather than point beyond end
      #endpoint needs to be black pixel
      if (window_edge > max(suture_sub$x)) {
        window_edge <- max(suture_sub$x)+1
      }
      
      #selects part of suture to the left of the window edge for analysis
      suture_part <- subset(suture_sub, suture_sub$x < window_edge)
      
      #finds start point
      start_pt <- which(suture_part$x==startx & suture_part$y==starty)
      
      pt_pairs <- data.frame(t(combn(seq(nrow(suture_part)), 2))) #creates list of all possible point pairs (vertex points of graph)
      
      names(pt_pairs) <- c("v1", "v2") #not sure if this is necessary
      
      pt_pairs$weight <- as.numeric(dist(suture_part)) #calculates Euclidean distance btwn pairs
      pt_pairs$weight[pt_pairs$weight>1.5] <- 100000 #sets distance weight to 10^5 for any distance 2 or greater
      #this forces the algorithm to only move between adjacent points (never jumping across gaps)
      
      ptgraph <- igraph::graph.data.frame(pt_pairs, directed=F) #creates graph network
      
      #there are multiple possible endpoints at a given x coordinate cutoff
      #the algorithm needs to try each one to find one that is actually connected to line segment
      for (i in which(suture_part$x == (window_edge-1))) {
        #uses Dijkstra's algorithm with weights to find shortest path between vertices
        sh_path <- igraph::get.shortest.paths(ptgraph, from=start_pt, to=i, weights=igraph::E(ptgraph)$weight)
        
        pt_order <- unlist(sh_path$vpath)
        
        #if the path is unconnected there will be only two points (the start and end)
        #otherwise it has found a path that works
        #it will add the x and y coordinates to the running list and stop looking
        if (length(pt_order) > 2) {
          
          path_x <- c(path_x, suture_part$x[pt_order]) #adds pixels corresponding to shortest path
          path_y <- c(path_y, suture_part$y[pt_order])
          
          #reset start point to end point of previous segment
          startx <- path_x[length(path_x)]
          starty <- path_y[length(path_y)]
          
          #creates a record of the points that have already been visited
          visited_points <- data.frame(x=path_x[-length(path_x)], y=path_y[-length(path_y)])
          
          #calculates distance from each suture point to the nearest visited point
          suture_part$dists <- apply(suture_part, 1, function(x) {
            sqrt(min((x[1]-visited_points[,1])^2 + (x[2]-visited_points[,2])^2))
          })
          
          #calculates distance from each suture point to the last visited point (new start point)
          suture_part$final_dists <- apply(suture_part, 1, function(x) {
            sqrt(min((x[1]-path_x[length(path_x)])^2 + (x[2]-path_y[length(path_y)])^2))
          })
          
          #creates list of points adjacent to line but not in contact with the last visited point (new start point)
          #extra step (distance to last visited point) avoids problem where gap in line may be created by deletion
          adjacent_pts <- paste(subset(suture_part$x, suture_part$dists < 1.5 & suture_part$final_dists > 1.5), subset(suture_part$y, suture_part$dists < 1.5 & suture_part$final_dists > 1.5))
          
          #removes points adjacent to line from data 
          suture_sub <- subset(suture_sub, !paste(suture_sub$x, suture_sub$y) %in% adjacent_pts)
          
          #remove points from early enough that they won't be involved in the line
          #it removes any point that is more than 0.4 suture widths prior to the current position
          #this is simply to reduce the number of points include in the calculation for efficiency
          suture_sub <- subset(suture_sub, suture_sub$x > max(path_x) - 0.4*max(suture_sub$x))
          
          #removes visited points (except last visited)
          #second command adds last visited point back
          suture_sub <- subset(suture_sub, !paste(suture_sub$x, suture_sub$y) %in% paste(path_x, path_y))
          suture_sub <- rbind(suture_sub, c(startx, starty))
          
          break #breaks out of for loop once solution found
        } else if (i == max(which(suture_part$x == (window_edge-1)))) {
          stop("There is a gap in the suture line")
          #if it has tried all possible end points without finding a solution,
          #there must be a gap, so the function stops with this error message
        }
      }
      
      #adds current points to graph (as progress and quality control check)
      points(path_x, path_y, col=rainbow(length(path_x)), cex=0.5)
      
      #increments window of interest by specified amount
      window_edge <- window_edge + window_inc
      
    }
    
    #formats path into data frame of visited points
    suture_pts <<- data.frame(x=path_x, y=path_y)
    
    } else {
    stop("Invalid method name") #if the name is not spline or dijkstra
    }
  
  #calculate the suture width (in pixels)
  SW <- max(suture_pts$x) - min(suture_pts$x) 
  
  #finds start and end points of suture
  x.final <- suture_pts[nrow(suture_pts), 1]
  x.start <- suture_pts[1, 1]
  y.final <- suture_pts[nrow(suture_pts), 2]
  y.start <- suture_pts[1, 2]
  
  dist.straight <- sqrt((x.final - x.start)^2 + (y.final - y.start)^2) #calculates straight-line distance from start to end
  
  dist.curved <- sum(sqrt(apply(data.frame(diff(suture_pts[, 1])^2, diff(suture_pts[, 2])^2), 1, sum))) #calculates curved distance along suture perimeter
  
  SI <- dist.curved/dist.straight #calculates suture index (SI)
  
  #save the suture points as a file
  write.csv(suture_pts, paste(file_name, ".csv", sep=""))  

  suture_results <- list("SI"=SI, "SW"=SW)
  
  suture_results
}
