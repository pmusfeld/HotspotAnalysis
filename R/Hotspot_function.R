
hotspot_analysis <- function(point_id, column_x_coord, column_y_coord, prop_nneighbours = 0.05, criterion_NNI = 0.8, number_cluster = NULL, prefer_more_cluster = TRUE){


  #create function to calculate the distance to a specified number of nearest neighbours
  distance_to_n_nearest_neighbours <- function(x_coord, y_coord, n_neighbours){
    distance_list <- list()
    for(i in c(1:n_neighbours)){
      distance <- spatstat::nndist(X = x_coord,  Y = y_coord, k = i)
      distance_list[[i]] <- distance
    }
    return(dplyr::bind_cols(distance_list))
  }


  #create data.frame from input_columns
  data <- data.frame(point = point_id, x_coord = column_x_coord, y_coord = column_y_coord)


  #calculate distance from each point to its n nearest neighbours. n is set so 10% of the total number of points
  distances <- data.frame(point = data$point, distance_to_n_nearest_neighbours(x_coord = data$x_coord, y_coord = data$y_coord, n_neighbours = floor(length(data$point)*prop_nneighbours)))


  #for each point, calculate the sum of all distances to the n nearest neighbours and devide by total number of points
  distances_agg <- distances %>%
    tidyr::pivot_longer(cols = -point, names_to = "names", values_to = "distance") %>%
    dplyr::group_by(point) %>%
    dplyr::summarize(sum_distance = sum(distance)) %>%
    dplyr::mutate(rel_distance = sum_distance / length(distances$point)) %>%
    dplyr::ungroup()



  #simulate data with the same length as the analysed data but randomly distributed
  simulate_x_coord <- stats::runif(n = length(data$point), min = min(data$x_coord), max = max(data$x_coord))
  simulate_y_coord <- stats::runif(n = length(data$point), min = min(data$y_coord), max = max(data$y_coord))


  #create plot putput for random data
  plot_random_data <- ggplot2::ggplot(data = NULL, ggplot2::aes(x = simulate_x_coord, y = simulate_y_coord)) + ggplot2::geom_point(alpha = 0.7)


  #calculate distances for random points
  random <- data.frame(distance_to_n_nearest_neighbours(x_coord = simulate_x_coord, y_coord = simulate_y_coord, n_neighbours = floor(length(data$point)*prop_nneighbours)))
  random$point = paste("point_", rownames(random), sep = "")


  #for each random point calculate the sum of all of his n distances and devide by total number of points
  random_agg <- random %>%
    tidyr::pivot_longer(cols = -point, names_to = "names", values_to = "distance_random") %>%
    dplyr::group_by(point) %>%
    dplyr::summarize(sum_distance_random = sum(distance_random)) %>%
    dplyr::mutate(rel_distance_random = sum_distance_random / length(distances$point)) %>%
    dplyr::ungroup() %>%
    #take the mean of the relative sums. This value serves as the expected value if all points were randomly distributed
    dplyr::mutate(mean_distance_random = mean(rel_distance_random))


  #for each point, calculate a nearest neighbour index by dividing the observed relative distance by the expected relative distance for randomly distributed points. A score of 1 indicates, that the observed point has the same distance to other points, as a randomly distributed point. A score below 1 indicates that this point has a lower distance to other points than expected and therefore is likely to be part of a hotspot area
  distances_NNI <- distances_agg %>%
    cbind(mean_distance_random = random_agg$mean_distance_random) %>%
    dplyr::mutate(NNI = rel_distance/mean_distance_random) %>%
    #select just points with a nearest neighbour index below the specified threshold
    dplyr::filter(NNI < criterion_NNI) %>%
    dplyr::select(point, NNI)


  #match coordinate information to selected points
  points_select <- distances_NNI %>%
    dplyr::inner_join(dplyr::select(data, point, x_coord, y_coord), by = "point")


  #save selected points (Hotpsotpoints) for plot-output
  plot_point_selection <- ggplot2::ggplot(points_select, ggplot2::aes(x=x_coord, y=y_coord)) + ggplot2::geom_point(alpha = 0.7)


  #assign selected hotspot points to clusters - method used: hierarchical clustering (starts by treating every point as a cluster and step by step merges the clusters with the smallest distance together. Distance measure is euclidean and linkage is average. Average has a tendency to form clusters with the same and small variance which could be usefull for our use case)

  cluster_analysis <- stats::hclust(dist(points_select[, c(3,4)]), method = "average")


  #if number of cluster is not specified: find optimal cluster solution v<
  ##calculating the silhouette width for a various number of cluster solutions. This value is an fit index based on compactness and seperation of the clusters. For each point i of a given cluster C, the average distance to the other points of the same cluster is calculated. Next, the average distance to all points of the nearest cluster Cn is calculated and both distance measures are related to each other. A silhouette index close to 1 implies a good fit of the point to its own cluster. A silhouette index close to -1 implies a good fit of the point to the "neighbor" cluster. The Average of all silhouette indexes can be used as an indicator for the fit of the cluster solution. This index favours solutions with more but compact clusters over solutions with less but scatter clusters and should therefore be appropriate for our use case.


  if(is.null(number_cluster)){

    silhouette_width <- purrr::map_dbl(2:20,  function(.x){
      cluster_solution <- stats::cutree(cluster_analysis, .x)
      solution_fit <- cluster::silhouette(cluster_solution, dist(points_select[, c(3,4)]), full = TRUE)
      mean(solution_fit[, 3])
    })


    #if prefer_more_cluster = TRUE, select best solutions in a specific range of silhouette width score and choose solution with the highest number of cluster
    if(prefer_more_cluster == TRUE){

      number_of_cluster <- data.frame(number_of_cluster = c(2:20), silhouette_width = silhouette_width) %>%
        dplyr::arrange(desc(silhouette_width)) %>%
        dplyr::filter(silhouette_width >= max(silhouette_width)-0.03) %>%
        dplyr::filter(number_of_cluster == max(number_of_cluster))
    }


    #if prefer_more_cluster = FALSE, choose the soluation with the max silhouette width score / the best fit
    else{

      number_of_cluster <- data.frame(number_of_cluster = c(2:20), silhouette_width = silhouette_width) %>%
        dplyr::filter(silhouette_width == max(silhouette_width))
    }


    ##clustering the data with the calculated number of cluster
    cluster <- stats::cutree(cluster_analysis, k = number_of_cluster$number_of_cluster)
  }


  #if number of cluster is specified: use given number of cluster for clustering
  else{
    cluster <- stats::cutree(cluster_analysis, k = number_cluster)
  }



  ##assign cluster to points
  points_select$cluster <- as.factor(cluster)


  ##combine original data with clustered data
  data_hotspots <- data %>%
    dplyr::left_join(dplyr::select(points_select, point, cluster), by = "point") %>%
    dplyr::mutate(Hotspot = as.factor(ifelse(is.na(cluster), "kein Hotspot", paste("Hotspot ", cluster))))


  ##visualise results
  plot_result <- ggplot2::ggplot(data_hotspots, ggplot2::aes(x = x_coord, y = y_coord, col = Hotspot)) +
    ggplot2::geom_point(alpha = 0.7)


  ##select point id and assigned cluster
  hotspot_assignment <- dplyr::select(data_hotspots, point, Hotspot)


  #save plot for selected points, plot for results, plot for random distribution and clustering result in output list
  output <- list(hotspot_assignment = hotspot_assignment, plot_point_selection = plot_point_selection, plot_result = plot_result, plot_random_data = plot_random_data)

  return(output)
}
