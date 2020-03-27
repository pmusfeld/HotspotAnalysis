\name{hotspot_analysis}
\alias{hotspot_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
function to find areas of high density (so called "Hotspots") in two-dimensional data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function takes a two-dimensional dataset as it's input and and detetcts areas of high density in the data. The distribution of the input data is compared to randomly distirbuted data of the same structure. Therefore, for every data-point the distance to a specified proportion of other nearest neighbor data-points is calculated and compared to the randomly expected value. Points which are closer to each other than randomly expected are selected and clustered to Hotpots.
}
\usage{
hotspot_analysis(point_id, column_x_coord, column_y_coord, prop_nneighbours = 0.05, criterion_NNI = 0.8, number_cluster = NULL, prefer_more_cluster = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{point_id}{
%%     ~~Describe \code{point_id} here~~
Input column which contains a point identifier for every point (x/y) in the data
}
  \item{column_x_coord}{
%%     ~~Describe \code{column_x_coord} here~~
Input column (a numeric vector) which contains the x-coordinates of the data
}
  \item{column_y_coord}{
%%     ~~Describe \code{column_y_coord} here~~
Input column (a numeric vector) which contains the y-coordinates of the data
}
  \item{prop_nneighbours}{
%%     ~~Describe \code{prop_nneighbours} here~~
A number < 1, specifying the proportion (percentage from the total number of data points) of nearest neighbors that should be considered when calculating the distances from every point to his sourrounded points. The default value is 5 percent (0.05) from the total number of points.
}
  \item{criterion_NNI}{
%%     ~~Describe \code{criterion_NNI} here~~
A number < 1, specifying the threshold for selecting points of high density. A Value of 1 indicates that a point is distributed as randomly expected. The closer the value is to 0, the stricter the selection of points. The default value is 0.8.
}
  \item{number_cluster}{
%%     ~~Describe \code{number_cluster} here~~
A number, manually specifying the number of cluster. By default, the number of cluster/hotspots is determined automatically by optimizing the silhouette width of different cluster solutions. Nevertheless, if the solution is not appropriate for personal needs or content-related issues, the number of needed Hotspots ca be specified manually.
}
  \item{prefer_more_cluster}{
%%     ~~Describe \code{prefer_more_cluster} here~~
A logical (TRUE/FALSE). By default, the cluster solution with the best fit (by optimizing silhouette width) is choosen. With regard to content-related issues, it sometimes can be more appropriate to prefer a higher number of smaller and close by cluster over less but bigger cluster. Therefore, all good cluster solutions in a small range from the optimum are selected and the solution with the highest number of cluster is chosen.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list object, containing a data frame with the Hotspot assignment and 3 graphs to control the selection process:
  \item{hotspot_assignment}{
  A data frame with 2 columns, containing the point identifier and the assigned Hotspot
  }
  \item{plot_result}{
  A scatter-plot of the data, showing the hotspot-assignment
  }
  \item{plot_point_selection}{
  A scatter plot, showing just the points which where selected to be distributed in areas of high density
  }
  \item{plot_random_data}{
  A scatter plot, showing the randomly distributed points which the data is compared to
  }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Philipp Musfeld
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##create data with Hotspots:
x_0_100_H1 <- rnorm(n = 100, mean = 20, sd = 5)
y_0_100_H1 <- rnorm(n = 100, mean = 20, sd = 5)

x_0_100_H2 <- rnorm(n = 100, mean = 50, sd = 5)
y_0_100_H2 <- rnorm(n = 100, mean = 50, sd = 5)

x_0_100_H3 <- rnorm(n = 100, mean = 90, sd = 5)
y_0_100_H3 <- rnorm(n = 100, mean = 90, sd = 5)

x_0_100_H4 <- rnorm(n = 50, mean = 25, sd = 2.5)
y_0_100_H4 <- rnorm(n = 50, mean = 95, sd = 2.5)

x_0_100_H5 <- rnorm(n = 50, mean = 70, sd = 5)
y_0_100_H5 <- rnorm(n = 50, mean = 30, sd = 5)

x_0_100_kH <- runif(n = 200, min = 0, max = 100)
y_0_100_kH <- runif(n = 200, min = 0, max = 100)

x_coord <- c(x_0_100_H1, x_0_100_H2,  x_0_100_H3,  x_0_100_H4,  x_0_100_H5,  x_0_100_kH)
y_coord <- c(y_0_100_H1, y_0_100_H2,  y_0_100_H3,  y_0_100_H4,  y_0_100_H5,  y_0_100_kH)

data <- data.frame(x_coord, y_coord)
data$point <- paste("point_", rownames(data), sep = "")

hotspot_analysis(point_id = data$point, column_x_coord = data$x_coord, column_y_coord = data$y_coord, prop_nneighbours = 0.05, criterion_NNI = 0.8, number_cluster = NULL, prefer_more_cluster = TRUE)
}


