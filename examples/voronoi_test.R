library(Rcpp)
library(ggplot2)

# Source the compiled C++ code
#sourceCpp("src/voronoi.cpp")

# Define Voronoi input points
points <- matrix(c(1, 2,
                   3, 4,
                   5, 1,
                   7, 6,
                   2, 8,
                   3,5,
                   1, -4 ),#,
#                   0,-14,
#                   -4,-12,
#                   4,-12),
                 ncol = 2, byrow = TRUE)


x_min <- -10
x_max <- 10
y_min <- -10
y_max <- 10

# Compute Voronoi diagram
voronoi_result <- generate_voronoi(points, x_min, x_max, y_min, y_max)

voronoi_result$vertices
voronoi_result$polygons
voronoi_result$site_points
sum(voronoi_result$polygon_areas)
(voronoi_result$x_max - voronoi_result$x_min)*(voronoi_result$y_max - voronoi_result$y_min)
voronoi_result$polygon_areas

voronoi_result$site3

# Convert finite edges to a data frame
finite_edges_df <- do.call(rbind, voronoi_result$finite_edges)
colnames(finite_edges_df) <- c("x1", "y1", "x2", "y2")

# Convert infinite edges to a data frame
infinite_edges_df <- do.call(rbind, voronoi_result$infinite_edges)
colnames(infinite_edges_df) <- c("x1", "y1", "x2", "y2")

# Convert boundary edges to a data frame
boundary_edges_df <- do.call(rbind, voronoi_result$boundary_edges)
colnames(boundary_edges_df) <- c("x1", "y1", "x2", "y2")

# Convert points to a data frame
points_df <- as.data.frame(points)
colnames(points_df) <- c("x", "y")

# Plot Voronoi edges, boundary, and input points
ggplot() +
  # Plot finite edges in blue
  geom_segment(data = finite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "blue", linewidth = 1) +
  # Plot infinite edges in green
  geom_segment(data = infinite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "green", linewidth = 1) +
  # Plot boundary edges in black
  geom_segment(data = boundary_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "black", linewidth = 1) +
  # Plot input points in red
  geom_point(data = points_df, aes(x = x, y = y), color = "red", size = 3) +
  coord_fixed() +
  theme_bw() +
  ggtitle("Boost Voronoi Diagram in R with Boundaries")+
  scale_x_continuous(breaks = seq(voronoi_result$x_min, voronoi_result$x_max, by = 2)) +   # Smaller 'by' means denser vertical grid lines
  scale_y_continuous(breaks = seq(voronoi_result$y_min, voronoi_result$y_max, by = 2))    # Smaller 'by' means denser horizontal grid lines


site_points_matrix <- do.call(rbind, voronoi_result$site_points)
labels_df <- as.data.frame(site_points_matrix)
colnames(labels_df) <- c("x", "y")

# Combine site points and polygon areas for labeling
labels_df$area <- voronoi_result$polygon_areas

# Format area (e.g., round to 2 decimal places)
labels_df$label <- sprintf("%.1f", labels_df$area)

# Plot with labels
ggplot() +
  # Plot finite edges in blue
  geom_segment(data = finite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "blue", linewidth  = 1) +
  # Plot infinite edges in green
  geom_segment(data = infinite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "green", linewidth  = 1) +
  # Plot boundary edges in black
  geom_segment(data = boundary_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "black", linewidth  = 1) +
  # Plot input points in red
  geom_point(data = points_df, aes(x = x, y = y), color = "red", size = 3) +
  # Add area labels
  geom_text(data = labels_df, aes(x = x, y = y, label = label), color = "purple", size = 2, vjust = 1, hjust = -0.3) +
  coord_fixed() +
  theme_bw() +
  ggtitle("Boost Voronoi Diagram in R with Polygon Areas") +
  scale_x_continuous(breaks = seq(voronoi_result$x_min, voronoi_result$x_max, by = 2)) +   # Smaller 'by' means denser vertical grid lines
  scale_y_continuous(breaks = seq(voronoi_result$y_min, voronoi_result$y_max, by = 2))    # Smaller 'by' means denser horizontal grid lines


#################################################################################


library(Rcpp)
library(ggplot2)

# Source the compiled C++ code
#sourceCpp("src/voronoi.cpp")

# Define Voronoi input points
set.seed(7)
points <- matrix(unique(rnorm(100,mean = 0,4)),#,
                 #                   0,-14,
                 #                   -4,-12,
                 #                   4,-12),
                 ncol = 2, byrow = TRUE)


x_min <- -10
x_max <- 10
y_min <- -10
y_max <- 10

# Compute Voronoi diagram
voronoi_result <- generate_voronoi(points, x_min, x_max, y_min, y_max)

voronoi_result$vertices
voronoi_result$polygons
voronoi_result$site_points
sum(voronoi_result$polygon_areas)
(voronoi_result$x_max - voronoi_result$x_min)*(voronoi_result$y_max - voronoi_result$y_min)
voronoi_result$polygon_areas

voronoi_result$site3

# Convert finite edges to a data frame
finite_edges_df <- do.call(rbind, voronoi_result$finite_edges)
colnames(finite_edges_df) <- c("x1", "y1", "x2", "y2")

# Convert infinite edges to a data frame
infinite_edges_df <- do.call(rbind, voronoi_result$infinite_edges)
colnames(infinite_edges_df) <- c("x1", "y1", "x2", "y2")

# Convert boundary edges to a data frame
boundary_edges_df <- do.call(rbind, voronoi_result$boundary_edges)
colnames(boundary_edges_df) <- c("x1", "y1", "x2", "y2")

# Convert points to a data frame
points_df <- as.data.frame(points)
colnames(points_df) <- c("x", "y")

# Plot Voronoi edges, boundary, and input points
ggplot() +
  # Plot finite edges in blue
  geom_segment(data = finite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "blue", size = 1) +
  # Plot infinite edges in green
  geom_segment(data = infinite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "green", size = 1) +
  # Plot boundary edges in black
  geom_segment(data = boundary_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "black", size = 1) +
  # Plot input points in red
  geom_point(data = points_df, aes(x = x, y = y), color = "red", size = 3) +
  coord_fixed() +
  theme_bw() +
  ggtitle("Boost Voronoi Diagram in R with Boundaries")+
  scale_x_continuous(breaks = seq(voronoi_result$x_min, voronoi_result$x_max, by = (voronoi_result$x_max - voronoi_result$x_min)/10)) +   # Smaller 'by' means denser vertical grid lines
  scale_y_continuous(breaks = seq(voronoi_result$y_min, voronoi_result$y_max, by =(voronoi_result$x_max - voronoi_result$x_min)/10))    # Smaller 'by' means denser horizontal grid lines


site_points_matrix <- do.call(rbind, voronoi_result$site_points)
labels_df <- as.data.frame(site_points_matrix)
colnames(labels_df) <- c("x", "y")

# Combine site points and polygon areas for labeling
labels_df$area <- voronoi_result$polygon_areas

# Format area (e.g., round to 2 decimal places)
labels_df$label <- sprintf("%.1f", labels_df$area)

# Plot with labels
ggplot() +
  # Plot finite edges in blue
  geom_segment(data = finite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "blue", size = 1) +
  # Plot infinite edges in green
  geom_segment(data = infinite_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "green", size = 1) +
  # Plot boundary edges in black
  geom_segment(data = boundary_edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "black", size = 1) +
  # Plot input points in red
  geom_point(data = points_df, aes(x = x, y = y), color = "red", size = 3) +
  # Add area labels
  geom_text(data = labels_df, aes(x = x, y = y, label = label), color = "purple", size = 2, vjust = 1, hjust = -0.3) +
  coord_fixed() +
  theme_bw() +
  ggtitle("Boost Voronoi Diagram in R with Polygon Areas") +
  scale_x_continuous(breaks = seq(voronoi_result$x_min, voronoi_result$x_max, by = (voronoi_result$x_max - voronoi_result$x_min)/10)) +   # Smaller 'by' means denser vertical grid lines
  scale_y_continuous(breaks = seq(voronoi_result$y_min, voronoi_result$y_max, by =(voronoi_result$x_max - voronoi_result$x_min)/10))    # Smaller 'by' means denser horizontal grid lines
