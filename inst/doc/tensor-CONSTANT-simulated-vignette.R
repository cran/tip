## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tip)

# ################## NOTE ##################
# Order 3 Tensor dimension: c(d1,d2,d3)
# d1: number of rows
# d2: number of columns
# d3: number of slices 
############################################

# Set a random seed for reproducibility 
set.seed(007)

# A function to generate an order-3 tensor
generate_gaussian_tensor <- function(.tensor_dimension, .mean = 0, .sd = 1){
  array(data = c(rnorm(n = prod(.tensor_dimension), 
                       mean = .mean, 
                       sd = .sd)), 
        dim = .tensor_dimension)
} 

# Define the tensor dimension
tensor_dimension <- c(256,256,3)

# Generate clusters of tensors
c1 <- lapply(1:10, function(x) generate_gaussian_tensor(.tensor_dimension = tensor_dimension,
                                                        .mean = 0,
                                                        .sd = 1))


# Generate clusters of tensors
c2 <- lapply(1:10, function(x) generate_gaussian_tensor(.tensor_dimension = tensor_dimension,
                                                        .mean = 5,
                                                        .sd = 1))


# Generate clusters of tensors
c3 <- lapply(1:10, function(x) generate_gaussian_tensor(.tensor_dimension = tensor_dimension,
                                                        .mean = -5,
                                                        .sd = 1))

# Make a list of tensors 
X <- c(c1, c2, c3)

# Compute the number of subjects for each cluster 
n1 <- length(c1)
n2 <- length(c2)
n3 <- length(c3)

# Create a vector of true labels. True labels are only necessary 
# for constructing network graphs that incorporate the true labels;
# this is often useful for research. 
true_labels <- c(rep("Cluster 1", n1),
                 rep("Cluster 2", n2),
                 rep("Cluster 3", n3))

# Compute the total number of subjects 
n <- length(X)

# Construct the distance matrix 
distance_matrix <- matrix(data = NA, nrow = n, ncol = n) 
for(i in 1:n){
  for(j in i:n){
    distance_matrix[i,j] <- sqrt(sum((X[[i]] - X[[j]])^2))
    distance_matrix[j,i] <- distance_matrix[i,j]
  }
}

# Compute the temperature parameter estiamte
temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])

# For each subject, compute the point estimate for the number of similar 
# subjects using  univariate multiple change point detection (i.e.)
init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)

# Set the number of burn-in iterations in the Gibbs samlper
# RECOMMENDATION: burn >= 1000
burn <- 1000

# Set the number of sampling iterations in the Gibbs sampler
# RECOMMENDATION: samples >= 1000
samples <- 1000

# Set the subject names
names_subjects <- paste(1:n)

# Run TIP clustering using only the prior
# --> That is, the likelihood function is constant
tip1 <- tip(.data = list(),
            .burn = burn,
            .samples = samples,
            .similarity_matrix = exp(-1.0*temperature*distance_matrix),
            .init_num_neighbors = init_num_neighbors,
            .likelihood_model = "CONSTANT",
            .subject_names = names_subjects,
            .num_cores = 1)

## -----------------------------------------------------------------------------
# Produce plots for the Bayesian Clustering Model
tip_plots <- plot(tip1)

## -----------------------------------------------------------------------------
# View the posterior distribution of the number of clusters
tip_plots$histogram_posterior_number_of_clusters

## -----------------------------------------------------------------------------
# View the trace plot with respect to the posterior number of clusters
tip_plots$trace_plot_posterior_number_of_clusters

## -----------------------------------------------------------------------------
# Extract posterior cluster assignments using the Posterior Expected Adjusted Rand (PEAR) index
cluster_assignments <- mcclust::maxpear(psm = tip1@posterior_similarity_matrix)$cl

# If the true labels are available, then show the cluster result via a contigency table
table(data.frame(true_label = true_labels,
                 cluster_assignment = cluster_assignments))

## -----------------------------------------------------------------------------
# Create the one component graph with minimum entropy
partition_list <- partition_undirected_graph(.graph_matrix = tip1@posterior_similarity_matrix,
                                             .num_components = 1,
                                             .step_size = 0.001)

## -----------------------------------------------------------------------------
# Associate class labels and colors for the plot
class_palette_colors <- c("Cluster 1" = "blue",
                          "Cluster 2" = 'green',
                          "Cluster 3" = "red")

# Associate class labels and shapes for the plot
class_palette_shapes <- c("Cluster 1" = 19,
                          "Cluster 2" = 18,
                          "Cluster 3" = 17)

# Visualize the posterior similarity matrix by constructing a graph plot of 
# the one-cluster graph. The true labels are used here (below they are not).
ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                .subject_names = NA,
                .subject_class_names = true_labels,
                .class_colors = class_palette_colors,
                .class_shapes = class_palette_shapes,
                .node_size = 2,
                .add_node_labels = FALSE)

## -----------------------------------------------------------------------------
# If true labels are not available, then construct a network plot
# of the one-cluster graph without any class labels.
# Note: Subject labels may be suppressed using .add_node_labels = FALSE.  
ggnet2_network_plot(.matrix_graph = partition_list$partitioned_graph_matrix,
                .subject_names = names_subjects,
                .node_size = 2,
                .add_node_labels = TRUE)

