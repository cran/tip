\donttest{

  # Import the tip library
  library(tip)

  # Import the iris dataset
  data(iris)

  # The first 4 columns are the data whereas
  # the 5th column refers to the true labels
  X <- data.matrix(iris[,c("Sepal.Length",
                           "Sepal.Width",
                           "Petal.Length",
                           "Petal.Width")])

  # Extract the true labels (optional)
  # True labels are only necessary for constructing network
  # graphs that incorporate the true labels; this is often
  # for research.
  true_labels <- iris[,"Species"]

  # Compute the distance matrix
  distance_matrix <- data.matrix(dist(X))

  # Compute the temperature parameter estiamte
  temperature <- 1/median(distance_matrix[upper.tri(distance_matrix)])

  # For each subject, compute the point estimate for the number of similar
  # subjects using  univariate multiple change point detection (i.e.)
  init_num_neighbors = get_cpt_neighbors(.distance_matrix = distance_matrix)

  # Set the number of burn-in iterations in the Gibbs samlper
  # RECOMMENDATION: burn >= 1000
  burn <- 1000

  # Set the number of sampling iterations in the Gibbs sampler
  # RECOMMENDATION: burn >= 1000
  samples <- 1000

  # Set the subject names
  names_subjects <- paste(1:dim(iris)[1])

  # Run TIP clustering using only the prior
  # --> That is, the likelihood function is constant
  tip1 <- tip(.data = data.matrix(X),
              .burn = burn,
              .samples = samples,
              .similarity_matrix = exp(-1.0*temperature*distance_matrix),
              .init_num_neighbors = init_num_neighbors,
              .likelihood_model = "NIW",
              .subject_names = names_subjects,
              .num_cores = 1)

  # Produce plots for the Bayesian Clustering Model
  tip_plots <- plot(tip1)
}
