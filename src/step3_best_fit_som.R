source('src/utils/library.R')

# Read in scaled cumulative timeseries data
cumulative_timeseries_scaled <- as.matrix(
  fread(
    file='output/python_output/full_factorial_cumulative_timeseries_scaled.csv',
    header=FALSE)
)

# Read in hyperparameter df
params_metrics_df <- read.csv('output/r_output/som_params_metrics.csv')

# Best parameters
idx = 442  # best som hyperparmeters chosen in step 2
n_epochs = 25  # best number of training epochs chosen in step 2

# Some repeat code from step 2. Create the neuron initialization matrix.
cov <- cov(cumulative_timeseries_scaled)
eigenvectors <- eigen(cov)
rm <- eigenvectors$vectors[, 1:2] # rotation matrix, eigenvectors of first two PCs
pcs <- cumulative_timeseries_scaled %*% rm

pc1range=range(pcs[,1]) # will use ranges to sample uniformly along PC1 and PC2
pc2range=range(pcs[,2])

d1 <- seq(
  from=pc1range[1],
  to=pc2range[2],
  length.out = params_metrics_df$x_dim[idx]
)

d2 <- seq(
  from=pc2range[1],
  to=pc2range[2],
  length.out=params_metrics_df$y_dim[idx]
)

pc_grid <- expand.grid(d1, d2)
init=as.matrix(pc_grid) %*% t(rm)

# Fit the best-fit SOM
best_som <- som(
  X=cumulative_timeseries_scaled,
  radius=params_metrics_df$init_radius[idx],
  dist.fcts=params_metrics_df$distance_fnc[idx],
  grid=somgrid(
    xdim=params_metrics_df$x_dim[idx],
    ydim=params_metrics_df$y_dim[idx],
    topo='hexagonal',
    toroidal=FALSE,
    neighbourhood.fct=params_metrics_df$neighborhood_fnc[idx]
  ),
  rlen=n_epochs,
  keep.data=TRUE,
  init=init,
  mode='pbatch'
)

# Save best_som object
saveRDS(object=best_som, file='output/r_output/best_som.rds')


# WRITE OUT BEST_SOM FIT INFORMATION FOR USE IN PYTHON

## 'Codes' are the weights of each neuron
som_codes <- best_som$codes[[1]]
colnames(som_codes) <- seq(2027, 2056)
write.csv(
  x=som_codes,
  file='output/r_output/som_codes_scaled.csv',
  row.names=FALSE
)

## Assignment of each SOW to a SOM neuron
neuron_ids <- matrix(
  data=1:nrow(cumulative_timeseries_scaled), 
  ncol=1
)
neuron_ids <- cbind(neuron_ids, best_som$unit.classif)

colnames(neuron_ids) <- c('SOW', 'Neuron')

write.csv(
  x=neuron_ids,
  file='output/r_output/som_sow_neuron_ids.csv',
  row.names=FALSE
)

## Coordinates of neurons
neuron_coordinates <- data.frame(best_som$grid$pts)

write.csv(
  x=neuron_coordinates,
  file='output/r_output/som_neuron_coordinates.csv',
  row.names=FALSE
)

## Rewrite the full_factorial_sow_info file to include neuron ids
full_factorial_sow_info <- fread(
  file='output/python_output/full_factorial_sow_info.csv',
  header=TRUE
)
full_factorial_sow_info$Neuron <- best_som$unit.classif

write.csv(
  x=full_factorial_sow_info,
  file='output/python_output/full_factorial_sow_info.csv',
  row.names=FALSE
)

