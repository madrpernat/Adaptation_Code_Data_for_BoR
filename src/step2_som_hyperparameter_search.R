# SOM hyperparameter search

# Originally coded by Nathan Bonham - March 2022 
# (https://doi.org/10.1016/j.envsoft.2022.105491)

# Adapted by Madeline Pernat - November 2023

source('src/utils/clear.R')
source('src/utils/som_library.R')

# Read in scaled cumulative timeseries data
cumulative_timeseries_scaled <- as.matrix(
  fread(
    file='output/full_factorial_cumulative_timeseries_scaled.csv',
    header=FALSE)
)

# LATIN HYPERCUBE SAMPLING

set.seed(3)

n_samples <- 500

unit_cube <- improvedLHS(
  n=n_samples,
  k=5
)
unit_cube <- as.data.frame(unit_cube)

colnames(unit_cube) <- c('radius', 'y_dim', 'x_y_ratio', 'distance_fnc',
                         'neighborhood_fnc')

## Define ranges of parameters
radius_range <- c(0, 1)
y_dim_range <- c(2, 6)
x_y_ratio_range <- c(1, 12)
distance_fnc <- c('euclidean', 'manhattan', 'sumofsquares')
neighborhood_fnc <- c('bubble', 'gaussian')

## Convert LHS samples (which are between 0-1) to uniform sample across the
## given range for each CONTINUOUS parameter
cube <- unit_cube

cube$radius <- runif(
  unit_cube$radius,
  min=radius_range[1],
  max=radius_range[2]
)

cube$y_dim <- round(runif(
  unit_cube$y_dim,
  min=y_dim_range[1] - 0.5,
  max=y_dim_range[2] + 0.5
))

cube$x_y_ratio <- runif(
  unit_cube$x_y_ratio,
  min=x_y_ratio_range[1],
  max=x_y_ratio_range[2]
)

## Convert LHS samples (which are between 0-1) to uniform sample across the
## given range for each DISCRETE parameter
distance_probs <- seq(0, 1, length.out = length(distance_fnc) + 1)
neighborhood_probs <- seq(0, 1, length.out = length(neighborhood_fnc) + 1)

cube$distance_fnc <- continuous_to_discrete(
  cube$distance_fnc,
  probs=distance_probs,
  categories=distance_fnc
)

cube$neighborhood_fnc <- continuous_to_discrete(
  cube$neighborhood_fnc,
  probs=neighborhood_probs,
  categories=neighborhood_fnc
)

## Calculate X dimension and total # of neurons given Y dimension and X/Y ratio
cube$x_dim <- round(cube$y_dim * cube$x_y_ratio)
cube$n_neurons <- cube$y_dim * cube$x_dim

## Test that each hyperparameter value is equally represented in the hypercube
hist(cube$radius)
hist(cube$y_dim)
hist(cube$x_y_ratio, breaks=seq(0, 14))
hist(cube$n_neurons)
hist(cube$x_dim)
table(cube$distance_fnc)
table(cube$neighborhood_fnc)



# TRAIN SOMS AND REPORT QUALITY OF EACH CONFIG

## Pre-allocate matrix to store quality of fit metrics
metrics_df <- matrix(
  nrow=n_samples,
  ncol=4
)
metrics_df <- data.frame(metrics_df)

## Pre-allocate matrix to store parameters
params_df <- matrix(
  nrow=n_samples,
  ncol=8
)
params_df <- data.frame(params_df)
colnames(params_df) <- c('x_dim', 'y_dim', 'x_y_ratio', 'n_neurons',
                         'init_radius', 'distance_fnc', 'neighborhood_fnc',
                         'radius_fraction')


## Calculate rotation matrix and PCs of data. Will be used to initialize SOMs.
cov <- cov(cumulative_timeseries_scaled)
eigenvectors <- eigen(cov)
rm <- eigenvectors$vectors[, 1:2] # rotation matrix, eigenvectors of first two PCs
pcs <- cumulative_timeseries_scaled %*% rm

pc1range=range(pcs[,1]) # will use ranges to sample uniformly along PC1 and PC2
pc2range=range(pcs[,2])

eigenvalues <- eigenvectors$values[1:2]
ratio=eigenvalues[1]/eigenvalues[2] # ratio~10

## Loop through all hyperparameter samples
for (i in 1:nrow(cube)){
  
  # Get parameters
  x <- cube$x_dim[i]
  
  y <- cube$y_dim[i]
  
  x_to_y <- x / y  # The rounding of x_dim above will make the ACTUAL ratio
                   # different than what is in cube, so we recalculate it here
  
  n_neurons <- cube$n_neurons[i]
  
  ## 'radius' in the LHS samples is a value between 0-1. We interpret the value 
  ## as the quantile of all neuron-to-neuron distances in the SOM. e.g., a 
  ## 'radius' value of 1 would return the largest possible distance between all 
  ## neurons. Here, we calculate this distance -- which we set as the initial
  ## radius of the neighborhood function.
  init_radius <- quantile2radius(
    fraction=cube$radius[i],
    x=x, 
    y=y,
    shape='hexagonal'
  )
  
  dist <- cube$distance_fnc[1]
  
  neighborhood <- cube$neighborhood_fnc[i]
  
  # Sequence of neurons along pc1
  d1 <- seq(
    from=pc1range[1],
    to=pc1range[2],
    length.out=x
  )
  
  # Sequence of neurons along pc2
  d2 <- seq(
    from=pc2range[1],
    to=pc2range[2],
    length.out=y
  )
  
  # Create rectangular matrix where d1 is repeated for every value of d2 (in pc 
  # space)
  pc_grid <- expand.grid(d1, d2)
  
  # Create neuron initialization matrix by transforming pc_grid into data space
  init <- as.matrix(pc_grid) %*% t(rm)
  
  
  # Fit SOM
  
  ## Kohonen 2013 recommends the batch algorithm vs the original stepwise 
  ## recursive approach "because it is faster and safer".
  temp_SOM <- som(
    X=cumulative_timeseries_scaled,
    radius=init_radius,
    dist.fcts=dist,
    grid=somgrid(
      xdim=x,
      ydim=y,
      topo='hexagonal',
      toroidal=FALSE,
      neighbourhood.fct = neighborhood
    ),
    rlen=24,
    keep.data=TRUE,
    init=init,
    mode='pbatch'
  )
  
  # Calculate and store fit metrics
  temp_metrics <- somQuality(
    som=temp_SOM,
    traindat=cumulative_timeseries_scaled
  )
  metrics_df[i,] <- temp_metrics[1:4]
  
  if (i == 1){
    colnames(metrics_df) <- names(temp_metrics)[1:4]
  }
  
  # Store parameters
  params_df[i,] <- c(
    x, y, x_to_y, n_neurons, init_radius, dist, neighborhood, cube$radius[i]
    )
  
}

# Combined param and metrics dataframes
params_metrics_df <- cbind(params_df, metrics_df)

## Convert numeric columns to numeric data type
params_metrics_df <- params_metrics_df %>% 
  mutate(x_dim=as.numeric(x_dim)) %>% 
  mutate(y_dim=as.numeric(y_dim)) %>% 
  mutate(x_y_ratio=as.numeric(x_y_ratio)) %>% 
  mutate(n_neurons=as.numeric(n_neurons)) %>% 
  mutate(init_radius=as.numeric(init_radius)) %>% 
  mutate(radius_fraction=as.numeric(radius_fraction)) %>% 
  mutate(err.quant=as.numeric(err.quant)) %>% 
  mutate(err.topo=as.numeric(err.topo)) %>% 
  mutate(err.varratio=as.numeric(err.varratio))

write.csv(
  x=params_metrics_df,
  file='output/som_params_metrics.csv',
  row.names=FALSE
)



# COMPARE SOM CONFIGS WITH PARALLEL COORDINATES PLOT

# Associate string parameters with numeric for plotting
params_metrics_df$distance_value <- revalue(
  x=params_metrics_df$distance_fnc,
  replace=c('manhattan' = 0, 'euclidean' = 1, 'sumofsquares' = 2)
)

params_metrics_df$neighborhood_value <- revalue(
  x=params_metrics_df$neighborhood_fnc,
  replace=c('bubble' = 0, 'gaussian' = 1)
)

# Create parallel coordinates plot
fig <-
  
  params_metrics_df %>% 
  
  plot_ly(
    type='parcoords',
    
    line=list(
      color = ~err.topo, 
      hoverinfo = 'y', 
      width = 3
    ),
    
    dimensions=list(
      list(range = c(3, 60),
           label = 'X Dim', values = ~x_dim),
      
      list(range = c(2, 5),
           label = 'Y Dim', values = ~y_dim),
      
      list(range = c(0, 1),
           label = 'Radius', values = ~radius_fraction),
      
      list(tickvals=c(0, 1, 2), 
           ticktext=c('Manhattan', 'Euclidean', 'Sum of Squares'),
           label='Distance Measure', values = ~distance_value),
      
      list(tickvals=c(0, 1), 
           ticktext=c('Bubble', 'Gaussian'),
           label='Neighborhood Function', values = ~neighborhood_value),
      
      list(range = c(0, 300),
           label = 'Number of Neurons', values = ~n_neurons),
      
      list(range = c(1, 12),
           label = 'X/Y Ratio', values = ~x_y_ratio),
      
      list(range = c(1,9),
           label = 'Quant. Error', values = ~err.quant),
      
      list(range = c(0, 1),
           label = 'Topo Error', values = ~err.topo)
    )
  )

# Edit fig style and save as interactive html widget
m <- list(
  l = 100,
  r = 100,
  b = 10,
  t = 10,
  pad = 3
)

fig <- fig %>% layout(margin=m, font=list(size=16))

saveWidget(as_widget(fig), 'output/figs/som_params_metrics_pc_plot.html')



# PERFORM EPOCH TEST ON SELECTED SOM CONFIG

idx = 442  # best config chosen from parallel coordinates plot

## Pre-allocate matrix to store quality of fit metrics
epochs_to_test <- 10:30
epochs_metrics_df <- matrix(
  nrow=length(epochs_to_test),
  ncol=4
)
epochs_metrics_df <- data.frame(epochs_metrics_df)

# Create neuron initialization matrix
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

# Loop through all epochs
for (i in 1:length(epochs_to_test)){
  
  epochs <- epochs_to_test[i]
  
  temp_som <- som(
    X=cumulative_timeseries_scaled,
    radius=params_metrics_df$init_radius[idx],
    dist.fcts = params_metrics_df$distance_fnc[idx],
    grid=somgrid(
      xdim=params_metrics_df$x_dim[idx],
      ydim=params_metrics_df$y_dim[idx],
      topo='hexagonal',
      toroidal=FALSE,
      neighbourhood.fct=params_metrics_df$neighborhood_fnc[idx]
    ),
    rlen=epochs,
    keep.data=TRUE,
    init=init,
    mode='pbatch'
  )
  
  # Calculate and store fit metrics
  temp_metrics <- somQuality(
    som=temp_SOM,
    traindat=cumulative_timeseries_scaled
  )
  epochs_metrics_df[i,] <- temp_metrics[1:4]
  
  if (i == 1){
    colnames(epochs_metrics_df) <- names(temp_metrics)[1:4]
  }
  
}

epochs_metrics_df$epochs <- epochs_to_test

write.csv(
  x=epochs_metrics_df,
  file='output/som_epochs_metrics.csv',
  row.names=FALSE
)

