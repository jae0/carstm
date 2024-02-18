
toy_spatial_data = function(seed = 123) {
  # generate some spatial data
  set.seed(seed)
  east = north = 1:10
  Grid = expand.grid(east, north)
  K = nrow(Grid)

  # set up distance and neighbourhood matrices
  distance = as.matrix(dist(Grid))
  W = array(0, c(K, K))
  W[distance == 1] = 1   

  # generate the covariates and response data
  x1 = rnorm(K)
  x2 = rnorm(K)
  theta = rnorm(K, sd = 0.05)
  phi = mvtnorm::rmvnorm(1, mean = rep(0, K), sigma = 0.4 * exp(-0.1 * distance) )
  eta = x1 + x2 + phi
  prob = exp(eta) / (1 + exp(eta))
  size = rep(50, K)
  y = rbinom(n = K, size = size, prob = prob)
  dat = data.frame(y, size, x1, x2)
  dat$ID = 1:K
  row.names(W) = dat$ID
  
  return(list(dat=dat, W=W))
}
