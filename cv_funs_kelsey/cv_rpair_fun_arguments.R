# internal deadline - package finished

glmnet(
  x,
  y,
  family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"), # loss
  # outcome.type = c("pair", "survival", "ordinal", "continuous")
  weights = NULL, # X
  offset = NULL, # X
  alpha = 1,
  nlambda = 100,
  lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL,
  standardize = TRUE,
  intercept = TRUE, # X
  thresh = 1e-07,
  dfmax = nvars + 1,
  pmax = min(dfmax * 2 + 20, nvars),
  exclude = NULL, # X
  penalty.factor = rep(1, nvars),
  lower.limits = -Inf,
  upper.limits = Inf,
  maxit = 1e+05,
  type.gaussian = ifelse(nvars < 500, "covariance", "naive"), # X
  type.logistic = c("Newton", "modified.Newton"),
  standardize.response = FALSE, # X
  type.multinomial = c("ungrouped", "grouped"), # X
  relax = FALSE, # X
  trace.it = 0, # X
  ...
)
