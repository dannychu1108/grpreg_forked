# Custom newXG function that skips standardization
custom_newXG <- function(X, group, group.multiplier, ...) {
  # Ensure X is a matrix
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # Custom processing to avoid standardization
  # (This is a placeholder; replace with the actual logic from newXG)
  
  # Assuming group is a factor or vector of consecutive integers
  g <- as.integer(group)
  
  # Return the processed design matrix and group info
  list(X = X, g = g, m = group.multiplier, scale = rep(1, ncol(X)), reorder = FALSE, names = colnames(X))
}

# Custom grpsurv function using the custom newXG
custom_grpsurv <- function(X, y, group=1:ncol(X), penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),
                           gamma=ifelse(penalty=="grSCAD", 4, 3), alpha=1, nlambda=100, lambda,
                           lambda.min={if (nrow(X) > ncol(X)) 0.001 else .05}, eps=.001, max.iter=10000,
                           dfmax=ncol(X), gmax=length(unique(group)), tau=1/3,
                           group.multiplier=rep(1, length(unique(group))), warn=TRUE, returnX=FALSE, ...) {
  
  penalty <- match.arg(penalty)
  if (gamma <= 1 & penalty %in% c("grMCP", "cMCP")) stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="grSCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0, 1]", call.=FALSE)
  
  Y <- newS(y)
  XG <- custom_newXG(X[Y$ind, , drop=FALSE], group, group.multiplier, 1)
  if (nrow(XG$X) != length(Y$fail)) stop("X and y do not have the same number of observations", call.=FALSE)
  
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XG$X, Y$time, Y$fail, XG$g, penalty, alpha, lambda.min, nlambda, XG$m)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  n <- length(Y$time)
  p <- ncol(XG$X)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  if (strtrim(penalty, 2) != "gr") {
    res <- .Call("lcdfit_cox", XG$X, Y$fail, penalty, K1, K0, lambda, alpha, eps, 0, gamma, tau, as.integer(max.iter),
                 XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
  } else {
    res <- .Call("gdfit_cox", XG$X, Y$fail, penalty, K1, K0, lambda, alpha, eps, as.integer(max.iter),
                 as.double(gamma), XG$m, as.integer(dfmax), as.integer(gmax), as.integer(warn), as.integer(user.lambda))
  }
  b <- matrix(res[[1]], p, nlambda)
  iter <- res[[2]]
  df <- res[[3]]
  loss <- -2*res[[4]]
  Eta <- matrix(res[[5]], n, nlambda)
  
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  loss <- loss[ind]
  Eta <- Eta[, ind, drop=FALSE]
  if (iter[1] == max.iter) stop("Algorithm failed to converge for any values of lambda.  This indicates a combination of (a) an ill-conditioned feature matrix X and (b) insufficient penalization.  You must fix one or the other for your model to be identifiable.", call.=FALSE)
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda", call.=FALSE)
  
  beta <- b
  
  dimnames(beta) <- list(XG$names, round(lambda, digits=4))
  colnames(Eta) <- round(lambda, digits=4)
  
  val <- structure(list(beta = beta,
                        group = factor(group),
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loss = loss,
                        n = n,
                        df = df,
                        iter = iter,
                        group.multiplier = XG$m,
                        time = Y$time,
                        fail = Y$fail,
                        order = Y$ind,
                        linear.predictors = sweep(Eta, 2, colMeans(Eta), '-')),
                   class = c("grpsurv", "grpreg"))
  if (returnX) val$XG <- XG
  val
}
