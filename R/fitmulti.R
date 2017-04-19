sample.uv <- function(old.v, sigma.sq.z,
                      Sigma.u.inv, Sigma.v.inv, XtX, Xty) {
  p <- length(old.v)

  # First sample new u given old v
  A.u.inv <- (XtX*(old.v%*%t(old.v)) + Sigma.u.inv)/sigma.sq.z
  A.u.inv.eig <- eigen(A.u.inv)
  A.u <- tcrossprod(tcrossprod(A.u.inv.eig$vectors[, A.u.inv.eig$values > 0],
                               diag(1/A.u.inv.eig$values[A.u.inv.eig$values > 0])),
                    A.u.inv.eig$vectors[, A.u.inv.eig$values > 0])
  A.u.rt <- tcrossprod(tcrossprod(A.u.inv.eig$vectors[, A.u.inv.eig$values > 0],
                                  diag(1/sqrt(A.u.inv.eig$values[A.u.inv.eig$values > 0]))),
                       A.u.inv.eig$vectors[, A.u.inv.eig$values > 0])
  b.u <- Xty*old.v/sigma.sq.z

  u <- crossprod(A.u, b.u) + crossprod(A.u.rt, rnorm(p))

  # Now sample new v given updated u
  A.v.inv <- (XtX*(u%*%t(u)) + Sigma.v.inv)/sigma.sq.z
  A.v.inv.eig <- eigen(A.v.inv)
  A.v <- tcrossprod(tcrossprod(A.v.inv.eig$vectors[, A.v.inv.eig$values > 0],
                               diag(1/A.v.inv.eig$values[A.v.inv.eig$values > 0])),
                    A.v.inv.eig$vectors[, A.v.inv.eig$values > 0])
  A.v.rt <- tcrossprod(tcrossprod(A.v.inv.eig$vectors[, A.v.inv.eig$values > 0],
                                  diag(1/sqrt(A.v.inv.eig$values[A.v.inv.eig$values > 0]))),
                       A.v.inv.eig$vectors[, A.v.inv.eig$values > 0])
  b.v <- Xty*u/sigma.sq.z

  v <- crossprod(A.v, b.v) + crossprod(A.v.rt, rnorm(p))

  return(cbind(u, v))

}
#' @export
mp.mcmc <- function(X, y, sigma.sq.z,
                      Sigma.u.inv = diag(ncol(X)), Sigma.v.inv = diag(ncol(X)), num.samp = 10000) {

  p <- ncol(X)
  eps <- 10^(-3)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)

  old.v <- crossprod(solve(XtX + eps*diag(p)), Xty)

  samples <- array(dim = c(num.samp, p))

  for (i in 1:num.samp) {
    s <- sample.uv(old.v, sigma.sq.z,
                   Sigma.u.inv, Sigma.v.inv, XtX, Xty)
    samples[i, ] <- s[, 1]*s[, 2]
    old.v <- s[, 2]
  }

  return(samples)

}

pv.bxy <- function(v, beta, sigma.sq.z = 1, Sigma.u = diag(length(beta)),
                   Sigma.v.inv = diag(length(beta)), log = TRUE) {
  log.lik <- -log(det(Sigma.u*(v%*%t(v))))/2 -
    sum(diag(t(v)%*%Sigma.v.inv%*%v))/(2*sigma.sq.z) -
    sum(diag(t(beta)%*%solve(Sigma.u*(v%*%t(v)))%*%beta))/(2*sigma.sq.z)
  if (log) {
    return(log.lik)
  } else {
    return(exp(log.lik))
  }
}
#' @export
mp.em <- function(Q, l, sigma.sq.z = 1,
                    Sigma.u.inv = diag(ncol(Q)), Sigma.v.inv = diag(ncol(Q)),
                    tol = 10^(-7),
                    mh = list("samps" = 10000, "burn" = 1000),
                    print.iter = FALSE, tune = 0.5, max.iter = 1000) {

  p <- length(l)
  beta.old <- rep(Inf, p)
  beta.new <- crossprod(solve(Q + 0.00001*diag(p)), l)
  j <- 1
  while (j <= max.iter & sum((beta.old - beta.new)^2) > tol) {
    v.old <- rep(1, p)
    samps <- mh[["samps"]]
    burn <- mh[["burn"]]
    vs <- matrix(nrow = samps, ncol = p)
    acc <- matrix(nrow = samps, ncol = p)
    for (i in 1:(samps + burn)) {
      accep <- rep(0, p)
      for (k in 1:p) {
        v.new <- v.old
        v.new[k] <- v.old[k] + tune*rnorm(1)
        pr.new <- pv.bxy(v.new, beta = beta, sigma.sq.z = 1, Sigma.u = solve(Sigma.u.inv),
                         Sigma.v.inv = solve(Sigma.v))
        pr.old <- pv.bxy(v.old, beta = beta, sigma.sq.z = 1, Sigma.u = solve(Sigma.u.inv),
                         Sigma.v.inv = solve(Sigma.v))
        if (pr.new > pr.old | runif(1) <= exp(pr.new - pr.old)) {
          accep[k] <- 1
          v.old <- v.new
        }
      }
      if (i > burn) {
        vs[i - burn, ] <- v.old
        acc[i - burn, ] <- accep
      }
    }
    if(print.iter) {
      cat("j = ", j, "\n")
      cat("Diff=", sum((beta.old - beta.new)^2), "\n")
      cat("Acc=", min(colMeans(acc)), "\n")
    }

    e.Sigma.u.vvt.inv <- matrix(rowMeans(apply(vs, 1, function(x) {as.vector(solve(Sigma.u*(x%*%t(x))))})), nrow = p, ncol = p)

    beta.old <- beta.new
    beta.new <- crossprod(solve(Q + e.Sigma.u.vvt.inv),l)
    j <- j + 1
  }
  return(list("beta" = beta.new, "conv" = sum((beta.old - beta.new)^2) <= tol,
              "vs" = vs))
}
