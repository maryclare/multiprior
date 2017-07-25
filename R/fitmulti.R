sample.sigma.z.inv <- function(y, X, old.u, old.v, S.u.i, S.v.i, pr.a = 1/2, pr.b = 1/2) {
  p <- nrow(old.u)
  n <- length(y)
  beta <- old.u*old.v
  ssr <- crossprod(y - crossprod(t(X), beta))
  r.u <- crossprod(old.u, crossprod(S.u.i, old.u))
  r.v <- crossprod(old.v, crossprod(S.v.i, old.v))
  b <- as.numeric(ssr + r.u + r.v)/2 + pr.b
  a <- n/2 + p/2 + p/2 + pr.a
  return(rgamma(1, shape = a, rate = b))
}

sample.Sigma.inv <- function(old, sigma.sq.z, pr.V.inv = diag(nrow(old)),
                             pr.df = nrow(old) + 2, str) {
  p <- nrow(old)
  if (str == "uns") {
    V.inv <- tcrossprod(old)/sigma.sq.z + pr.V.inv
    df <- ncol(old) + p + 2
    return(rWishart(1, df, solve(V.inv))[, , 1])
  } else if (str == "het") {
    b <- apply(old, 1, function(x) {sum(x^2)})/(2*sigma.sq.z) + diag(pr.V.inv)/2
    a <- rep(ncol(old), nrow(old))/2 + pr.df/2
    return(diag(rgamma(p, shape = a, rate = b)))
  } else if (str == "con") {
    b <- sum(apply(old, 1, function(x) {sum(x^2)}))/(2*sigma.sq.z) + sum(diag(pr.V.inv))/2
    a <- sum(rep(ncol(old), nrow(old)))/2 + p*pr.df/2
    return(rgamma(1, shape = a, rate = b)*diag(p))
  }
}

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
                      Sigma.u.inv = NULL, Sigma.v.inv = NULL, num.samp = 10000,
                    str = "uns", burn.in = 500) {

  p <- ncol(X)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)

  eps <- -min(eigen(XtX)$values) + 1
  ridge.est <- crossprod(solve(XtX + eps*diag(p)), Xty)
  old.v <- sqrt(abs(ridge.est))
  old.u <- sign(ridge.est)*sqrt(abs(ridge.est))

  samples.beta <- array(dim = c(num.samp + burn.in, p))
  samples.Sigma <- array(dim = c(num.samp + burn.in, p^2))
  samples.sigma.sq.z <- array(dim = c(num.samp + burn.in, 1))

  if (!is.null(Sigma.u.inv)) {
    S.u.i <- Sigma.u.inv
  }
  if (!is.null(Sigma.v.inv)) {
    S.v.i <- Sigma.v.inv
  }
  if (!is.null(sigma.sq.z)) {
    s.s.z <- sigma.sq.z
  } else {
    s.s.z <- 1
  }

  for (i in 1:(num.samp + burn.in)) {
    if (is.null(Sigma.u.inv)) {
      S.u.i <- sample.Sigma.inv(old = old.u, sigma.sq.z = s.s.z, str = str)
    }
    if (is.null(Sigma.v.inv)) {
      S.v.i <- sample.Sigma.inv(old = old.v, sigma.sq.z = s.s.z, str = str)
    }
    s <- sample.uv(old.v, s.s.z,
                   S.u.i, S.v.i, XtX, Xty)
    samples.beta[i, ] <- s[, 1]*s[, 2]
    samples.Sigma[i, ] <- as.vector(solve(S.u.i)*solve(S.v.i))
    old.u <- s[, 1, drop = FALSE]
    old.v <- s[, 2, drop = FALSE]
    if (is.null(sigma.sq.z)) {
      s.s.z <- 1/sample.sigma.z.inv(y = y, X = X, old.u = old.u, old.v = old.v,
                                    S.u.i = S.u.i, S.v.i = S.v.i)
    }
    samples.sigma.sq.z[i, ] <- s.s.z

  }

  return(list("beta" = samples.beta[(burn.in + 1):(burn.in + num.samp), ],
              "Sigma" = samples.Sigma[(burn.in + 1):(burn.in + num.samp), ],
              "sigma.sq.z" = samples.sigma.sq.z[(burn.in + 1):(burn.in + num.samp)]))

}
