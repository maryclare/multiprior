sample.rho <- function(old, sigma.sq.z, tau.sq,
                       rho.old, pr, tune, C.inv.old) {

  p <- nrow(old)

  # Draw new value of rho
  z.old <- log(((rho.old + 1)/2)/(1 - (rho.old + 1)/2))
  z.new <- z.old + tune*rnorm(1)
  rho.new <- -1 + 2*(exp(z.new)/(1 + exp(z.new)))

  # Compute proposal probabilities
  pr.old <- dnorm(z.old, z.new, sd = tune, log = TRUE) - log(1 - rho.old^2)
  pr.new <- dnorm(z.new, z.old, sd = tune, log = TRUE) - log(1 - rho.new^2)

  # .b gives the more standard way of computing the acc. prob, same answer
  # pr.old.b <- dnorm(z.old, z.new, sd = tune, log = TRUE)
  # pr.new.b <- dnorm(z.new, z.old, sd = tune, log = TRUE)

  C.inv.new <- diag(p)
  for (i in 1:p) {
    if (i %in% c(1, p)) {
      C.inv.new[i, i] <- (1 - rho.new^2)^(p - 2)
    } else {
      C.inv.new[i, i] <- (-1)^(p - 2)*(-1 + rho.new)^(p - 2)*(rho.new + 1)^(p - 2)*(1 + rho.new^2)
    }
    if (i < p) {
      C.inv.new[i, i + 1] <- -rho.new*(1 - rho.new^2)^(p - 2)
      C.inv.new[i + 1, i] <- -rho.new*(1 - rho.new^2)^(p - 2)
    }
  }
  C.inv.new <- C.inv.new/(1 - rho.new^2)^(p - 1)

  # Compute likelihood
  ll.old <- -(p - 1)*log((1 - rho.old^2))/2 - tcrossprod(crossprod(old, C.inv.old), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.old + 1)/2, pr, pr, log = TRUE)
  ll.new <- -(p - 1)*log((1 - rho.new^2))/2 - tcrossprod(crossprod(old, C.inv.new), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.new + 1)/2, pr, pr, log = TRUE)

  # .b gives the more standard way of computing the acc. prob, same answer
  # ll.old.b <- -(p - 1)*log((1 - rho.old^2))/2 - tcrossprod(crossprod(old, C.inv.old), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.old + 1)/2, pr, pr, log = TRUE) + log(exp(z.old)/(1 + exp(z.old))^2)
  # ll.new.b <- -(p - 1)*log((1 - rho.new^2))/2 - tcrossprod(crossprod(old, C.inv.new), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.new + 1)/2, pr, pr, log = TRUE) + log(exp(z.new)/(1 + exp(z.new))^2)

  ratio <- min(1, exp(ll.new + pr.old - (ll.old + pr.new)))

  if (!runif(1) < ratio) {
    rho.new <- rho.old
    C.inv.new <- C.inv.old
  }

  return(list("rho" = rho.new,
              "C.inv" = C.inv.new,
              "acc" = (rho.new != rho.old)))

}

sample.tau.sq.inv <- function(old.u, old.v, sigma.sq.z, C.inv.u, C.inv.v,
                              pr.a, pr.b) {
  p <- nrow(old.u)
  ssr <- tcrossprod(crossprod(old.u, C.inv.u), t(old.u))/sigma.sq.z + tcrossprod(crossprod(old.v, C.inv.v), t(old.v))/sigma.sq.z
  b <- as.numeric(ssr)/2 + pr.b
  a <- 2*p/2 + pr.a
  return(rgamma(1, shape = a, rate = b))
}

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
    # I'm a little worried about the code below if 'old' is a matrix wtih more than 1 column,
    # Should check
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
  A.u <- tcrossprod(tcrossprod(A.u.inv.eig$vectors[, A.u.inv.eig$values > 0, drop = FALSE],
                               diag(1/A.u.inv.eig$values[A.u.inv.eig$values > 0], nrow = sum(A.u.inv.eig$values > 0), ncol = sum(A.u.inv.eig$values > 0))),
                    A.u.inv.eig$vectors[, A.u.inv.eig$values > 0, drop = FALSE])
  A.u.rt <- tcrossprod(tcrossprod(A.u.inv.eig$vectors[, A.u.inv.eig$values > 0, drop = FALSE],
                                  diag(1/sqrt(A.u.inv.eig$values[A.u.inv.eig$values > 0]), nrow = sum(A.u.inv.eig$values > 0), ncol = sum(A.u.inv.eig$values > 0))),
                       A.u.inv.eig$vectors[, A.u.inv.eig$values > 0, drop = FALSE])
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
  e.XtX <- eigen(XtX)

  eps <- -min(e.XtX$values) + 1
  D <- tcrossprod(tcrossprod(e.XtX$vectors, diag(1/(e.XtX$values + eps))), e.XtX$vectors)
  ridge.est <- crossprod(D, Xty)
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

#' @export
nd.mcmc <- function(X, y, sigma.sq.z,
                    Sigma.inv = NULL, num.samp = 10000, str = "uns",
                    burn.in = 500) {

  p <- ncol(X)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  e.XtX <- eigen(XtX)

  eps <- -min(e.XtX$values) + 1
  D <- tcrossprod(tcrossprod(e.XtX$vectors, diag(1/(e.XtX$values + eps))), e.XtX$vectors)
  ridge.est <- crossprod(D, Xty)
  old <- ridge.est

  samples.beta <- array(dim = c(num.samp + burn.in, p))
  samples.Sigma <- array(dim = c(num.samp + burn.in, p^2))
  samples.sigma.sq.z <- array(dim = c(num.samp + burn.in, 1))

  if (!is.null(Sigma.inv)) {
    S.i <- solve(Sigma.inv)
  }
  if (!is.null(Sigma.v.inv)) {
    S.v.i <- diag(p)
  }
  if (!is.null(sigma.sq.z)) {
    s.s.z <- sigma.sq.z
  } else {
    s.s.z <- 1
  }

  for (i in 1:(num.samp + burn.in)) {
    if (is.null(Sigma.inv)) {
      S.i <- sample.Sigma.inv(old = old, sigma.sq.z = s.s.z, str = str)
    }
    s <- sample.uv(rep(1, p), s.s.z,
                   S.u.i = S.i, S.v.i = diag(p), XtX, Xty)
    samples.beta[i, ] <- s[, 1]
    samples.Sigma[i, ] <- as.vector(solve(S.i))
    old <- s[, 1, drop = FALSE]
    if (is.null(sigma.sq.z)) {
      s.s.z <- 1/sample.sigma.z.inv(y = y, X = X, old.u = old, old.v = rep(0, p),
                                    S.u.i = S.i, S.v.i = diag(p))
    }
    samples.sigma.sq.z[i, ] <- s.s.z

  }

  return(list("beta" = samples.beta[(burn.in + 1):(burn.in + num.samp), ],
              "Sigma" = samples.Sigma[(burn.in + 1):(burn.in + num.samp), ],
              "sigma.sq.z" = samples.sigma.sq.z[(burn.in + 1):(burn.in + num.samp)]))

}

#' @export
mp.ar.mcmc <- function(X, y, num.samp = 10000, burn.in = 500,
                       sig.sq.inv.shape = 1/2,
                       sig.sq.inv.rate = 1/2,
                       tau.sq.inv.shape = 1/2,
                       tau.sq.inv.rate = 1/2,
                       rho.a = 2, tune = 1, samp.rho = TRUE, print.iter = FALSE) {

  p <- ncol(X)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  e.XtX <- eigen(XtX)

  eps <- -min(e.XtX$values) + 1
  D <- tcrossprod(tcrossprod(e.XtX$vectors, diag(1/(e.XtX$values + eps))), e.XtX$vectors)
  ridge.est <- crossprod(D, Xty)
  old.v <- sqrt(abs(ridge.est))
  old.u <- sign(ridge.est)*sqrt(abs(ridge.est))

  samples.beta <- array(dim = c(num.samp + burn.in, p))
  samples.vpar <- array(dim = c(num.samp + burn.in, 3))
  accs <- array(dim = c(num.samp + burn.in, 2))

  s.s.z <- 1
  rho.old.u <- 0
  rho.old.v <- 0
  C.inv.u <- diag(p)
  C.inv.v <- diag(p)
  acc.u <- 0
  acc.v <- 0

  for (i in 1:(num.samp + burn.in)) {

    if (print.iter) {cat("i = ", i, "\n")}

    t.sq <- 1/sample.tau.sq.inv(old = old.u, old.v = old.v, sigma.sq.z = s.s.z,
                                C.inv.u = C.inv.u, C.inv.v = C.inv.v,
                                pr.a = tau.sq.inv.shape, pr.b = tau.sq.inv.rate)

    if (samp.rho) {
      samp.rho.u <- sample.rho(old = old.u, sigma.sq.z = s.s.z, tau.sq = t.sq,
                               rho.old = rho.old.u, pr = rho.a, tune = tune,
                               C.inv.old = C.inv.u)
      rho.old.u <- samp.rho.u$rho
      C.inv.u <- samp.rho.u$C.inv
      acc.u <- samp.rho.u$acc

      samp.rho.v <- sample.rho(old = old.v, sigma.sq.z = s.s.z, tau.sq = t.sq,
                               rho.old = rho.old.v, pr = rho.a, tune = tune,
                               C.inv.old = C.inv.v)

      rho.old.v <- samp.rho.v$rho
      C.inv.v <- samp.rho.v$C.inv
      acc.v <- samp.rho.v$acc
    }

    S.u.i <- C.inv.u/(t.sq)
    S.v.i <- C.inv.v/(t.sq)

    s <- sample.uv(old.v, s.s.z,
                   S.u.i, S.v.i, XtX, Xty)
    samples.beta[i, ] <- s[, 1]*s[, 2]
    old.u <- s[, 1, drop = FALSE]
    old.v <- s[, 2, drop = FALSE]
    s.s.z <- 1/sample.sigma.z.inv(y = y, X = X, old.u = old.u, old.v = old.v,
                                  S.u.i = S.u.i, S.v.i = S.v.i,
                                  pr.a = sig.sq.inv.shape, pr.b = sig.sq.inv.rate)
    samples.vpar[i, ] <- c(s.s.z, t.sq, rho.old.v*rho.old.u)
    accs[i, ] <- c(acc.u, acc.v)

  }

  return(list("beta" = samples.beta[(burn.in + 1):(burn.in + num.samp), ],
              "var" = samples.vpar[(burn.in + 1):(burn.in + num.samp), ],
              "accs" = accs[(burn.in + 1):(burn.in + num.samp), ]))

}

#' @export
nd.ar.mcmc <- function(X, y, num.samp = 10000, burn.in = 500,
                       sig.sq.inv.shape = 1/2,
                       sig.sq.inv.rate = 1/2,
                       tau.sq.inv.shape = 1/2,
                       tau.sq.inv.rate = 1/2,
                       rho.a = 2, tune = 1, samp.rho = TRUE, print.iter = FALSE) {

  p <- ncol(X)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  e.XtX <- eigen(XtX)

  eps <- -min(e.XtX$values) + 1
  D <- tcrossprod(tcrossprod(e.XtX$vectors, diag(1/(e.XtX$values + eps))), e.XtX$vectors)
  ridge.est <- crossprod(D, Xty)
  old <- ridge.est

  samples.beta <- array(dim = c(num.samp + burn.in, p))
  samples.vpar <- array(dim = c(num.samp + burn.in, 3))
  accs <- array(dim = c(num.samp + burn.in, 1))

  s.s.z <- 1
  rho.old <- 0
  C.inv <- diag(p)
  acc <- 0

  for (i in 1:(num.samp + burn.in)) {

    if (print.iter) {cat("i = ", i, "\n")}

    t.sq <- 1/sample.tau.sq.inv(old.u = old, old.v = rep(0, p), sigma.sq.z = s.s.z,
                                C.inv.u = C.inv, C.inv.v = diag(p),
                                pr.a = tau.sq.inv.shape, pr.b = tau.sq.inv.rate)

    if (samp.rho) {
      samp.rho <- sample.rho(old = old, sigma.sq.z = s.s.z, tau.sq = t.sq,
                             rho.old = rho.old, pr = rho.a, tune = tune,
                             C.inv.old = C.inv)
      rho.old <- samp.rho$rho
      C.inv <- samp.rho$C.inv
      acc <- samp.rho$acc

    }

    S.i <- C.inv/(t.sq)

    s <- sample.uv(old.v = rep(1, p), s.s.z,
                   S.u.i = S.i, S.v.i = diag(p), XtX, Xty)
    samples.beta[i, ] <- s[, 1]
    old <- s[, 1, drop = FALSE]
    s.s.z <- 1/sample.sigma.z.inv(y = y, X = X, old.u = old, old.v = rep(0, p),
                                  S.u.i = S.i, S.v.i = diag(p),
                                  pr.a = sig.sq.inv.shape, pr.b = sig.sq.inv.rate)
    samples.vpar[i, ] <- c(s.s.z, t.sq, rho.old)
    accs[i, ] <- c(acc)

  }

  return(list("beta" = samples.beta[(burn.in + 1):(burn.in + num.samp), ],
              "var" = samples.vpar[(burn.in + 1):(burn.in + num.samp), ],
              "accs" = accs[(burn.in + 1):(burn.in + num.samp), ]))

}

