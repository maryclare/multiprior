sample.rho.normal <- function(old, tau.sq, rho.other,
                              rho.old, pr, tune) {

  p <- length(old)

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
  C.inv.both.new <- diag(p)
  C.inv.both.old <- diagp(p)
  for (i in 1:p) {
    if (i %in% c(1, p)) {
      C.inv.new[i, i] <- (1 - rho.new^2)^(p - 2)
      C.inv.both.new[i, i] <- (1 - (rho.other*rho.new)^2)^(p - 2)
      C.inv.both.old[i, i] <- (1 - (rho.other*rho.old)^2)^(p - 2)
    } else {
      C.inv.new[i, i] <- (-1)^(p - 2)*(-1 + rho.new)^(p - 2)*(rho.new + 1)^(p - 2)*(1 + rho.new^2)
      C.inv.both.new[i, i] <- (-1)^(p - 2)*(-1 + rho.other*rho.new)^(p - 2)*(rho.new*rho.other + 1)^(p - 2)*(1 + (rho.other*rho.new)^2)
      C.inv.both.old[i, i] <- (-1)^(p - 2)*(-1 + rho.other*rho.old)^(p - 2)*(rho.old*rho.other + 1)^(p - 2)*(1 + (rho.other*rho.old)^2)

    }
    if (i < p) {
      C.inv.new[i, i + 1] <- -rho.new*(1 - rho.new^2)^(p - 2)
      C.inv.new[i + 1, i] <- -rho.new*(1 - rho.new^2)^(p - 2)
      C.inv.both.new[i, i + 1] <- -rho.new*rho.other*(1 - (rho.new*rho.other)^2)^(p - 2)
      C.inv.both.new[i + 1, i] <- -rho.new*rho.other*(1 - (rho.new*rho.other)^2)^(p - 2)
      C.inv.both.old[i, i + 1] <- -rho.old*rho.other*(1 - (rho.old*rho.other)^2)^(p - 2)
      C.inv.both.old[i + 1, i] <- -rho.old*rho.other*(1 - (rho.old*rho.other)^2)^(p - 2)
    }
  }
  C.inv.new <- C.inv.new/(1 - rho.new^2)^(p - 1)
  C.inv.both.new <- C.inv.both.new/(1 - (rho.other*rho.new)^2)^(p - 1)
  C.inv.both.old <- C.inv.both.old/(1 - (rho.other*rho.old)^2)^(p - 1)

  # Compute likelihood
  ll.old <- -(p - 1)*log((1 - (rho.other*rho.old)^2))/2 - tcrossprod(crossprod(old, C.inv.both.old), t(old))/(2*tau.sq) + dbeta((rho.old + 1)/2, pr, pr, log = TRUE)
  ll.new <- -(p - 1)*log((1 - (rho.other*rho.new)^2))/2 - tcrossprod(crossprod(old, C.inv.both.new), t(old))/(2*tau.sq) + dbeta((rho.new + 1)/2, pr, pr, log = TRUE)

  # .b gives the more standard way of computing the acc. prob, same answer
  # ll.old.b <- -(p - 1)*log((1 - rho.old^2))/2 - tcrossprod(crossprod(old, C.inv.old), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.old + 1)/2, pr, pr, log = TRUE) + log(exp(z.old)/(1 + exp(z.old))^2)
  # ll.new.b <- -(p - 1)*log((1 - rho.new^2))/2 - tcrossprod(crossprod(old, C.inv.new), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.new + 1)/2, pr, pr, log = TRUE) + log(exp(z.new)/(1 + exp(z.new))^2)

  ratio <- min(1, exp(ll.new + pr.old - (ll.old + pr.new)))

  if (!runif(1) < ratio) {
    rho.new <- rho.old
    C.inv.new <- C.inv.old
    C.inv.both.new <- C.inv.both.old
  }

  return(list("rho" = rho.new,
              "C.inv" = C.inv.new,
              "acc" = (rho.new != rho.old),
              "C.inv.both" = C.inv.both.new))

}

sample.rho <- function(old, tau.sq,
                       rho.old, pr, tune, C.inv.old) {

  p <- length(old)

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
  ll.old <- -(p - 1)*log((1 - rho.old^2))/2 - tcrossprod(crossprod(old, C.inv.old), t(old))/(2*tau.sq) + dbeta((rho.old + 1)/2, pr, pr, log = TRUE)
  ll.new <- -(p - 1)*log((1 - rho.new^2))/2 - tcrossprod(crossprod(old, C.inv.new), t(old))/(2*tau.sq) + dbeta((rho.new + 1)/2, pr, pr, log = TRUE)

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

sample.tau.sq.inv <- function(old, C.inv,
                              pr.a, pr.b) {
  p <- length(old)
  ssr <- tcrossprod(crossprod(old, C.inv), t(old))
  b <- as.numeric(ssr)/2 + pr.b
  a <- p/2 + pr.a
  return(rgamma(1, shape = a, rate = b))
}

sample.sigma.z.inv <- function(beta, y, X, pr.a = 1/2 + 1, pr.b = 1/2) {
  n <- length(y)
  ssr <- crossprod(y - crossprod(t(X), beta))
  b <- as.numeric(ssr)/2 + pr.b
  a <- n/2 + pr.a
  return(rgamma(1, shape = a, rate = b))
}

sample.Sigma.inv <- function(old, pr.V.inv = diag(nrow(old)),
                             pr.df = nrow(old) + 2, str) {
  p <- nrow(old)
  if (str == "uns") {
    V.inv <- tcrossprod(old) + pr.V.inv
    df <- ncol(old) + p + 2
    return(rWishart(1, df, solve(V.inv))[, , 1])
  } else if (str == "het") {
    b <- apply(old, 1, function(x) {sum(x^2)})/(2) + diag(pr.V.inv)/2
    a <- rep(ncol(old), nrow(old))/2 + pr.df/2
    return(diag(rgamma(p, shape = a, rate = b)))
  } else if (str == "con") {
    b <- sum(apply(old, 1, function(x) {sum(x^2)}))/(2) + sum(diag(pr.V.inv))/2
    # I'm a little worried about the code below if 'old' is a matrix wtih more than 1 column,
    # Should check
    a <- sum(rep(ncol(old), nrow(old)))/2 + p*pr.df/2
    return(rgamma(1, shape = a, rate = b)*diag(p))
  }
}

sample.beta <- function(sigma.sq.z, Sigma.inv, XtX, Xty) {

  A.inv <- (XtX/sigma.sq.z + Sigma.inv)

  A.inv.eig <- eigen(A.inv)
  A <- tcrossprod(tcrossprod(A.inv.eig$vectors[, A.inv.eig$values > 0],
                               diag(1/A.inv.eig$values[A.inv.eig$values > 0])),
                    A.inv.eig$vectors[, A.inv.eig$values > 0])
  A.rt <- tcrossprod(tcrossprod(A.inv.eig$vectors[, A.inv.eig$values > 0],
                                  diag(1/sqrt(A.inv.eig$values[A.inv.eig$values > 0]))),
                       A.inv.eig$vectors[, A.inv.eig$values > 0])
  b <- Xty/sigma.sq.z

  v <- crossprod(A, b) + crossprod(A.rt, rnorm(p))

  return(v)

}

sample.uv <- function(old.v, sigma.sq.z,
                      Sigma.u.inv, Sigma.v.inv, XtX, Xty) {

  u <- sample.beta(sigma.sq.z = sigma.sq.z, Sigma.inv = Sigma.u.inv, XtX = XtX*tcrossprod(old.v),
                   Xty = Xty*old.v)
  v <- sample.beta(sigma.sq.z = sigma.sq.z, Sigma.inv = Sigma.v.inv, XtX = XtX*tcrossprod(u),
                   Xty = Xty*u)

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
      S.u.i <- sample.Sigma.inv(old = old.u, str = str)
    }
    if (is.null(Sigma.v.inv)) {
      S.v.i <- sample.Sigma.inv(old = old.v, str = str)
    }
    s <- sample.uv(old.v, s.s.z,
                   S.u.i, S.v.i, XtX, Xty)
    samples.beta[i, ] <- s[, 1]*s[, 2]
    samples.Sigma[i, ] <- as.vector(solve(S.u.i)*solve(S.v.i))
    old.u <- s[, 1, drop = FALSE]
    old.v <- s[, 2, drop = FALSE]
    if (is.null(sigma.sq.z)) {
      s.s.z <- 1/sample.sigma.z.inv(beta = samples.beta[i, ], y = y, X = X)
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
  if (!is.null(sigma.sq.z)) {
    s.s.z <- sigma.sq.z
  } else {
    s.s.z <- 1
  }

  for (i in 1:(num.samp + burn.in)) {
    if (is.null(Sigma.inv)) {
      S.i <- sample.Sigma.inv(old = old, str = str)
    }
    s <- sample.uv(old = rep(1, p), s.s.z,
                   Sigma.u.inv = S.i, Sigma.v.inv = diag(p), XtX, Xty)
    samples.beta[i, ] <- s[, 1]
    samples.Sigma[i, ] <- as.vector(solve(S.i))
    old <- s[, 1, drop = FALSE]
    if (is.null(sigma.sq.z)) {
      s.s.z <- 1/sample.sigma.z.inv(beta = samples.beta[i, ], y = y, X = X)
    }
    samples.sigma.sq.z[i, ] <- s.s.z

  }

  return(list("beta" = samples.beta[(burn.in + 1):(burn.in + num.samp), ],
              "Sigma" = samples.Sigma[(burn.in + 1):(burn.in + num.samp), ],
              "sigma.sq.z" = samples.sigma.sq.z[(burn.in + 1):(burn.in + num.samp)]))

}

#' @export
mp.ar.mcmc <- function(X, y, num.samp = 10000, burn.in = 500,
                       sig.sq.inv.shape = 1/2 + 1,
                       sig.sq.inv.rate = 1/2 + 1,
                       tau.sq.inv.shape = 1/2 + 1,
                       tau.sq.inv.rate = 1/2 + 1,
                       rho.a = 1, tune = 1, samp.rho = TRUE, print.iter = FALSE) {

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

    t.sq.u <- 1/sample.tau.sq.inv(old = old.u,
                                C.inv = C.inv.u,
                                pr.a = tau.sq.inv.shape, pr.b = tau.sq.inv.rate)
    t.sq.v <- 1/sample.tau.sq.inv(old = old.v,
                                  C.inv = C.inv.v,
                                  pr.a = tau.sq.inv.shape, pr.b = tau.sq.inv.rate)

    if (samp.rho) {
      samp.rho.u <- sample.rho(old = old.u, tau.sq = t.sq.u,
                               rho.old = rho.old.u, pr = rho.a, tune = tune,
                               C.inv.old = C.inv.u)
      rho.old.u <- samp.rho.u$rho
      C.inv.u <- samp.rho.u$C.inv
      acc.u <- samp.rho.u$acc

      samp.rho.v <- sample.rho(old = old.v, tau.sq = t.sq.v,
                               rho.old = rho.old.v, pr = rho.a, tune = tune,
                               C.inv.old = C.inv.v)

      rho.old.v <- samp.rho.v$rho
      C.inv.v <- samp.rho.v$C.inv
      acc.v <- samp.rho.v$acc
    }

    S.u.i <- C.inv.u/(t.sq.u)
    S.v.i <- C.inv.v/(t.sq.v)

    s <- sample.uv(old.v, s.s.z,
                   S.u.i, S.v.i, XtX, Xty)
    samples.beta[i, ] <- s[, 1]*s[, 2]
    old.u <- s[, 1, drop = FALSE]
    old.v <- s[, 2, drop = FALSE]
    s.s.z <- 1/sample.sigma.z.inv(beta = samples.beta[i, ], y = y, X = X,
                                  pr.a = sig.sq.inv.shape, pr.b = sig.sq.inv.rate)
    samples.vpar[i, ] <- c(s.s.z, t.sq.u*t.sq.v, rho.old.v*rho.old.u)
    accs[i, ] <- c(acc.u, acc.v)

  }

  return(list("beta" = samples.beta[(burn.in + 1):(burn.in + num.samp), ],
              "var" = samples.vpar[(burn.in + 1):(burn.in + num.samp), ],
              "accs" = accs[(burn.in + 1):(burn.in + num.samp), ]))

}

#' @export
nd.ar.mcmc <- function(X, y, num.samp = 10000, burn.in = 500,
                       sig.sq.inv.shape = 1/2 + 1,
                       sig.sq.inv.rate = 1/2 + 1,
                       tau.sq.inv.shape = 1/2 + 1,
                       tau.sq.inv.rate = 1/2 + 1,
                       rho.a = 1, tune = 1, samp.rho = TRUE, print.iter = FALSE) {

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
  accs <- array(dim = c(num.samp + burn.in, 2))

  s.s.z <- t.sq.u <- t.sq.v <- 1
  rho.old.u <- 0
  rho.old.v <- 0
  C.inv.u <- diag(p)
  C.inv.v <- diag(p)
  acc.u <- 0
  acc.v <- 0

  for (i in 1:(num.samp + burn.in)) {

    if (print.iter) {cat("i = ", i, "\n")}

    t.sq.u <- 1/sample.tau.sq.inv(old = old,
                                  C.inv = t.sq.v*C.inv.u*C.inv.v,
                                  pr.a = tau.sq.inv.shape, pr.b = tau.sq.inv.rate)
    t.sq.v <- 1/sample.tau.sq.inv(old = old,
                                  C.inv = t.sq.u*C.inv.u*C.inv.v,
                                  pr.a = tau.sq.inv.shape, pr.b = tau.sq.inv.rate)

    if (samp.rho) {
      samp.rho.u <- sample.rho.normal(old = old, tau.sq = t.sq.u*t.sq.v, rho.other = rho.old.v,
                                      rho.old = rho.old.u, pr = rho.a, tune = tune,
                                      C.inv.old = C.inv.u)
      rho.old.u <- samp.rho.u$rho
      C.inv.u <- samp.rho.u$C.inv
      acc.u <- samp.rho.u$acc

      samp.rho.v <- sample.rho.normal(old = old.v, tau.sq = t.sq.v*t.sq.u, rho.other = rho.old.u,
                                      rho.old = rho.old.v, pr = rho.a, tune = tune,
                                      C.inv.old = C.inv.v)

      rho.old.v <- samp.rho.v$rho
      C.inv.v <- samp.rho.v$C.inv
      acc.v <- samp.rho.v$acc
      C.inv.both <- samp.rho.v$C.inv.both
    }

    S.i <- C.inv.both/(t.sq.u*t.sq.v)



    samples.beta[i, ] <- sample.beta(sigma.sq.z = sigma.sq.z, Sigma.inv = S.i,
                                     XtX = XtX, Xty = Xty)
    old <- matrix(samples.beta[i, ], nrow = p, ncol = 1)
    s.s.z <- 1/sample.sigma.z.inv(beta = samples.beta[i, ], y = y, X = X,
                                  pr.a = sig.sq.inv.shape, pr.b = sig.sq.inv.rate)
    samples.vpar[i, ] <- c(s.s.z, t.sq.u*t.sq.v, rho.old.v*rho.old.u)
    accs[i, ] <- c(acc.u, acc.v)

  }

  return(list("beta" = samples.beta[(burn.in + 1):(burn.in + num.samp), ],
              "var" = samples.vpar[(burn.in + 1):(burn.in + num.samp), ],
              "accs" = accs[(burn.in + 1):(burn.in + num.samp), ]))

}

