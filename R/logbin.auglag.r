elog <- function(x) {
    ans <- rep_len(-Inf, length(x))
    b <- x > 0
    ans[b] <- log(x[b])
    ans
}


set_defaults <- function(x, val) {
    stopifnot(is.list(x), is.list(val))
    for (i in seq_along(val)) {
        key <- names(val)[i]
        if (!key %in% names(x)) {
            x[[key]] <- val[[i]]
        }
    }
    x
}


is_valid_ceps <- function(x) isTRUE(is.numeric(x)) & length(x) == 1L


auglag_outer_control <- function(control = list()) {
    stopifnot(is.list(control))
    default <- list(lam0 = 10, sig0 = 100, eps = 1e-07, itmax = 10000L, 
                    method = "BFGS", trace = FALSE, NMinit = FALSE, ilack.max = 6, 
                    i.scale = 1, e.scale = 1, kkt2.check = TRUE)
    modifyList(default, control[names(control) %in% names(default)])
}


auglag_optim_control <- function(control = list()) {
    stopifnot(is.list(control))
    default <- list(trace = 0, fnscale = 1, parscale = rep.int(1,
                    length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100L,
                    abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
                    beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
                    factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
    modifyList(default, control[names(control) %in% names(default)])
}


logbin.auglag.fit <- function(x, y0, y1, start, offset, ceps = 1e-7, 
                              control.outer = list(), control.optim = list()) {
    y <- y0 + y1
    nloglike <- function(beta) {
        eta <- drop(tcrossprod(beta, x)) + offset
        b0 <- y0 > 0
        b1 <- y1 > 0
        -drop(y1[b1] %*% eta[b1]) - drop(y0[b0] %*% elog(1 - exp(eta[b0])))
    }

    gradient <- function(beta) {
        eta <- drop(tcrossprod(beta, x)) + offset
        mu <- pmin(exp(eta), 1 - .Machine$double.neg.eps)
        -drop(crossprod(y1, x)) + drop(crossprod(x, (y0 * mu) / (1 - mu)))
    }

    constraint <- function(beta) -drop(x %*% beta) - ceps
    hin.jac <- function(beta) -x
    s <- auglag(par = start, fn = nloglike, gr = gradient, 
                hin = constraint, hin.jac = hin.jac, 
                control.outer = control.outer, control.optim = control.optim)
    s
}

logbin.auglag <- function(mt, mf, Y, offset, mono, start, control, control.method, warn) {
    x <- as.matrix(model.matrix(mt, mf))
    xnames <- dimnames(x)[[2L]]
    y <- Y
    ynames <- if (is.matrix(y)) rownames(y) else names(y)

    nvars <- ncol(x)
    nobs <- NROW(y)

    n <- weights <- rep(1, nobs)
    if (is.null(offset)) offset <- rep.int(0, nobs)

    fam <- binomial(link = log)
    eval(fam$initialize)

    y1 <- round(n * y)

    if (is.null(start)) {
        start <- c(-1, 1e-3 + double(ncol(x) - 1L))
    }

    ceps <- if (is_valid_ceps(control.method$ceps)) control.method$ceps else 1e-7
    fit <- logbin.auglag.fit(x, y0 = n - y1, y1 = y1, start = start, 
                             offset = offset, ceps = ceps,
                             control.outer = auglag_outer_control(control.method), 
                             control.optim = auglag_optim_control(control.method$control.optim))

    coefficients <- fit$par
    names(coefficients) <- xnames
    eta <- drop(x %*% coefficients) + offset
    mu <- n * fam$linkinv(eta)
    residuals <- (y - (mu / n)) / fam$mu.eta(eta)
    deviance <- sum(fam$dev.resids(y, mu / n, n))
    loglik <- -neg_log_likelihood(coefficients, y1, n, x, offset)
    aic.model <- fam$aic(y, n, mu / n, weights, dev.new) + 2 * nvars
    aic.c <- aic.model + 2 * nvars * (nvars + 1) / (nobs - nvars - 1)
    wtdmu <- sum(n * y) / sum(n)
    null.deviance <- sum(fam$dev.resids(y, wtdmu, n))
    iter <- fit$outer.iterations
    names(iter) <- NULL
  
    converged <- as.logical(fit$convergence == 0)
    boundary <- any(eta >= -control$bound.tol)

    list(coefficients = coefficients, residuals = residuals, fitted.values = mu / n,
         linear.predictors = eta, deviance = deviance,
         loglik = loglik, aic = aic.model, aic.c = aic.c, null.deviance = null.deviance,
         iter = iter, prior.weights = n, df.residual = nobs - nvars, df.null = nobs - 1, 
         y = y, x = x, converged = converged, boundary = boundary)
}
