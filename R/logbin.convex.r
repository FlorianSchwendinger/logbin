
stdm <- simple_triplet_diag_matrix
stzm <- simple_triplet_zero_matrix


convex_control <- function(solver = "ecos", tol = 1e-8, ceps = 1e-7, dry_run = FALSE, ...) {
    solver_control <- list(...)

}

neg_log_likelihood <- function(beta, y1, n, x, offset = 0) {
  eta <- drop(tcrossprod(beta, x))
  -sum(dbinom(y1, size = n, prob = exp(eta + offset), log = TRUE))
}


logbin.convex.fit <- function(x, y0, y1, solver = "ecos", ceps = 1e-7,
                              roi_control = list(), dry_run = FALSE) {
    b0 <- y0 > 0L
    x0 <- as.simple_triplet_matrix(x[ b0, , drop = FALSE])
    x1 <- as.simple_triplet_matrix(x[!b0, , drop = FALSE])
    n <- nrow(x0) + nrow(x1)

    L1 <- cbind(rbind(x0, x1), stzm(n, 2 * nrow(x0)))
    K1 <- K_lin(n)
    rhs1 <- rep(-ceps, nrow(x))

    k <- ncol(x)
    sv <- as.vector(rbind(t(-x[b0, ]), -1, -1, 1))
    si <- c(rep.int(1, k), 3L, 4L, 6L)
    si <- si + rep(6L * (seq_len(nrow(x0)) - 1L), each = length(si))
    sj <- c(k, k + nrow(x0), k)
    sj <- c(rbind(matrix(seq_len(k), nrow = k, ncol = nrow(x0)),
            matrix(sj + rep(seq_len(nrow(x0)), each = 3), nrow = 3)))
    L2 <- simple_triplet_matrix(si, sj, sv, nrow = 6L * nrow(x0), 
                                ncol = (k + 2 * nrow(x0)))
    K2 <- K_expp(2 * nrow(x0))
    rhs2 <- rep(c(c(0, 1, 0), c(0, 1, 1)), nrow(x0))

    obj <- c(y1 %*% x, double(nrow(x0)), y0[b0])
    o <- OP(objective = obj, 
            constraints = C_constraint(rbind(L1, L2), c(K1, K2), c(rhs1, rhs2)),
            bounds = V_bound(ld = -Inf, nobj = length(obj)),
            maximum = TRUE)

    if ( dry_run ) return(o)
    so <- ROI_solve(o, solver = solver, control = roi_control)
    rval <- head(solution(so, force = TRUE), NCOL(x))
    if ( solution(so, "status_code") != 0L ) {
        warning("Solution: ", paste(solution(so, "status")$msg$message, collapse = " "))
    }
    attributes(rval)$solution <- so
    return(rval)
}



# logbin.ab  <- function(mt, mf, Y, offset, mono, start, control, control.method, warn)
# logbin.cem <- function(mt, mf, Y, offset, mono, start, control, accelerate, control.method, warn)

# param mono currently not supported by function logbin.convex
# param start not needed for function logbin.convex
# param offset TODO
logbin.convex <- function(mt, mf, Y, offset = NULL, mono = NULL, start = NULL, control, control.method, warn) {
    x <- as.matrix(model.matrix(mt, mf))
    xnames <- dimnames(x)[[2L]]
    y <- Y
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    
    nvars <- ncol(x)
    nobs <- NROW(y)
    
    n <- weights <- rep(1, nobs)
    if (!is.null(offset)) {
        stop("method convex does not support offset")
    }
    
    fam <- binomial(link = log)
    eval(fam$initialize)
    
    y1 <- round(n * y)

    if (is.null(control.method$method)) {
      solver <- "ecos"
    } else {
      solver <- control.method$method
      control.method$method <- NULL
    }

    ceps <- if (is_valid_ceps(control$ceps)) control$ceps else 1e-7
    control.method <- set_defaults(control.method, list(tol = 1e-8))
    fit <- logbin.convex.fit(x, y0 = n - y1, y1 = y1, solver = solver, ceps = ceps,
                             roi_control = control.method)

    coefficients <- as.numeric(fit)
    names(coefficients) <- xnames
    eta <- drop(x %*% coefficients)
    mu <- n * fam$linkinv(eta)
    residuals <- (y - (mu / n)) / fam$mu.eta(eta)
    deviance <- sum(fam$dev.resids(y, mu / n, n))
    loglik <- -neg_log_likelihood(coefficients, y1, n, x)
    aic.model <- fam$aic(y, n, mu / n, weights, dev.new) + 2 * nvars
    aic.c <- aic.model + 2 * nvars * (nvars + 1) / (nobs - nvars - 1)
    wtdmu <- sum(n * y) / sum(n)
    null.deviance <- sum(fam$dev.resids(y, wtdmu, n))
    iter <- tryCatch(unname(attr(fit, "solution")$message$retcodes["iter"]), 
                     error = function(e) NA_integer_)

    converged <- isTRUE(solution(attr(fit, "solution"), "status_code") == 0L)
    boundary <- any(eta >= -control$bound.tol)

    list(coefficients = coefficients, residuals = residuals, fitted.values = mu / n,
         linear.predictors = eta, deviance = deviance,
         loglik = loglik, aic = aic.model, aic.c = aic.c, null.deviance = null.deviance,
         iter = iter, prior.weights = n, df.residual = nobs - nvars, df.null = nobs - 1, 
         y = y, x = x, converged = converged, boundary = boundary)
}

