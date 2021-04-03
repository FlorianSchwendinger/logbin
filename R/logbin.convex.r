
stdm <- simple_triplet_diag_matrix
stzm <- simple_triplet_zero_matrix


convex_control <- function(solver = "ecos", tol = 1e-8, ceps = 1e-7, dry_run = FALSE, ...) {
    solver_control <- list(...)

}

neg_log_likelihood <- function(beta, y1, n, x, offset = 0) {
  eta <- drop(tcrossprod(beta, x))
  -sum(dbinom(y1, size = n, prob = exp(eta + offset), log = TRUE))
}

logbin.convex.fit_dense <- function(x, y0, y1, solver = "ecos", ceps = 1e-7, 
                              roi_control = list(), dry_run = FALSE) {
    if (!"tol" %in% names(roi_control)) {
      roi_control$tol <- 1e-8  
    }

    b0 <- y0 > 0L
    b1 <- y1 > 0L
    x0 <- x[b0,]
    x1 <- x[b1,]
    n <- nrow(x0) + nrow(x1)
    
    L1 <- cbind(rbind(x0, x1), matrix(0, n, 2 * nrow(x0)))
    log1exp <- function(xi, j, n_y_is_0) {
        k <- length(xi)        
        si <- c(rep.int(1, k),     3,                4,     6)
        sj <- c(seq_along(xi), k + j, k + n_y_is_0 + j, k + j)
        sv <- c(          -xi,    -1,               -1,     1)
        simple_triplet_matrix(si, sj, sv, nrow = 6L, ncol = (k + 2 * n_y_is_0))
    }
    L2 <- mapply(log1exp, split(x0, seq_len(nrow(x0))), seq_len(nrow(x0)), 
                 MoreArgs = list(n_y_is_0 = nrow(x0)), SIMPLIFY = FALSE)
    rhs <- c(c(0, 1, 0), c(0, 1, 1))
    rhs <- c(rep(-ceps, n), rep(rhs, nrow(x0)))
    cones <- c(K_lin(n), K_expp(2 * nrow(x0)))
    L <- do.call(rbind, c(list(L1), L2))

    obj <- c(y1 %*% x, double(nrow(x0)), y0[b0])
    o <- OP(objective = obj, 
            constraints = C_constraint(L, cones, rhs),
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


logbin.convex.fit <- function(x, y0, y1, solver = "ecos", ceps = 1e-7, 
                              roi_control = list(), dry_run = FALSE) {
    if (!"tol" %in% names(roi_control)) {
      roi_control$tol <- 1e-8  
    }

    b0 <- y0 > 0L
    b1 <- y1 > 0L
    x0 <- x[b0,]
    x1 <- x[b1,]
    n <- nrow(x0) + nrow(x1)
    
    L1 <- cbind(rbind(x0, x1), matrix(0, n, 2 * nrow(x0)))
    log1exp <- function(xi, j, n_y_is_0) {
        M <- matrix(0, nrow = 6, ncol = length(xi) + 2 * n_y_is_0)
        M[1, seq_along(xi)] <- -xi
        M[3, length(xi) + j] <- -1
        M[4, length(xi) + n_y_is_0 + j] <- -1
        M[6, length(xi) + j] <- 1
        M
    }
    L2 <- mapply(log1exp, split(x0, seq_len(nrow(x0))), seq_len(nrow(x0)), 
                 MoreArgs = list(n_y_is_0 = nrow(x0)), SIMPLIFY = FALSE)
    rhs <- c(c(0, 1, 0), c(0, 1, 1))
    rhs <- c(rep(-ceps, n), rep(rhs, nrow(x0)))
    cones <- c(K_lin(n), K_expp(2 * nrow(x0)))
    L <- do.call(rbind, c(list(L1), L2))

    obj <- c(y1 %*% x, double(nrow(x0)), y0[b0])
    o <- OP(objective = obj, 
            constraints = C_constraint(L, cones, rhs),
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

