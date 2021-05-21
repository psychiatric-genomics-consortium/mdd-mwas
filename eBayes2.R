
# Updated ebayes function by Rosie to get round integer overflow problem (essentially adds as.numeric() below)
ebayes2 <- function (fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) 
{
    coefficients <- fit$coefficients
    stdev.unscaled <- fit$stdev.unscaled
    sigma <- fit$sigma
    df.residual <- fit$df.residual
    if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || 
        is.null(df.residual)) 
        stop("No data, or argument is not a valid lmFit object")
    if (all(df.residual == 0)) 
        stop("No residual degrees of freedom in linear model fits")
    if (all(!is.finite(sigma))) 
        stop("No finite residual standard deviations")
    if (trend) {
        covariate <- fit$Amean
        if (is.null(covariate)) 
            stop("Need Amean component in fit to estimate trend")
    }
    else {
        covariate <- NULL
    }
    out <- squeezeVar(sigma^2, df.residual, covariate = covariate, 
        robust = robust, winsor.tail.p = winsor.tail.p)
    out$s2.prior <- out$var.prior
    out$s2.post <- out$var.post
    out$var.prior <- out$var.post <- NULL
    out$t <- coefficients/stdev.unscaled/sqrt(out$s2.post)
    df.total <- df.residual + out$df.prior
    df.pooled <- sum(as.numeric(df.residual), na.rm = TRUE)
    df.total <- pmin(df.total, df.pooled)
    out$df.total <- df.total
out$p.value <- 2 * pt(-abs(out$t), df = df.total)
    var.prior.lim <- stdev.coef.lim^2/median(out$s2.prior)
    out$var.prior <- tmixture.matrix(out$t, stdev.unscaled, df.total, 
        proportion, var.prior.lim)
    if (any(is.na(out$var.prior))) {
        out$var.prior[is.na(out$var.prior)] <- 1/out$s2.prior
        warning("Estimation of var.prior failed - set to default value")
    }
    r <- rep(1, NROW(out$t)) %o% out$var.prior
    r <- (stdev.unscaled^2 + r)/stdev.unscaled^2
    t2 <- out$t^2
    Infdf <- out$df.prior > 10^6
    if (any(Infdf)) {
        kernel <- t2 * (1 - 1/r)/2
        if (any(!Infdf)) {
            t2.f <- t2[!Infdf]
            r.f <- r[!Infdf]
            df.total.f <- df.total[!Infdf]
            kernel[!Infdf] <- (1 + df.total.f)/2 * log((t2.f + df.total.f)/(t2.f/r.f + df.total.f))
        }
    }
    else kernel <- (1 + df.total)/2 * log((t2 + df.total)/(t2/r + df.total))
	out$lods <- log(proportion/(1 - proportion)) - log(r)/2 + kernel
    out
}

eBayes2 <- function (fit, proportion=0.01, stdev.coef.lim=c(0.1, 4), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05, 0.1)) 
{
    if (trend) 
        if (is.null(fit$Amean)) 
            stop("Need Amean component in fit to estimate trend")
    eb <- ebayes2(fit=fit, proportion=proportion, stdev.coef.lim=stdev.coef.lim, trend=trend, robust=robust, winsor.tail.p=winsor.tail.p)
    fit$df.prior <- eb$df.prior
    fit$s2.prior <- eb$s2.prior
    fit$var.prior <- eb$var.prior
    fit$proportion <- proportion
    fit$s2.post <- eb$s2.post
    fit$t <- eb$t
    fit$df.total <- eb$df.total
    fit$p.value <- eb$p.value
    fit$lods <- eb$lods
    if (!is.null(fit$design) && is.fullrank(fit$design)) {
        F.stat <- classifyTestsF(fit, fstat.only = TRUE)
        fit$F <- as.vector(F.stat)
        df1 <- attr(F.stat, "df1")
        df2 <- attr(F.stat, "df2")
        if (df2[1] > 1e+06) 
            fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)
        else fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
    }
    fit
}
