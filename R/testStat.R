#' Significance test of an additive term
#' @description Significance test of an additive term relying on the methodology
#' in Wood (Biometrika 2013). It is extracted from a hidden function
#' in the 'mgcv' package. The additive term is estimated using the product of a matrix <X> and a vector <p>.
#' @param p Vector of spline parameters used to estimate f(x)
#' @param X (Centered) B-spline basis evaluated on a fine regular grid on the support of variable <x>
#' @param V Posterior variance-covariance matrix of parameter <p>
#' @param rank Effective dimension of <p>
#' @param type 0 value by default
#' @param res.df -1 indicates that the scale is fixed (cf. ordinal response)
#' @return Returns a list with following elements:
#' \itemize{
#' \item{stat : \verb{ }}{Value of the test statistics}
#' \item{pval : \verb{ }}{P-value of the test for the null hypothesis Ho: p=0}
#' \item{rank : \verb{ }}{Effective dimension of <p>}
#' }
#' @keywords internal
#' @export
testStat = function (p, X, V, rank = NULL, type = 0, res.df = -1){
    qrx <- qr(X, tol = 0)
    R <- qr.R(qrx)
    V <- R %*% V[qrx$pivot, qrx$pivot, drop = FALSE] %*% t(R)
    V <- (V + t(V))/2
    ed <- eigen(V, symmetric = TRUE)
    siv <- sign(ed$vectors[1, ])
    siv[siv == 0] <- 1
    ed$vectors <- sweep(ed$vectors, 2, siv, "*")
    k <- max(0, floor(rank))
    nu <- abs(rank - k)
    if (type == 1) {
        if (rank > k + 0.05 || k == 0)
            k <- k + 1
        nu <- 0
        rank <- k
    }
    if (nu > 0)
        k1 <- k + 1
    else k1 <- k
    r.est <- sum(ed$values > max(ed$values) * .Machine$double.eps^0.9)
    if (r.est < k1) {
        k1 <- k <- r.est
        nu <- 0
        rank <- r.est
    }
    vec <- ed$vectors
    if (k1 < ncol(vec))
        vec <- vec[, 1:k1, drop = FALSE]
    if (nu > 0 && k > 0) {
        if (k > 1)
            vec[, 1:(k - 1)] <- t(t(vec[, 1:(k - 1)])/sqrt(ed$val[1:(k -
                1)]))
        b12 <- 0.5 * nu * (1 - nu)
        if (b12 < 0)
            b12 <- 0
        b12 <- sqrt(b12)
        B <- matrix(c(1, b12, b12, nu), 2, 2)
        ev <- diag(ed$values[k:k1]^-0.5, nrow = k1 - k + 1)
        B <- ev %*% B %*% ev
        eb <- eigen(B, symmetric = TRUE)
        rB <- eb$vectors %*% diag(sqrt(eb$values)) %*% t(eb$vectors)
        vec1 <- vec
        vec1[, k:k1] <- t(rB %*% diag(c(-1, 1)) %*% t(vec[, k:k1]))
        vec[, k:k1] <- t(rB %*% t(vec[, k:k1]))
    }
    else {
        vec1 <- vec <- if (k == 0)
            t(t(vec) * sqrt(1/ed$val[1]))
        else t(t(vec)/sqrt(ed$val[1:k]))
        if (k == 1)
            rank <- 1
    }
    d <- t(vec) %*% (R %*% p)
    d <- sum(d^2)
    d1 <- t(vec1) %*% (R %*% p)
    d1 <- sum(d1^2)
    rank1 <- rank
    if (nu > 0) {
        if (k1 == 1)
            rank1 <- val <- 1
        else {
            val <- rep(1, k1)
            rp <- nu + 1
            val[k] <- (rp + sqrt(rp * (2 - rp)))/2
            val[k1] <- (rp - val[k])
        }
        if (res.df <= 0)
            pval <- (mgcv::psum.chisq(d, val) + mgcv::psum.chisq(d1, val))/2
        else {
            k0 <- max(1, round(res.df))
            pval <- (mgcv::psum.chisq(0, c(val, -d/k0), df = c(rep(1,
                length(val)), k0)) + mgcv::psum.chisq(0, c(val, -d1/k0),
                df = c(rep(1, length(val)), k0)))/2
        }
    }
    else {
        pval <- 2
    }
    if (pval > 1) {
        if (res.df <= 0)
            pval <- (pchisq(d, df = rank1, lower.tail = FALSE) +
                pchisq(d1, df = rank1, lower.tail = FALSE))/2
        else pval <- (pf(d/rank1, rank1, res.df, lower.tail = FALSE) +
            pf(d1/rank1, rank1, res.df, lower.tail = FALSE))/2
    }
    list(stat = d, pval = min(1, pval), rank = rank)
}
## End testStat  (Wood, Biometrika 2013)
