\name{residuals}
\alias{residuals.gel}
\alias{residuals.gmm}
\title{Residuals of GEL or GMM}
\description{
 Method to extract the residuals of the model estimated by \code{gmm} or \code{gel}.
}
\usage{
\method{residuals}{gel}(object, ...)
\method{residuals}{gmm}(object, ...)
}
\arguments{
 \item{object}{An object of class \code{gmm} or \code{gel} returned by the function \code{\link{gmm}} or \code{\link{gel}}}
\item{...}{Other arguments when \code{residuals} is applied to an other classe object}
}

\value{
It returns the  matrix of residuals \eqn{(y-\hat{y})} in \code{g=y~x} as it is done by  \code{residuals.lm}.
}

\examples{

# GEL can deal with endogeneity problems

n = 200
phi<-c(.2,.7)
thet <- 0.2
sd <- .2
set.seed(123)
x <- matrix(arima.sim(n = n, list(order = c(2,0,1), ar = phi, ma = thet, sd = sd)), ncol = 1)

y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]
H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
g <- y ~ ym1 + ym2
x <- H

res <- gel(g, x, c(0,.3,.6))
e <- residuals(res)
plot(e, type = 'l', main = "Residuals from an ARMA fit using GEL")

# GMM is like GLS for linear models without endogeneity problems

set.seed(345)
n = 200
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n = n, list(order = c(2,0,1), ar = phi, ma = thet, sd = sd)), ncol = 1)
y <- 10 + 5*rnorm(n) + x

res <- gmm(y ~ x, x)
plot(x, residuals(res), main = "Residuals of an estimated model with GMM")


}

