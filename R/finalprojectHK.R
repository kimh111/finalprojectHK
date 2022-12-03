#' @title Rounding to a specified decimal point
#' @description \code{rounding} Rounding to a specified decimal point
#'
#' @param x numeric value which want to be rounded
#' @param d decimal points to be rounded
#' @export
#' @examples
#' # rounding(2.55, 1)
#'
rounding <- function(x, d) {
  f <- 10^d
  a <- (x*f*2+1) %/% 2 / f
  return(a)
}




#' @title Huber function
#' @description \code{huber} output Huber function value from x axis value
#'
#' @param x numeric value of x axis
#' @export
#' @examples
#' # huber(1)
#' # huber(c(1,2,3))
huber <- function(x){
  if(x >= -1 & x <= 1) {
    y <- x^2
    return(y)
  }

  else {
    y <- 2*abs(x)-1
    return(y)
  }
}
huber <- Vectorize(huber)





#' @title type I error rate control
#' @description \code{error} output type I error rate control of Wald, Score, and Likelihood ratio test for Poisson distribution with specified sample size with bootstrapping method
#'
#' @param alpha Nominal type I error rate α
#' @param lambda0 λ0 value in Poisson distribution
#' @param B Number of repeats
#' @param n Limited sample size
#' @export
#' @examples
#' # error(0.05, 5, 100, 20)
#' # error(0.05, 10, 1000, 2000)
error <- function(alpha,lambda0,B,n){

  counts <- 0
  s.counts <- 0
  l.counts <- 0

  for(b in 1:B){
    set.seed(b)

    lambda <- mean(rpois(n, lambda=lambda0))

    w <-(lambda-lambda0)^2*n/lambda #statistical value of wald test
    w.pvalue <- 1-pchisq(w,1)

    s <- n*lambda/lambda0 - n #statistical value of score test
    ss <- s^2 * (lambda0/n)
    s.pvalue <- 1-pchisq(ss,1)

    l <- n * (lambda * log(lambda) - lambda) #statistical value of likelihood ratio test
    l0 <- n * (lambda * log(lambda0) - lambda0)
    ll <- -2*(l0 - l)
    l.pvalue <- 1-pchisq(ll,1)

    if(w.pvalue < alpha){
      counts <- counts + 1
    }

    if(s.pvalue < alpha){
      s.counts <- s.counts + 1
    }

    if(l.pvalue < alpha){
      l.counts <- l.counts + 1
    }
  }

  wald.test <-counts/B
  Score.test <-s.counts/B
  likelihood.ratio.test <-l.counts/B

  return(c(wald.test, Score.test, likelihood.ratio.test))
}




#' @title Statistical power of Wald, Score, and Likelihood ratio test in specified sample size
#' @description \code{power} Output statistical power of Wald, Score, and Likelihood ratio test for Poisson distribution with specified sample size with bootstrapping method
#'
#' @param alpha Nominal type I error rate α
#' @param lambda0 λ0 value in Poisson distribution
#' @param lambda1 λ1 value in Poisson distribution
#' @param B Number of repeats
#' @param n Limited sample size
#' @export
#' @examples
#' # power(0.05, 4, 5, 100, 20)
#' # power(0.05, 9, 10, 1000, 2000)
power <- function(alpha,lambda0,lambda1,B,n){

  counts <- 0
  s.counts <- 0
  l.counts <- 0

  for(b in 1:B){
    set.seed(b)

    lambda <- mean(rpois(n, lambda=lambda1))

    w <-(lambda-lambda0)^2*n/lambda #statistical value of wald test
    w.pvalue <- 1-pchisq(w,1)

    s <- n*lambda/lambda0 - n #statistical value of score test
    ss <- s^2 * (lambda0/n)
    s.pvalue <- 1-pchisq(ss,1)

    l <- n * (lambda * log(lambda) - lambda) #statistical value of likelihood ratio test
    l0 <- n * (lambda * log(lambda0) - lambda0)
    ll <- -2*(l0 - l)
    l.pvalue <- 1-pchisq(ll,1)


    if(w.pvalue < alpha){
      counts <- counts + 1
    }

    if(s.pvalue < alpha){
      s.counts <- s.counts + 1
    }

    if(l.pvalue < alpha){
      l.counts <- l.counts + 1
    }
  }

  wald.test <-counts/B
  Score.test <-s.counts/B
  likelihood.ratio.test <-l.counts/B

  return(c(wald.test, Score.test, likelihood.ratio.test))
}


