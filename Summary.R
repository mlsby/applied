# General function, good2know ----
  # Reading txt - files 
  data <- read.table('filename.txt', header=T)
  
  sapply(data,FUN = sd) # returns sd (or other) of each column
  data.st <- scale(data) # returns standardized data
  ones <- rep(1, n) # n isthe length
  
  # General and Total Variance
  var.gen <- det(S)
  var.tot <- sum( diag(S) )
  
  # Visualization functions
  coler.latter <- rainbow(13)

#
# PCA ----
  # Computing PCA, data should be standardized
    pc.data <- princomp(data, scores = T)
    summary(pc.data) # Good first step
    pc.data$sd # standard div of components
    pc.data$loadings # loadings = the linear compination vector for each component
    cumsum(pc.data$sd^2)/sum(pc.data$sd^2) # cumulative proportion of explained variance
  
  # Vizualisation:
    # Explained variance (basically plots the cumsum)
    yl <- c(0,3*max(sapply(data,sd))^2)
    x11()
    layout(matrix(c(2,3,1,3),2,byrow=T))
    plot(pc.data, las=2, main='Principal components', ylim=yl)
    barplot(sapply(data,sd)^2, las=2, main='Original Variables', ylim=yl, ylab='Variances')
    plot(cumsum(pc.data$sd^2)/sum(pc.data$sd^2), type='b', axes=F, xlab='number of components', 
         ylab='contribution to the total variance', ylim=c(0,1))
    abline(h=1, col='blue')
    abline(h=0.8, lty=2, col='blue')
    box()
    axis(2,at=0:10/10,labels=0:10/10)
    axis(1,at=1:ncol(data),labels=1:ncol(data),las=2)
    
    # barplot of components where N is the number of selected PCs
    x11()
    par(mar = c(1,4,0,2), mfrow = c(N,1))
    for(i in 1:N) barplot(pc.data$loadings[,i], ylim = c(-1, 1))
    
    # Variability of the original variables / scores
    x11()
    layout(matrix(c(1,2),2))
    boxplot(data, las=2, col='gold', main='Original variables')
    scores.data <- data.frame(pc.data$scores)
    boxplot(scores.data, las=2, col='gold', main='Principal components')
    
    # Plotting scores = outcome after linear transformation in relation to the two first PCAs
    x11()
    plot(scores.tourists[,1:2])
    abline(h=0, v=0, lty=2, col='grey')
    
    # Scores with old variable directions
    x11()
    biplot(pc.tourists)



#
# Testing multivariate Gaussianity ----
  # Tests:
  shapiro.test(X[,1]) # (one-dim) H0: X~N, W close to 1 = gaussian, close to 0 is not.
  load("/Users/molsby/Documents/R/Applied_Stat/fun/mcshapiro.test.RData") # Two dim
  mcshapiro.test(X) # Simoultaineous test on X returning min(W.a)
  
  # Generating data 
  mu <- c(1,2)
  sig <- matrix(c(1,1,1,2), 2)
  n   <-  150
  X <- rmvnorm(n, mu, sig)
  
  # One-dimensional plots, for 2 dim X
  x <- X[,1]
  y <- X[,2]
  s <- seq(-2,6,len = 1000)
  x11()
  par(mfrow=c(2,2))
  
  hist(x, prob=T,col='grey85')
  lines(s, dnorm(s,mean(x),sd(x)), col='blue', lty=2)
  
  hist(y, prob=T,col='grey85')
  lines(s, dnorm(s,mean(y),sd(y)), col='blue', lty=2)
  
  qqnorm(x, main='QQplot of x')
  qqline(x)
  
  qqnorm(y, main='QQplot of y')
  qqline(y)


# Making data Gaussian ----
  # Methods:
  #  - Identifying clusters
  #  - Removing outliers
  #  - Transformation (e.g. Box - Cox)
  
  # UNIVARIATE
    # Choosing lambda
    # For lambda>1: observations <1 are "shrinked", observations >1 are "spread"
    # For lambda<1: observations <1 are "spread", observations >1 are "shrinked"
    lambda.X <- powerTransform(X) # Returns the optimal lambda
    
    # Box - Cox transformation 
    bc.X <- bcPower(X, lambda.X$lambda)
    
  # MULTIVARIATE 
    x <- X[,1]
    y <- X[,2]
    # Choosing lambda
    lambda <- powerTransform(X)
    
    # Box - Cox transformation 
    BC.x <- bcPower(x, lambda$lambda[1])
    BC.y <- bcPower(y, lambda$lambda[2]) # Obviously just keep going if dim > 2
    
    # Visualization of 2-dim data (X) with projection on each dim
    xl <- c(-10,10)
    yl <- c(-10,10)
    x11()
    plot(X,pch=19,main='Data',xlim=xl, ylim=yl)
    points(X[,1], rep(yl[1],dim(X)[1]), col='red', pch=19)
    points(rep(xl[1],dim(X)[1]), X[,2], col='blue', pch=19)
  
  # Summary: 
    # Box Cox is good as a start, if lambda is close to 1, 
    # consider the hypothesis of no transformation. On the other hand,
    # if lambda close to 0 consider the hypothesis of log-transformation.


# Probability Regions ----
  # c(x,y,z) = X ~ N(mu,X.Sigma)
  # Probability region A: P((x y)' \in A) = alpha
  # This region is obviusly an ellips with the following features
    alpha <- 0.9
    Sigma <- X.Sigma[1:2,1:2]
    # Direction of the axes:
    eigen(Sigma)$vectors
    # Center:
    M <- mu
    # Radius of the ellipse:
    r <- sqrt(qchisq(p = alpha, df = 2))
    # Length of the semi-axes:
    r*sqrt(eigen(Sigma)$values)
    
  # Conditional probability region A: P(X1' \in A2 | X2=x2) = alpha, X1 and X2 can be multi dim
    mu.cond <- function(mu1,mu2,Sig12,Sig22,x2){
      return(mu1+Sig12%*%solve(Sig22)%*%(x2-mu2))
    }
    Sig.cond <- function(Sig11,Sig12,Sig22){
      Sig21=t(Sig12)
    return(Sig11-Sig12%*%solve(Sig22)%*%Sig21)
    }
    # Example, c(x,y,z) = X ~ N(mu,X.Sigma), X1 = c(x,y), X2 = z, x2 = b
    M.c <- mu.cond(mu1=mu[1:2],mu2=mu[3],Sig11=Sigma[1:2,1:2],Sig12=Sigma[1:2,3],Sig22=Sigma[3,3],x2=b)
    Sigma.c <- Sig.cond(Sig11=Sigma[1:2,1:2],Sig12=Sigma[1:2,3],Sig22=Sigma[3,3])
    # Features of ellips is calculated as before....
      
    # Clearification of SigXX:
      # Sig11 = sigma of X1 X.Sigma[1:2,1:2]
      # Sig12 = correlation between X1 & X2 X.Sigma[1:2,3]
      # Sig22 = sigma of X2 X.Sigma[3,3]
      
    
    
# Confidence Regions ----
    # Test for the mean of a multivariate Gaussian with level alpha
    # H0: mu=mu0 vs H1: mu!=mu0
    # Assumption: Gaussianity
      mcshapiro.test(x) # test assumption
    # Precomputations on x
      alpha <- 0.01
      mu0 <- c(1,0)
      n <- dim(x)[1]
      p <- dim(x)[2]
      x.mean   <- sapply(x,mean)
      x.cov    <- cov(x)
      x.invcov <- solve(x.cov)
    # Test
      # T2 Statistics  n(x-mu)'*S2
      x.T2 <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
      # Radius of the ellipsoid
      cfr.fisher <- ((n-1)*p/(n-p)) * qf(1-alpha,p,n-p)
      # Test: 
      x.T2 < cfr.fisher   # Rejection region: x.T2>cfr.fisher
      # P-value
      P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
    # Visualizing ellipse
      S <- x.cov # Note that (S/n)^-1 = n*S^-1
      ellipse(x.mean, S/n, sqrt(cfr.fisher), col = 'red' , lty = 2, lwd=2, center.cex=1) # CR
      ellipse(   mu0, S/n, sqrt(cfr.fisher), col = 'blue', lty = 2, lwd=2, center.cex=1) # AR
      # ellipse(mean, cov, radius, graphics...)
    # CR vs AR
    # Confidence Region is centerd in mu, Acceptance Region is centered in mu0, otherwise same.
    
    # Test for the mean with level alpha without Gaussianity (Assymptotic)
    # H0: mu=mu0 vs H1: mu!=mu0
    # Assumption: n is large (n -> inf)
      # Radius of the ellipsoid
      cfr.chisq <- qchisq(1-alpha, p) # qchisq(probability, df)
      # Test:
      x.T2 < cfr.chisq
      # Compute the p-value
      PA <- 1-pchisq(x.T2A, p)
  
    # General rule to perform a test:
      ###
      # 1)  Formulate the test (and test the Gaussian assumption, if needed)
      # 2)  Compute the test statistics 
      # 3a) Having set the level of the test, verify whether the test statistics 
      #     belongs to the region of rejection (i.e., if there is statistical  
      #     evidence to reject H0)
      # 3b) Compute the p-value of the test
      ###
  
  # Under assumption of Gaussianity
    # Simultaneous CIs
      # The 100(1-alpha) % cofidence interval that holds simoultaneously for all directions:
        a <- c(1,0)
        c1 <- a%*%x.mean
        l1 <- a%*%x.mean-sqrt(cfr.fisher*(t(a)%*%x.cov%*%a)/n)
        r1 <- a%*%x.mean+sqrt(cfr.fisher*(t(a)%*%x.cov%*%a)/n)
        T2 <- cbind(inf = l1, center = c1, sup = r1)
        T2
  
    # Bonferroni CIs
      # The 100(1-alpha) % cofidence inteerval that holds simoultaneously for the m following directions:
        a <- c(0,1) # There are m number of a but just as an example we have a <- c(1,0) and c(0,1) 
        m <- 2
        cfr.t <- qt(1 - alpha/(m*2), n-1)
        c1 <- a%*%x.mean
        l1 <- a%*%x.mean - cfr.t*sqrt(t(a)%*%x.cov%*%a/n)
        r1 <- a%*%x.mean + cfr.t*sqrt(t(a)%*%x.cov%*%a/n)
        T2 <- cbind(inf = l1, center = c1, sup = r1)
        T2
      
  # Without Gaussian assumption simpultaneous CIs crf.fisher becomes:
      cfr.nonG <- qchisq(1-alpha,p)
        
      
      
      
# ----