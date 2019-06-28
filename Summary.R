# -------------------- Topics --------------------
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
      plot(scores.data[,1:2])
      abline(h=0, v=0, lty=2, col='grey')
      
      # Scores with old variable directions
      x11()
      biplot(pc.data)
  
  
  
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
        M   <- sapply(x,mean)
        S    <- cov(x)
        Sinv <- solve(x.cov)
      # Test
        # T2 Statistics  n(x-mu)'*S2
        x.T2 <- n * (M-mu0) %*% Sinv %*% (M-mu0) 
        # Radius of the ellipsoid
        cfr.fisher <- ((n-1)*p/(n-p)) * qf(1-alpha,p,n-p)
        # Test: 
        x.T2 < cfr.fisher   # Rejection region: x.T2>cfr.fisher
        # P-value
        P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
      # Visualizing ellipse
        # Note that (S/n)^-1 = n*S^-1
        ellipse(M, S/n, sqrt(cfr.fisher), col = 'red' , lty = 2, lwd=2, center.cex=1) # CR
        ellipse(   mu0, S/n, sqrt(cfr.fisher), col = 'blue', lty = 2, lwd=2, center.cex=1) # AR
        # ellipse(mean, cov, radius, graphics...)
      # CR vs AR
      # Confidence Region is centerd in mu, Acceptance Region is centered in mu0, otherwise same.
        # Radius of the ellipsoid
        cfr.fisher <- ((n-1)*p/(n-p)) * qf(1-alpha,p,n-p)
        r <- sqrt(cfr.fisher)
        # Center 
        M
        # Direction of axises
        eigen(S)$vec
        # Length of semi
        r*sqrt(eigen(S/n)$values) 
      
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
          c1 <- a%*%M
          l1 <- a%*%M-sqrt(cfr.fisher*(t(a)%*%S%*%a)/n)
          r1 <- a%*%M+sqrt(cfr.fisher*(t(a)%*%S%*%a)/n)
          T2 <- cbind(inf = l1, center = c1, sup = r1)
          T2
    
      # Bonferroni CIs 
      # Note that for small k, BfCI is smaller than T2 but if k is large the contrary might be true
        # The 100(1-alpha) % cofidence inteerval that holds simoultaneously for the m following directions:
          a <- c(0,1) # There are m number of a but just as an example we have a <- c(1,0) and c(0,1) 
          k <- 2
          cfr.t <- qt(1 - alpha/(k*2), n-1)
          c1 <- a%*%M
          l1 <- a%*%M - cfr.t*sqrt(t(a)%*%S%*%a/n)
          r1 <- a%*%M + cfr.t*sqrt(t(a)%*%S%*%a/n)
          T2 <- cbind(inf = l1, center = c1, sup = r1)
          T2
        
    # Without Gaussian assumption simpultaneous CIs crf.fisher becomes Chi^2:
        cfr.nonG <- qchisq(1-alpha,p)
          
    # Common CI: (for each varible) where CI = M ± deltaCI
        # Bonferoni 
        deltaCI <- cfr.t*sqrt(diag(S)/n)
        # Sim T2
        deltaCI <- sqrt(cfr.fisher*diag(S)/n)
        
    # If we want to do multiple a at once we can do the following. 
      # Note "diag" in order to only get 3 values and not a 3 by 3 matrix
        A <- rbind(c(1,0),c(0,1),c(-1,1))
        c1 <- A%*%M
        l1 <- A%*%M-sqrt(cfr.fisher*diag(A%*%S%*%t(A))/n)
        r1 <- A%*%M+sqrt(cfr.fisher*diag(A%*%S%*%t(A))/n)
        T2 <- cbind(inf = l1, center = c1, sup = r1)
        T2
        # NOTE also that p_new = length(A%*%M), keep in mind if calculating cfr.fisher
        
        
  # Test mean on paired data (diff)----
    # The data in this case consists of D which is the difference of the paired data. 
    # Assuming Gaussianity 
      mcshapiro.test(D)
    # Test on H0: delta = 0 (no difference) vs ..
      n <- dim(D)[1]
      p <- dim(D)[2]
      dM   <- sapply(D,mean)
      dS    <-  cov(D)
      dSinv <- solve(dS)
      alpha   <- .05
      delta.0 <- rep(0,p)
      D.T2 <- n * (dM-delta.0) %*% dSinv %*% (dM-delta.0)
      cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
      D.T2 < cfr.fisher # We reject/accept at level alpha 
      P <- 1-pf(D.T2*(n-p)/((n-1)*p), p, n-p)
      # "worst" i.e. the direction in which delta0 is furthest away from the interval is:
      worst <- D.invcov %*% (D.mean-delta.0)
      worst <- worst/sqrt(sum(worst^2))
      worst # but why.....?
      theta.worst <- atan(worst[2]/worst[1])+pi
      
    # General case where we work with X as usual
      # Contrast matrix: in this case: x2-x1, x3-x1, x4-x1 
      C <- matrix(c(-1, 1, 0, 0,
                    -1, 0, 1, 0,
                    -1, 0, 0, 1), 3, 4, byrow=T)
      delta.0 = c(0,0,0) # Hypothesis!! (q = 3 differences in this case )
      dM <- C%*%colMeans(X)
      dS <- C%*%cov(X)%*%t(C)
      dSinv <- solve(dS)
      q <- length(dM)
      dT2 <- n*t(dM-delta.0)%*%dSinv%*%(dM-delta.0)
      cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha,(q-1),n-(q-1)) # Changed because of C
      dT2 < cfr.fisher # Test!
      P <- 1-pf(dT2*(n-q)/((n-1)*q), q, n-q)
      
  # 2 different populations (g1 and g2) with same sigma (possibly different n)
      # Test: H0: mu.1-mu.2==0 vs H1: mu.1-mu.2!=0
      p  <- 2 # or other obviously
      n1 <- dim(g1)[1]
      n2 <- dim(g2)[1]
      alpha <- 0.10
      
      mean1 <- sapply(g1,mean)
      mean2 <- sapply(g2,mean)
      cov1  <-  cov(g1)
      cov2  <-  cov(g2)
      Sp      <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)
      
      delta.0 <- c(0,0)
      Spinv   <- solve(Sp)
      
      T2 <- n1*n2/(n1+n2) * (mean1-mean2-delta.0) %*% Spinv %*% (mean1-mean2-delta.0)
      
      cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
      
      pvalue <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
      pvalue
    # with lin comb:
      dm <- (mean1-mean2)
      A  <- rbind(c(1,0), c(0,1), c(1,1), c(1,-1))
      k  <- dim(A)[1]
      
      A.s2 <- diag(A%*%Sp%*%t(A))
      A.dm <- A%*%(mean1-mean2)
      
      Bonf <- cbind(inf=A.dm - qt(1-(alpha/(2*k)), n1+n2-2) * sqrt( A.s2*(1/n1+1/n2) ), 
                    center=A.dm, 
                    sup=A.dm + qt(1-(alpha/(2*k)), n1+n2-2) * sqrt( A.s2*(1/n1+1/n2) ))
      Bonf
      
      
  # 2 populations -----
    # With probability 1-alpha the following will cover t(a)*(mu1-mu2)  
    # t(a)*(X1-X2) ± * c * qf(1-alpha, p, n1+n2-p-1) * sqrt(t(a)*(1/n1 + 1/n2)*Sp*a), in particualar
    # t(a)(X1_i-X2_i)                                            (1/n1 + 1/n2)*s_ii)
    # for i = 1...p
    # where Sp <- (n1 - 1)/(n1 + n2 -2)*S1 + (n2 - 1)/(n1 + n2 -2)*S2 
    # and c = (n1 + n2 -1)*p/(n1 + n2 - p - 1)
      # Other comment
    # to compare  mean(X1)-mean(X2) ( two population ) we estimate sigma (witch is assumed to be the same in both populations)
    # Sigma is estimated by (1/n1 + 1/n2) * Spooled
    # Spooled 
      
      
      mcshapiro.test(G1)
      mcshapiro.test(G2)
      
      # Test for independent Gaussian populations
      m1 <- sapply(G1,mean)
      m2 <- sapply(G2,mean)
      S1  <-  cov(G1)
      S2  <-  cov(G2)
      Sp      <- ((n1-1)*S1 + (n2-1)*S2)/(n1+n2-2)
      
      # Test: H0: mu.1-mu.2==0 vs H1: mu.1-mu.2!=0
      delta.0 <- c(0,0)
      Spinv   <- solve(Sp)
      T2 <- n1*n2/(n1+n2) * (m1-m2-delta.0) %*% Spinv %*% (m1-m2-delta.0)
      cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
      P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
      P
      
      
      # Bonferroni conf interval for difference in mean 
      k <- 2
      df <- n1+n2-2
      alpha <- 0.1
      IC <- cbind(m2-m1 - sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k*2), df),
                  m2-m1,
                  m2-m1 + sqrt(diag(Sp)*(1/n1+1/n2)) * qt(1 - alpha/(k*2), df))
      
      
      # Bonferroni conf interval var
      Bfvar <- cbind(df*diag(Sp)/qchisq(1 - alpha / (2*k), df),
                  df*diag(Sp),
                  df*diag(Sp)/qchisq(1 - alpha/(p*2), df))
      
      
      # Bonferroni for all for means (two in G1 and 2 in G2)
      k <- 4
      df <- n1+n2-2
      alpha <- 0.1
      IC <- cbind(c(m1,m2) - sqrt(diag(Sp)) * qt(1 - alpha/(k*2), df),
                  c(m1,m2),
                  c(m1,m2) + sqrt(diag(Sp))* qt(1 - alpha/(k*2), df))
      IC
      
      
  # One - way ANOVA ----
    # Test:
    ### Model: weigth.ij = mu + tau.i + eps.ij; eps.ij~N(0,sigma^2)
    ### H0: tau.1 = tau.2 = tau.3 = tau.4 = tau.5 = tau.6 = 0
    ### H1: (H0)^c
    ### We reject H0 if (SStr/g-1)/(SSres/n-g) is large. Distr as F(1-alpha, g-1,n-g)
    ### because SStr and SSres are both sum of square normdist <=> Chi^2
      
      # Setteing up enviroment data = c(X = data points, cat = categories)
        n       <- length(cat)      # total number of obs.
        ng      <- table(cat)       # number of obs. in each group
        treat    <- levels(cat)      # levels of the treatment
        g       <- length(types)     # number of levels (i.e., of groups)
      
      # Testing assumptions
        # Gaussianity in each group
        Ps <- NULL
        for(i in 1:g) Ps <- c(Ps,c(shapiro.test(X[ cat==treat[i] ])$p))
        
        # Same covariance structure (= same sigma^2)
        Var <- NULL
        for(i in 1:g) Var <- c(Var,c(var(X[ cat==treat[i] ])))
        # Test of homogeneity of variances
        # H0: sigma.1 = sigma.2 = sigma.3 = sigma.4 = sigma.5 = sigma.6 
        # H1: there exist i,j s.t. sigma.i!=sigma.j
        bartlett.test(X, cat)
      
      # One way ANOVA
        fit <- aov(X ~ cat)
        summary(fit)
        SSres <- sum(fit$residuals^2) # = W
        
      # Bonferroni intervalls for all combinations:
      # t_ik-t_jk = xbar_ik-xbar_jk ~ N(t_ik-t_jk, sigma_kk(1/n1 + 1/n2))
      # w/(n-g) is an estimator for sigma_kk 
        k <- g*(g-1)/2 
        S <- SSres/(n-g)
        M   <- mean(X)
        Mg  <- tapply(X, cat, mean)
        ICrange=NULL
        for(i in 1:(g-1)) {
          for(j in (i+1):g) {
            print(paste(cat[i],"-",cat[j]))        
            print(as.numeric(c(Mg[i]-Mg[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                               Mg[i]-Mg[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
            ICrange=rbind(ICrange,as.numeric(c(Mg[i]-Mg[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                                               Mg[i]-Mg[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
          }}
  # One - way MANOVA ----
        # prepare data
        data <- client[,2:3]
        groups <- client[,1]
        group.levels <- levels((groups))
        p <- dim(data)[2]
        g <- length(group.levels)
        
        
        i1 <- which(groups == group.levels[1])
        i2 <- which(groups == group.levels[2])
        i3 <- which(groups == group.levels[3])
        # ...
        ng <- c(length(i1),length(i2),length(i3)) # ...
        N <- sum(ng)
        # One-way MANOVA
        ### Model: X.ij = mu + tau.i + eps.ij; eps.ij~N_p(0,Sigma), [p=2]
        ###       X.ij, mu, tau.i in R^2, i=1,2,3 
        
        # 1) Test Gaussianity
        for(i in 1:g) message(mcshapiro.test(var.risp[ which(groups == group.levels[i]), ])$p)
        
        # 2) homogeneity in variance
        S1 <-  cov(data[ i1, ])
        S2 <-  cov(data[ i2, ])
        S3 <-  cov(data[ i3, ])
        
        x11()
        par(mfrow=c(1,3))
        image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
        image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
        image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
        
        
        # 3) fit the model
        fit <- manova(as.matrix(data) ~ groups)
        summary.manova(fit,test="Wilks")
        # Here we are looking for a small P - value so that we reject H0.
        
        # In order to investigate what effected the test above 
        # we look at the anova for each variable (p number of anovas)
        summary.aov(fit,test="Wilks")
        
        
        # Bonferroni CIs
        alpha <- 0.10
        k <- p*g*(g-1)/2
        qT <- qt(1-alpha/(2*k), N-g)
        W <- diag(t(fit$res) %*% fit$res)/(N-g)   
        
        m1 <- colMeans(data[i1,]) # mean for g1
        m2 <- colMeans(data[i2,]) # mean for g2
        m3 <- colMeans(data[i3,]) # ...
        
        qT <- qt(1 -alpha/(2*k), N-g)
        
        # Contrast BFs
        Bf12 <- cbind(m1-m2 - qT * sqrt((1/ng[1]+1/ng[2])*W), m1-m2, m1-m2 + qT * sqrt((1/ng[1]+1/ng[2])*W))
        Bf23 <- cbind(m2-m3 - qT * sqrt((1/ng[2]+1/ng[3])*W), m2-m3, m2-m3 + qT * sqrt((1/ng[2]+1/ng[3])*W))
        Bf31 <- cbind(m3-m1 - qT * sqrt((1/ng[3]+1/ng[1])*W), m3-m1, m3-m1 + qT * sqrt((1/ng[3]+1/ng[1])*W))
  # Two - way ANOVA (Alltså två stycken gruppeingar) + awesome bonferroni för anova ----
        ### Two-ways ANOVA
        ### Model without interaction (additive model): 
        ### X.ijk = mu + tau.i + beta.j + eps.ijk; eps.ijk~N(0,sigma^2), 
        ###     i=1,2 (effect direction centre-aero/aero-centre), 
        ###     j=1,2 (effect day weekday/weekend)
        # data is in this example a N x 3 matrix, column 2 and 3 corresponds to groups called AR and FF
        data <- euros
        names(data) <- c('X','G','B')
        attach(data)
        glev <- levels(G)
        blev <- levels(B)
        g <- length(glev) # nmbr of groups in first
        b <- length(blev) # nmbr of groups in second
        p <- 1 # dim = 1 since we are doing ANOVA
        n <- length(X[ G==glev[1] & B == blev[1] ]) # in this example n_ij = n_kl for all groupcombos
        N <- n*g*b 
        
        # 1) normality (univariate) in each group
        Ps <- c(shapiro.test(X[ G==glev[1] & B == blev[1] ])$p,
                shapiro.test(X[ G==glev[1] & B == blev[2] ])$p,
                shapiro.test(X[ G==glev[2] & B == blev[1] ])$p,
                shapiro.test(X[ G==glev[2] & B == blev[2] ])$p)
        Ps
        
        for(i in blev){
          for (j in glev) {
            message(i,j,shapiro.test(X[ G==j & B == i ])$p," ",length(X[ G==j & B == i ]))
          }
        }
        # 2) homogeneity of variances
        bartlett.test(list(X[ G==glev[1] & B == blev[1] ],
                           X[ G==glev[1] & B == blev[2] ],
                           X[ G==glev[2] & B == blev[1] ],
                           X[ G==glev[2] & B == blev[2] ]))
        
        # 3) Fit the model:
        fit <- aov(X ~ B + G + B:G)
        summary(fit)
        
        # Estimate variances (for everything) 
        W <- sum(fit$residuals^2)  # SS_res
        var <- W/(g*b*n-g-b+1)     # SS_res/gdl(res)
        
        # Estimate the great mean mu:
        m <- mean(X)
        
        # Estimate tau.i, beta.j:
        tau1  <- mean(data[G==glev[1],1]) - m  # tau.1
        tau2  <- mean(data[G==glev[2],1]) - m  # tau.2
        
        beta1 <- mean(data[B==blev[1],1]) - m  # beta.1
        beta2 <- mean(data[B==blev[2],1]) - m  # beta.2
        
        # Point-wise estimates of mean, m_ij, i in G and j in B
        # (model without interaction!)
        m11 <- m + tau1 + beta1
        m12  <- m + tau1 + beta2
        m21 <- m + tau2 + beta1
        m22  <- m + tau2 + beta2
        
        # Estimate variances (for everything) 
        W <- sum(fit$residuals^2)  # SS_res
        Sp <- W/(N-b)     # SS_res/gdl(res)
        m1 <- mean(X[B == blev[1]])
        m2 <- mean(X[B == blev[2]])
        m3 <- mean(X[B == blev[3]])
        n1 = length(X[B == blev[1]])
        n2 = length(X[B == blev[2]])
        n3 = length(X[B == blev[3]])
        # Bonferroni conf interval mean 
        k <- 3
        df <- N-3
        alpha <- 0.01
        IC <- cbind(c(-m1+m2,-m2+m3,m1-m3) - c(sqrt(Sp*(1/n1+1/n2)),sqrt(Sp*(1/n2+1/n3)),sqrt(Sp*(1/n1+1/n3))) * qt(1 - alpha/(k*2), df),
                    c(-m1+m2,-m2+m3,m1-m3),
                    c(-m1+m2,-m2+m3,m1-m3) + c(sqrt(Sp*(1/n1+1/n2)),sqrt(Sp*(1/n2+1/n3)),sqrt(Sp*(1/n1+1/n3))) * qt(1 - alpha/(k*2), df))
        IC
  # LDA/QDA (Classifiers)----
    library(MASS)
    # Assumptions:
        # 1) X.i ~ N(mu.i, sigma.i^2), i=A,B (label)
        # 2) sigma.A=sigma.B (only LDA)
        # 3) c(A|B)=c(B|A) (equal misclassification costs)
      # verify assumptions 1) e 2): 
      # 1) normality (univariate) within the groups
      shapiro.test(G1)
      shapiro.test(G2)
      
      # 2.1) Test variance  (univariate) 2 groups
      var.test(G1, G2)
      # 2.2) equal variance (univariate) >2 groups
      bartlett.test(values ~ groups)
      # 2.3) Test by visualizing (multivariate)
      S1 <- cov(G1)
      S2 <- cov(G2)
      x11()
      par(mfrow=c(1,2))
      image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
      image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE, breaks = quantile(rbind(S1,S2), (0:100)/100, na.rm=TRUE))
      
      # 3) LDA/QDA, note that for LDA the test above has to hold
      # If we have a cost of missclassification we use the following priors
      ct1p2 <- 10 # cost of error: true = 1, predicted 2
      ct2p1 <- 1  # cost of error: true = 2, predicted 1
      prior.c <- c(prior[1]*ct2p1,prior[2]*ct1p2)/(ct1p2*prior[2]+ct2p1*prior[1])
      
      
      fit <- lda(X, Xg) #, CV = T) for cross-validation
      x  <- seq(min(X[,1]), max(X[,1]), length=200)
      y  <- seq(min(X[,2]), max(X[,2]), length=200)
      xy <- expand.grid(X1=x, X2=y)
      LDA  <- predict(fit, xy)
      
      # LDA$posterior returns the posterior probability of x belonging to each group
      
      # Visualization
      x11()
      plot(X, main='Data', xlab='X1', ylab='X2', pch=20)
      points(G1, col='red', pch=20)
      points(G2, col='blue', pch=20)
      legend('bottomleft', legend=levels(Xg), fill=c('blue','red'), cex=.7)
      
      points(fit$means, pch=4,col=c('blue','red') , lwd=2, cex=1.5)
      
      lines(LDA$post)
      
      z  <- LDA$post  
      
      z1 <- z[,1] - z[,2] 
      z2 <- z[,2] - z[,1]  
      contour(x, y, matrix(z1, 200), levels=0, drawlabels=F, add=T)  
      contour(x, y, matrix(z2, 200), levels=0, drawlabels=F, add=T)
      
      
      
  # Knn ----
      # Returns a class which it would assign x to
      output <- knn(train = X, test = x, cl = X.class, k = k)
  # Error rates ----
      misc <- table(trueC = X.class, predictedC=predict(fit)$class)
      # AER (Actuall Error Rate) 
      # AER = p1*n12/(n11+n12) + p2*n21/(21+22) (confusion matrix aka misc)
      # for 2 x 2: 
      AER <- p1*misc[1,2]/(misc[1,1]+misc[1,2]) + p2*misc[2,1]/(misc[2,1]+misc[2,2])
      # Analytically calculated AER
      misclass <- function(x){
        prior[1]*(1 - pnorm(x, M1, SD)) + prior[2]*pnorm(x, M2, SD)
      }
      AER <- optimize(f=misclass, lower=min(X), upper=max(X))$objective
      # APER (Apperent Error Rate)
      # In confusion matrix (misc), sum of errrors divided by number of players
      # (n12 + n21)/n
      # fit <- lda() or qda() or similar
      G <- 2
      APER <- 0
      for(g in 1:G)
        APER <- APER + sum(misc[g,-g])/sum(sum(misc[g,])) * prior[g]
  # SVM ----
    library(e1071)
    # dat is a data.frame with one column named y 
    # cost and gamma are tuning parameters
    fitlin = svm(y~., data=dat, kernel ='linear', cost = 10, scale =FALSE )
    fitrad = svm(y~., data=dat, kernel ='radial', cost = 1 , gamma =1)
  # Clustering ----
    # Hierarchical clustering
    dat.dist <- dist(iris4, method='euclidean')# or 'manhattan' or 'canberra'
    dat.clust <- hclust(dat.dist, method='single')# or 'average' or 'complete' OR 'ward.D2' 
    # single = smallest dist between 2 groups, complete = largest, ...
    # average - obvious and 'ward.D2' minimizes loss of information when merging 2 sets
    plot(dat.clust, main='Hierarchical clustering', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
    cluster <- cutree(dat.clust, k=2) # Choose k
    table(cluster) #how many in each cluster?
    coph <- cophenetic(dat.clust)# compute the cophenetic matrices 
    
    # compute cophenetic coefficients (the closer to one the better)
    cc <- cor(dat.dist, coph)# compute cophenetic coefficient(s...?)
    
    # Visualization
    x11()
    plot(dat, col = as.vector(cluster)+1, pch=19)
  # K - means clustering ----
    # to choose k:
    # 1) evaluate the variability between the groups w.r.t. the variability 
    #    withing the groups - plot for different values of k and choose one w.r.t 
    #    large gaps between values of k (tänk typ hierarchisk cluster)
    # 2) evaluate the result of hierarchical clustering (not recommended,
    #    quite computationally expensive)
    
    result.k <- kmeans(dat, centers=2) # Centers: fixed number of clusters
    
    names(result.k)
    
    result.k$cluster      # labels of clusters
    result.k$centers      # centers of the clusters
    result.k$totss        # tot. sum of squares
    result.k$withinss     # sum of squares within clusters
    result.k$tot.withinss # sum(sum of squares within cluster)
    result.k$betweenss    # sum of squares between clusters
    result.k$size         # dimention of the clusters
    
    # Visualize
    result.k <- kmeans(dat, k)
    
    x11()
    plot(Q, col = result.k$cluster+1)
    
    open3d()
    plot3d(dat, size=3, col=result.k$cluster+1, aspect = F) 
    points3d(result.k$centers, pch = 4, cex = 2, lwd = 4)
    
  # Linear model ----
    dat <- data.frame(x = NA,
                      y = NA,
                      z = NA,
                      s = NA)
    n <- length(x)
    # linmod: y = beta0 + beta1*x + beta2*y + beta3*z + beta4*x*y
    fit <- lm(s ~ x + y + z + x:y , data=dat)
    
    # linmod: 3 på varandra följande grupper 
    xG <- c(x,y,z)
    sG <- rep(s,3)
    # g levels, g-1 dummies
    datG <- data.frame(sG = sG, 
                       xG2 = xG^2, 
                       dumx = rep(c(1,0), c(n,2*n)),   # ones only on x
                       dumy = rep(c(0,1,0), c(n,n,n))) # ones only on y
    # Linmod: both of the following are equivalent
    # Y | (X = x, G = g) = beta0.g + beta1.g * x^2 + eps
    # Sg = b.0 + b.1*dumx + b.2*dumy + b.3*datG^2 + b.4*dumx*datG^2 + b.5*dumy*datG^2
    fit <- lm(sG ~ dumx + dumy + xG2 + dumx:xG2 + dumy:xG2, data=datG) 
    # note that z is not present, thus intercept for x, y, z would be:
    # coef(fit)[1] + coef(fit)[2], coef(fit)[1] + coef(fit)[3] and coef(fit)[1]   
    fit2 <- lm(sG ~ xG2 + xG2:dumx + xG2:dumy, data=datG) # reduced alternative
    
    var(fit$residuals) # sigma^2 estimate
    
    
    # Residual investigation:
    shapiro.test(residuals(fit)) # model assumptions!
    shapiro.test(rstudent(fit))  # model assumptions!
    x11()
    par(mfrow=c(2,2))
    plot(fit)
    
    # Linear hypothesis:
    # H0: A*beta = c, beta and c vectors, A matrix
    # ex: H0: c(b3,b4,b5) = c(2,3,4) vs ....
    A <- rbind(c(0,0,0,1,0,0),
               c(0,0,0,0,1,0),
               c(0,0,0,0,0,1))
    c <- c(2,3,4)
    linearHypothesis(fit,A,c) # Basically just look at p-value
    # the F-test would be c = (0,0,0)
    
    
    # Prediction: 
    predict(fit, newdata = ndat,interval = 'prediction', level = 1-alpha)
    # ndat is a data.frame with one value for each beta in fit
    # The last two objects are just if you want an prediction interval 
    # Example of prediction with dummy with a global 95% conf interval
    alpha <- 0.01
    k <- 2
    new <- data.frame(Vi2 = 25^2, Va2 = 35^2, dum=1)
    IP <- predict(fit.opt, newdata=new, interval='prediction', level=1-alpha/k)
    IP
    
    
# -------------------- Other --------------------
  # General function, good2know ----
    load("/Users/molsby/Documents/R/Applied_stat/fun/mcshapiro.test.RData")
      # Reading txt - files 
      data <- read.table('filename.txt', header=T)
      
      # d*dist* gives the density function (pdf)
      # p*dist* gives the cummulative distribution function (cdf)
      # q*dist* gives the quantile function aka for what value it covers p percent aka inverse of cdf 
      # r*dist* generates random deviates.
      # Examples: df, dchisq, dnorm 
      
      sapply(data,FUN = sd) # returns sd (or other) of each column
      data.st <- scale(data) # returns standardized data
      ones <- rep(1, n) # n isthe length
      
      # General and Total Variance
      var.gen <- det(S)
      var.tot <- sum( diag(S) )
      
      # Visualization functions
      coler.latter <- rainbow(13)
      
      #
  # Go2: Multivariate Test and CIs with lin-comb ----
    # General setting, data is given by X
      X <- as.matrix(dat)
      # Test Gaussianity
      mcshapiro.test(X)
      # Set up
      alpha <- 0.01
      n <- dim(X)[1]
      p <- dim(X)[2]
      A <- matrix(c(-1, 1, 0, 0,
                     0, -1, 1, 0,
                     0, 0, -1, 1), 3, 4, byrow=T)
      M.A <- A %*% sapply(X,mean) #colMeans(X)
      S.A <- A %*% cov(X) %*% t(A)
      Sinv.A <- solve(S.A)
      p.new <- length(M.A)
      
      # Test
      mu0 <- rep(0,p.new)
      T2.A <- n * t(M.A - mu0) %*% Sinv.A %*% (M.A - mu0)
      cfr.fisher.new <- ((n-1)*p.new/(n-p.new)) * qf(1-alpha,p.new,n-p.new)
      T2.A < cfr.fisher.new
      P <- 1-pf(T2.A * (n-p.new)/(p.new*(n-1)), p.new, n-p.new)
      P
      
      # T2 sim CI (mean)
      T2 <- cbind( "Inf"= M.A-sqrt(cfr.fisher.new*diag(S.A)/n), 
                   'Center' = M.A, 
                   'Sup'= M.A+sqrt(cfr.fisher.new*diag(S.A)/n))
      
      # Bonferroni CI
      k <- p.new
      cfr.t.new <- qt(1 - alpha/(k*2), n-1)
      BF <- cbind( "Inf"= M.A-cfr.t.new*sqrt(diag(S.A)/n), 
                   'Center' = M.A, 
                   'Sup'= M.A+cfr.t.new*sqrt(diag(S.A)/n))
      
  
  # Go2: Multivariate Test and CIs without lin comb ---- 
    # General setting, data is given by X
      #X <- as.matrix(data)
      # Test Gaussianity
      mcshapiro.test(X)
      # Set up
      alpha <- 0.01
      n <- dim(X)[1]
      p <- dim(X)[2]
      M <- sapply(X,mean)
      S <- cov(X)
      Sinv <- solve(S)
      
      # Test
      mu0 <- rep(0,p)
      x.T2 <- n * t(M - mu0) %*% Sinv %*% (M - mu0)
      cfr.fisher <- ((n-1)*p/(n-p)) * qf(1-alpha,p,n-p)
      x.T2 < cfr.fisher
      P <- 1-pf(x.T2 * (n-p)/(p*(n-1)), p, n-p)
      P
      
      # T2 sim CI (mean)
      T2 <- cbind( "Inf"= M-sqrt(cfr.fisher*diag(S)/n), 
                   'Center' = M, 
                   'Sup'= M+sqrt(cfr.fisher*diag(S)/n))
      
      # Bonferroni CI
      k <- p.new
      cfr.t.new <- qt(1 - alpha/(k*2), n-1)
      BF <- cbind( "Inf"= M-cfr.t*sqrt(S/n), 
                   'Center' = M.A, 
                   'Sup'= M+cfr.t*sqrt(S/n))
      
      
      
  # Univariate t-test ----
      a <- c(-1,1,0) # CHange this
      alpha <- 0.02
      t.test(as.matrix(X)%*%a, alternative = 'two.sided', mu = delta.0, conf.level = 1-alpha)
      # Testa om 2 mu är olika från 2 data (H0: mu1 = mu2)
      t.test(x = group1, y = group2, var.eq=T)
  # Testing if variance is equal ----
      # Få till den här
      bartlett.test(Infg ~ group)
      # Och denna i univar
      var.test(group1,group2)
      
      
      
  # GOOD QUESTIONS----
      # LDA/QDA on multivariate with g >=3
  # Estimate an elliptical region containing 99% of X ----
      mcshapiro.test(X)
      
      M <- sapply(X, mean)
      S <- cov(X)
      alpha <- 0.01
      cfr.chisq <- qchisq(1-alpha,p)
      
      # Characterize the ellipse:
      # Axes directions:
      eigen(S)$vectors
      # Center:
      M
      # Radius of the ellipse:
      r <- sqrt(cfr.chisq)
      # Length of the semi-axes:
      r*sqrt(eigen(S)$values)  
      
      x11()
      plot(D, asp = 1, col='gold', pch=19, xlim=c(-10,50))
      ellipse(center=M, shape=S, radius=sqrt(cfr.chisq), col = 'black', lty = 2, center.pch = 4)
  
  
  # Geo ----
      # data is dat with 1 and 2 col x, y, then other data i.e. zinc 
      coordinates(dat) <- c('x','y')
      v <- variogram(log(zinc) ~ 1, dat) 
      #, alpha = c(0, 45, 90, 135)) or , cutoff = 1000, width = 1000/15) or 
      # , boundaries = c(0,200,seq(400,1500,100)))
      v.fit <- fit.variogram(v, vgm(1, "Sph", 800, 1)) # vgm(psill, model, range, nugget)
      vgm() 
      plot(v, v.fit, pch = 19)
      # error:
      attr(v.fit, 'SSErr')
      # example where dirr: aplha = 45 gives us the major range. 
      # i.e. the range in the direction with strongest correlation 
      v.dir <- variogram(log(zinc)~1,meuse,alpha=(0:3)*45) 
      v.anis <- vgm(.6, "Sph", 1600, .05, anis=c(45, 0.3))
      
      # kriging
      g.tr <- gstat(formula = log(zinc) ~ 1, data = meuse, model = v.fit)
      s0=as.data.frame(matrix(c(x0,y0),1,2))
      names(s0)=c('X','Y')
      coordinates(s0)=c('X','Y')
      # Make the ordinary kriging prediction with the function:
      predict(s0, g.tr, BLUE=FALSE) # Predicts a observatio given a distance
      predict(s0, g.tr, BLUE=TRUE)  # Predicts a mean given a distance
      
      # Or if we have a model dependent of distance for example:
      g.tr <- gstat(formula = zinc ~ D, data = meuse, model = v.fit)
      D.s0=# some value 
      s0=as.data.frame(matrix(c(x0,y0,D.s0),1,3))
      names(s0)=c('X','Y','D')
      coordinates(s0)=c('X','Y')
      
      #(a) Discontinuity at the origin: The process presents a highly irregular 
      #    behavior and it is not L2-continuous. Such a discontinuity occurs in 
      #    case of nugget effect, due to measurement errors, or microscale variability.
      
      #(b) Linear at the origin: Such a behavior is common in continuous but non- differentiable 
      #    processes. The linear, exponential and spherical models are instances of such structures.
      
      #(c) Quadratic at the origin: The associated random field is very regular with smooth 
      #    realizations, and often presents a drift term.
      
  # Ridge & Lasso ----
    # Ridge
      ls <- seq(0,10,0.01)
      fit.ridge <- lm.ridge(Y ~ X1 + X2 , lambda = ls) 
      lambda.opt <- ls[which.min(fit.ridge$GCV)]
    # Lasso
      library(glmnet)
      x <- model.matrix(Y~X1+X2)[,-1]
      # Build the vector of response
      y <- Y
      # Let's set a grid of candidate lambda's for the estimate
      lambda.grid <- 10^seq(10,-2,length=100)
      fit.lasso <- glmnet(x,y, lambda = lambda.grid) # default: alpha=1 -> lasso [note: if alpha=0 -> ridge regression]
      cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid) # with cross val
      bestlam2 <- cv.lasso$lambda.min
      predict(fit.lasso, s=bestlam2, type = 'coefficients')
      
      
      
      
      
      
  # Functional Data ----
      
  # Spooled ----
      # to compare  mean(X1)-mean(X2) ( two population ) we estimate sigma (witch is assumed to be the same in both populations)
      # Sigma is estimated by (1/n1 + 1/n2) * Spooled
      # Spooled 
      
      
      
      
      