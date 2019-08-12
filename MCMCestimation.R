mcmcS <- function(X,Y,alpha.p=2,beta.p=2,r.g=5,alpha.g=200,n.iter=500){
    #### parameters
    n.A <- nrow(X)   #number of potential curves
    n.T <- length(Y) #length of the process
    sum.trace <- 1:n.iter #evolution of the error to see the convergence
    #### Initialization of the algo 
    AA <- matrix(0,n.iter,n.A)
    AA[1,] <- rbinom(n.A,1,0.25)
    pp <- rep(0.3, n.iter)
    gg <- rep(0, n.iter)
    sum.A <- rep(0,ncol(X))
    for(k in 1:n.A){
        sum.A <- sum.A + AA[1,k]*X[k,]
    }
    RMSS <- mean((sum.A-Y)^2)
    gg[1] <- 1/RMSS
    #####let us run MCMC
    for(k in 2:n.iter){
        print(k)
        idx <- sample(1:n.A) ##sample the curves
        for(j in idx){
            sum.B <- sum.A-AA[k-1,j]*X[j,]
            sum.C <- sum.B + X[j,]
            w0 <- -gg[k-1]*sum((Y-sum.B)^2)
            w1 <- -gg[k-1]*sum((Y-sum.C)^2)
            prb <-  1/(1+(1-pp[k-1])/pp[k-1]*exp(w0-w1))
            AA[k,j] <- (runif(1) < prb)
            sum.A <- sum.A + (AA[k,j] - AA[k-1,j]) * X[j,]
        }
        ssA <- sum(AA[k,]) ##number of curves involve
        pp[k] <- rbeta(1, ssA+alpha.p, (n.A-ssA)+beta.p)
        RSS  <- sum((sum.A-Y)^2)#/2
        gg[k] <- rgamma(1, r.g + n.T/2, alpha.g + RSS )
        sum.trace[k] <- RSS
}
    return(list(tr=sum.trace,COEF=AA,yfit=sum.A,N=ssA))
}

baseline <- function(X,A,nn=10){
    n.A <- nrow(X)
    las <- nrow(A)
    RES <- matrix(0,nrow=nn,ncol(X))
    for(k in 1:nn){
        sel <- which(A[las-nn+k,]==1)
        RES[k,] <- apply(X[sel,],2,sum)
       }
    return(apply(RES,2,mean))
}

IC <- function(X,A){
    RES <- matrix(0,nrow=100,ncol=ncol(X))
    las <- nrow(A)
    n.A <- nrow(X)
    for(j in 1:100){
        sum.A <- rep(0,ncol(X))
        for(k in 1:n.A){
            sum.A <- sum.A + A[las-100+j,k]*X[k,]
        }
       RES[j,] <- sum.A
   }
    return(IC=RES)
}

IC2 <- function(X,A){
    RES <- matrix(0,nrow=100,ncol=ncol(X))
    las <- nrow(A)
    n.A <- nrow(X)
    for(j in 1:100){
        sel <- which(A[las-100+j,]==1)
        RES[j,] <- apply(X[sel,],2,sum)
   }
    return(IC=RES)
}    

sequential <- function(X,Ym,p=2,m=2){
    #### parameters
    #Ym the mean of the DR program
    n.T <- length(Ym) #length of the process
    YY <- matrix(Ym,byrow=T,nrow=nrow(X),ncol=n.T)
    n.A <- round(nrow(X)/m)
    selected <- NULL
    err <- rep(0,n.A)
    #### Initialization
    res <- apply(abs(X-YY)^p,1,sum)
    who <- names(which.min(res))
    selected <- c(selected,who)
    err[1] <- min(res)
    #####
    #iter <- 2
    #repeat{
    for(iter in 2:n.A){
          print(iter)
          ind <- which(rownames(X)%in%selected)
          XX <- X[-ind,]/iter
          if(length(selected)==1) Xs <- matrix(X[ind,],nrow=1) else Xs <- X[ind,]
          Z <- apply(Xs,2,mean)
          Z <- Z/iter*(iter-1)
          YY <- matrix(Ym-Z,byrow=T,nrow=nrow(XX),ncol=n.T)
          res <- apply(abs(XX-YY)^p,1,sum)
          err[iter] <- min(res)
          #if(err>err1){break}
          who <- names(which.min(res))
          selected <- c(selected,who)
          #iter <- iter+1
          #err1 <- err
          res <- which.min(err)
      }
    return(list(tr=err,sel=selected[1:res]))
}

