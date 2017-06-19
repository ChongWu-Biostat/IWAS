### SSU #################

SumSqU<-function(U, CovS){
    if (is.null(dim(CovS))) {# only one-dim:
        Tscore<- sum(U^2 /CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
        pTg1<-as.numeric(1-pchisq(Tscore, 1))
    }
    else {
        #it's possible U=0 and Cov(U)=0:
        if (all(abs(U)<1e-20)) pTg1<-1 else{
            Tg1<- t(U) %*% U
            ##distr of Tg1 is sum of cr Chisq_1:
            cr<-eigen(CovS, only.values=TRUE)$values
            ##approximate the distri by alpha Chisq_d + beta:
            alpha1<-sum(cr*cr*cr)/sum(cr*cr)
            beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
            d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
            alpha1<-as.double(alpha1)
            beta1<-as.double(beta1)
            d1<-as.double(d1)
            pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))
        }
    }
    return(pTg1)
}

##########SumTest########################
Sum<-function(U, CovS){
    #it's possible U=0 and Cov(U)=0:
    if (all(abs(sum(U))<1e-20)) pTsum<-1 else{
        a<-rep(1, length(U))
        Tsum<- sum(U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
        pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )
    }
    pTsum
}

##########UminP Test########################
UminPd<-function(U, CovS){
    
    if (is.null(dim(CovS))) {# only one-dim:
        Tu<- sum(U^2 /CovS)
        if (is.na(Tu) || is.infinite(Tu) || is.nan(Tu)) Tu<-0
        pTu<-as.numeric(1-pchisq(Tu, 1))
    }
    else{
        ####it's POSSIBLR Ui=0 and CovS[i,i]=0!
        Tu<-as.vector(abs(U)/(sqrt(diag(CovS)) + 1e-20) )
        k<-length(U)
        V <- matrix(0,nrow=k, ncol=k)
        for(i in 1:k){
            for(j in 1:k){
                if (abs(CovS[i,j])>1e-20)
                V[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
                else   V[i,j] <- 1e-20
            }
        }
        pTu <- as.numeric(PowerUniv(Tu,V))
    }
    
    pTu
}


PowerUniv <- function(U,V){
    n <- dim(V)[1]
    
    x <- as.numeric(max(abs(U)))
    TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))
    
    return(TER)
}

