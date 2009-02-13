
efron.weightmat <- function(time,status) {
    n <- length(time)
    uncens <- which(status == 1)
    weightmat <- matrix(0,n,length(uncens))

    rept <- rep(0,n)
    for (i in 1:n) rept[i] <- sum(time[i:n]==time[i] & status[i:n]==1)

    for (i in 1:length(uncens)) {
        weightmat[time >= time[uncens][i] | (status != 0 & status != 1), i] <- 1
        tie <- time == time[uncens[i]] & status==1
        di <- max(rept[tie])
        weightmat[tie, i] <- weightmat[tie, i] - (di - rept[uncens[i]])/di
    }

    #   check for competing risks scenario
    if (any(status != 0 & status != 1)) {
        cens.ind <- ifelse(status == 0,1,0)
        invcensprob <- rep(NA,length(time))
        invcensprob[order(time)] <- summary(survival::survfit(Surv(time,cens.ind) ~ 1),times=sort(time))$surv
                
        for (i in 1:length(uncens)) {
            current.invcens <- invcensprob[uncens][i]
            weightmat[,i] <- weightmat[,i] * 
                current.invcens/ifelse(time < time[uncens][i],invcensprob,current.invcens)
        }
    }
    
    weightmat
}

CoxBoost <- function(time,status,x,unpen.index=NULL,standardize=TRUE,stepno=100,penalty=9*sum(status==1),
                     stepsize.factor=1,sf.scheme=c("sigmoid","linear"),pendistmat=NULL,connected.index=NULL,
                     trace=FALSE)
{
    sf.scheme <- match.arg(sf.scheme)
    
    object <- list()
    object$stepno <- stepno
    object$unpen.index <- unpen.index
    pen.index <- 1:ncol(x)
    if (!is.null(unpen.index)) pen.index <- pen.index[-unpen.index]
    if (length(penalty) < length(pen.index)) penalty <- rep(penalty[1],length(pen.index))
    object$penalty <- NULL
    clusterno <- 1

    if (is.null(colnames(x))) {
        object$xnames <- paste("V",1:ncol(x),sep="")
    } else {
        object$xnames <- colnames(x)
    }

    if (!is.null(connected.index) && any(connected.index %in% unpen.index)) {
        stop("sets of unpenalized and connected covariates may not overlap")
    }

    if (!is.null(pendistmat) && is.null(connected.index)) {
        if (ncol(pendistmat) == ncol(x) - length(unpen.index)) {
            if (!is.null(unpen.index)) {
                connected.index <- (1:ncol(x))[-unpen.index]
            } else {
                connected.index <- 1:ncol(x)
            }
        } else {
            stop("'connected.index' is missing and cannot be guessed")
        }
    }

    if (!is.null(unpen.index)) {
        if (!is.null(connected.index)) connected.index <- match(connected.index,(1:ncol(x))[-unpen.index])
        unpen.x <- x[,unpen.index,drop=FALSE]
        x <- x[,-unpen.index,drop=FALSE]
    }

    uncens <- which(status == 1)
    n <- nrow(x)
    n.uncens <- length(uncens)
    p <- ncol(x)

    penpos <- match(1:p,connected.index)

    object$n <- n
    object$p <- p
    object$time <- time
    object$status <- status
    object$event.times <- sort(unique(time[uncens]))

    object$coefficients <- matrix(NA,stepno+1,p)
    if (!is.null(unpen.index)) {
        unpen.coefficients <- matrix(NA,stepno+1,ncol(unpen.x)) 
    } else {
        unpen.coefficients <- NULL
    }
    object$linear.predictor <- matrix(NA,stepno+1,n)

    object$meanx <- rep(0,length(object$xnames))
    object$sdx <- rep(1,length(object$xnames))

    if (standardize) {
        x <- scale(x)
        object$meanx[pen.index] <- attr(x,"scaled:center")
        object$sdx[pen.index] <- attr(x,"scaled:scale")
    }
    
    object$Lambda <- matrix(NA,stepno+1,length(object$event.times))

    #   create clusters

    if (clusterno > 1) {
        object$cluster.assignment <- kmeans(t(x),centers=clusterno,iter.max=50)$cluster
    } else {
        object$cluster.assignment <- NULL
    }
    
    #   Efron handling of ties

    weightmat <- efron.weightmat(time,status)

    actual.beta <- rep(0,p)
    if (!is.null(unpen.index)) actual.unpen.beta <- rep(0,ncol(unpen.x))
    actual.linear.predictor <- rep(0,n)
    actual.risk.score <- rep(1,n)

    ml.fraction <- rep(0,p)

    #   boosting iterations

    x.double.vec <- as.double(x)
    weight.double.vec <- as.double(weightmat)
    uncens.C <- as.integer(uncens - 1)
    

    for (actual.step in 0:stepno) {
        if (actual.step > 0) object$penalty <- rbind(object$penalty,penalty)
        weightmat.times.risk <- weightmat*actual.risk.score
        weightmat.times.risk.sum <- colSums(weightmat.times.risk)

        #   update unpenalized covariates by one estimation step
        
        if (!is.null(unpen.index)) {
            if (actual.step == 1) {
                #   calculations from step zero do not have to be repeated
                unpen.coefficients[actual.step+1,] <- actual.unpen.beta
            } else {
                x.bar <- (t(weightmat.times.risk) %*% unpen.x) / weightmat.times.risk.sum
                U <- colSums(unpen.x[uncens,] - x.bar)

                I <- matrix(0,ncol(unpen.x),ncol(unpen.x))
                for (i in 1:n.uncens) {
                    x.minus.bar <- t(t(unpen.x) - x.bar[i,])        
                    I <- I + (t(x.minus.bar*(weightmat.times.risk[,i])) %*% x.minus.bar)/weightmat.times.risk.sum[i]
                }

                unpen.beta.delta <- drop(solve(I) %*% U)
                actual.unpen.beta <- actual.unpen.beta + unpen.beta.delta
                unpen.coefficients[actual.step+1,] <- actual.unpen.beta

                actual.linear.predictor <- actual.linear.predictor + drop(unpen.x %*% unpen.beta.delta)
                actual.risk.score <- exp(drop(actual.linear.predictor))
                weightmat.times.risk <- weightmat*actual.risk.score
                weightmat.times.risk.sum <- colSums(weightmat.times.risk)
            }
        }

        if (actual.step == 0) {
            actual.Lambda <- rep(NA,length(object$event.times))
            for (i in seq(along=object$event.times)) {
                actual.mask <- time[uncens] <= object$event.times[i]
                actual.Lambda[i] <- sum(1/weightmat.times.risk.sum[actual.mask])
            }
        
            object$coefficients[actual.step+1,] <- actual.beta
            object$linear.predictor[actual.step+1,] <- actual.linear.predictor
            object$Lambda[actual.step+1,] <- actual.Lambda
            
            next
        }           
                
        if (clusterno > 1) {
            actual.cluster.mask <- object$cluster.assignment == ifelse(actual.step %% clusterno == 0,clusterno,actual.step %% clusterno)
            
            res <- .C(find_best,
                      as.double(x[,actual.cluster.mask]),
                      as.integer(nrow(x)),
                      as.integer(sum(actual.cluster.mask)),
                      uncens.C,
                      as.integer(length(uncens)),
                      as.double(actual.beta[actual.cluster.mask]),
                      as.double(actual.risk.score),
                      as.double(actual.linear.predictor),
                      weight.double.vec,
                      as.double(weightmat.times.risk),
                      as.double(weightmat.times.risk.sum),
                      as.double(penalty[actual.cluster.mask]),
                      min.index=integer(1),
                      min.deviance=double(1),
                      min.beta.delta=double(1),
                      DUP=FALSE
                      )

            #cat("C result:",res$min.index,"(",res$min.deviance,")\n")

            min.index <- (1:length(actual.cluster.mask))[actual.cluster.mask][res$min.index]
            min.deviance <- res$min.deviance
            min.beta.delta <- res$min.beta.delta            
        } else {
            #   C implementation

            res <- .C(find_best,
                      x.double.vec,
                      as.integer(nrow(x)),
                      as.integer(ncol(x)),
                      uncens.C,
                      as.integer(length(uncens)),
                      as.double(actual.beta),
                      as.double(actual.risk.score),
                      as.double(actual.linear.predictor),
                      weight.double.vec,
                      as.double(weightmat.times.risk),
                      as.double(weightmat.times.risk.sum),
                      as.double(penalty),
                      min.index=integer(1),
                      min.deviance=double(1),
                      min.beta.delta=double(1),
                      DUP=FALSE
                      )

            #cat("C result:",res$min.index,"(",res$min.deviance,")\n")

            min.index <- res$min.index
            min.deviance <- res$min.deviance
            min.beta.delta <- res$min.beta.delta

            # #   R implementation using a loop
            # 
            # min.index <- NULL
            # min.deviance <- NULL
            # min.beta <- NULL
            # 
            # for (i in 1:ncol(x)) {
            #     actual.x.bar <- apply(weightmat*actual.risk.score*x[,i],2,sum)/apply(weightmat*actual.risk.score,2,sum)
            #     U <- x[uncens,i] - actual.x.bar
            #     V <- apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[,i],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum)
            # 
            #     candidate.beta.delta <- sum(U)/(sum(V)+penalty)    
            #     candidate.beta <- actual.beta
            #     candidate.beta[i] <- candidate.beta[i] + candidate.beta.delta
            #     candidate.loglik <- sum(drop(x[uncens,,drop=FALSE] %*% candidate.beta) - log(apply(weightmat*exp(drop(x %*% candidate.beta)),2,sum)))
            #     candidate.deviance <- -2*candidate.loglik
            # 
            #     if (i == 1 || candidate.deviance < min.deviance) {
            #         min.index <- i
            #         min.deviance <- candidate.deviance
            #         min.beta <- candidate.beta
            #     }
            # }            
        }


        #cat("selected:",min.index,"(",min.deviance,")\n")
        if (trace) cat(object$xnames[pen.index][min.index]," ",sep="")

        #   update the maximum likelihood fractions needed for penalty distribution
        if (!is.null(pendistmat)) { 
            actual.x.bar <- apply(weightmat*actual.risk.score*x[,min.index],2,sum)/apply(weightmat*actual.risk.score,2,sum)
            I <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[,min.index],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))
            nu <- I / (I + penalty[min.index])
            
            ml.fraction[min.index] <- ml.fraction[min.index] + (1-ml.fraction[min.index])*nu
        }

        actual.beta[min.index] <- actual.beta[min.index] + min.beta.delta
        #print(actual.beta)
        actual.linear.predictor <- drop(x %*% actual.beta)
        if (!is.null(unpen.index)) actual.linear.predictor <- actual.linear.predictor + drop(unpen.x %*% actual.unpen.beta)

        actual.risk.score <- exp(drop(actual.linear.predictor))
        weightmat.times.risk <- weightmat*actual.risk.score
        weightmat.times.risk.sum <- colSums(weightmat.times.risk)

        actual.Lambda <- rep(NA,length(object$event.times))
        for (i in seq(along=object$event.times)) {
            actual.mask <- time[uncens] <= object$event.times[i]
            actual.Lambda[i] <- sum(1/weightmat.times.risk.sum[actual.mask])
        }
        
        object$coefficients[actual.step+1,] <- actual.beta
        object$linear.predictor[actual.step+1,] <- actual.linear.predictor
        object$Lambda[actual.step+1,] <- actual.Lambda

        #   update the penalty if the user has chosen any other value than the default
        actual.stepsize.factor <- ifelse(length(stepsize.factor) >= min.index,stepsize.factor[min.index],stepsize.factor[1])
        if (actual.stepsize.factor != 1 && ml.fraction[min.index] < 1) { 
            if (is.null(pendistmat) || (min.index %in% connected.index && any(pendistmat[penpos[min.index],] != 0))) {
                actual.x.bar <- apply(weightmat*actual.risk.score*x[,min.index],2,sum)/apply(weightmat*actual.risk.score,2,sum)
                I <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[,min.index],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))
                old.penalty <- penalty[min.index]
                if (sf.scheme == "sigmoid") {
                    new.nu <- max(1 - (1-(I/(I+penalty[min.index])))^actual.stepsize.factor,0.00001)  # prevent penalty -> Inf
                    penalty[min.index] <- (1/new.nu - 1)*I
                } else {
                    penalty[min.index] <- (1/actual.stepsize.factor - 1)*I + penalty[min.index]/actual.stepsize.factor
                }
                if (penalty[min.index] < 0) penalty[min.index] <- 0
                
                if (!is.null(pendistmat)) {
                    if (trace) {
                        cat("\npenalty update for ",object$xnames[pen.index][min.index]," (mlf: ",round(ml.fraction[min.index],3),"): ",old.penalty," -> ",penalty[min.index],"\n",sep="")
                        cat("connected:\n")
                    }

                    connected <- connected.index[which(pendistmat[penpos[min.index],] != 0)]
                    
                    for (actual.target in connected) {
                        if (ml.fraction[actual.target] < 1) {
                            if (trace) cat(object$xnames[pen.index][actual.target]," (mlf: ",round(ml.fraction[actual.target],3),"): ",penalty[actual.target]," -> ",sep="")
                            actual.x.bar <- apply(weightmat*actual.risk.score*x[,actual.target],2,sum)/apply(weightmat*actual.risk.score,2,sum)
                            I.target <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[,actual.target],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))                            
                            new.target.penalty <- pendistmat[penpos[min.index],penpos[actual.target]]*(1 - ml.fraction[actual.target])*I.target/
                                                      ((1-actual.stepsize.factor)*pendistmat[penpos[min.index],penpos[actual.target]]*(1-ml.fraction[min.index])*I/(I+old.penalty) +
                                                       (1-ml.fraction[actual.target])*I.target/(I.target+penalty[actual.target])) - 
                                                      I.target
                            if (new.target.penalty > 0) penalty[actual.target] <- new.target.penalty
                            if (trace) cat(penalty[actual.target],"\n")
                        }
                    }
                }
            }
        }
    }
    if (trace) cat("\n")
    
    #   combine penalized and unpenalized covariates
    if (!is.null(object$unpen.index)) {
        object$p <- object$p + length(object$unpen.index)
        combined.coefficients <- matrix(NA,nrow(object$coefficients),object$p)
        combined.coefficients[,pen.index] <- object$coefficients
        combined.coefficients[,object$unpen.index] <- unpen.coefficients
        object$coefficients <- combined.coefficients
    }
    
    class(object) <- "CoxBoost"
    object$logplik <- predict(object,type="logplik")
    
    object
}

print.CoxBoost <- function(x,...) {
    cat(x$stepno,"boosting steps resulting in",
        sum(x$coefficients[x$stepno,] != 0),
        "non-zero coefficients",ifelse(is.null(x$unpen.index),"",paste("(with",length(x$unpen.index),"being mandatory)")),
        "\n")
    cat("partial log-likelihood:",x$logplik,"\n")
}

summary.CoxBoost <- function(object,...) {
    print(object)
    cat("\n")
    if (!is.null(object$unpen.index)) {
        cat("Parameter estimates for mandatory covariates at boosting step ",object$stepno,":\n",sep="")
        print(matrix(signif(object$coefficients[object$stepno+1,object$unpen.index],4),length(object$unpen.index),1,dimnames=list(object$xnames[object$unpen.index],c("Estimate"))))
        cat("\n")
    }

    cat("Optional covariates with non-zero coefficients at boosting step ",object$stepno,":\n",sep="")
    cat("parameter estimate > 0:\n",paste(object$xnames[object$coefficients[object$stepno+1,] > 0],collapse=", "),"\n")
    cat("parameter estimate < 0:\n",paste(object$xnames[object$coefficients[object$stepno+1,] < 0],collapse=", "),"\n")
}

predict.CoxBoost <- function(object,newdata=NULL,newtime=NULL,newstatus=NULL,at.step=NULL,times=NULL,type=c("lp","logplik","risk","CIF"),...) {
    if (is.null(at.step)) at.step <- object$stepno

    if (is.null(newdata)) {
        linear.predictor <- object$linear.predictor[at.step+1,,drop=FALSE]
    } else {
        linear.predictor <- t(scale(newdata,center=object$meanx,scale=object$sdx) %*% t(object$coefficients[at.step+1,,drop=FALSE]))
    }
    
    type <- match.arg(type)
    if (type == "lp") return(linear.predictor)
    
    if (type == "logplik") {
        if (is.null(newdata)) {
            newtime <- object$time
            newstatus <- object$status
        } else {
            if (is.null(newtime) || is.null(newstatus)) stop("'newtime' and 'newstatus' required for prediction on new data")
        }
        uncens <- which(newstatus == 1)
        weightmat <- efron.weightmat(newtime,newstatus)
        
        logplik <- c()
        
        for (i in seq(along=at.step)) {
            logplik <- c(logplik,sum(linear.predictor[i,uncens] - log(apply(weightmat*exp(linear.predictor[i,]),2,sum))))
        }
        return(logplik)
    }
    
    if (type == "risk" || type == "CIF") {
        if (length(at.step) > 1) warning("predicted risk is only calculated for a single step (the first in 'at.step')")
        if (is.null(times)) times <- unique(object$time)

        breslow.Lambda <- unlist(lapply(times,function(x) ifelse(x < object$event.times[1],0,object$Lambda[at.step[1]+1,rev(which(object$event.times <= x))[1]])))
        pred.risk <- exp(exp(linear.predictor[1,]) %*% -t(breslow.Lambda))
        if (type == "risk") {
            return(pred.risk)
        } else {
            return(1-pred.risk)
        }
    }
    
    NULL
}

cv.CoxBoost <- function(time,status,x,maxstepno=100,K=10,type=c("verweij","naive"),parallel=FALSE,upload.x=TRUE,folds=NULL,
                        trace=FALSE,...) {
    type <- match.arg(type)

    if (!is.null(folds) && length(folds) != K) stop("'folds' has to be of length 'K'")

    if (is.null(folds)) {
        if (K >= sum(status != 0)) {    #   leave-one-out cross-validation
            if (type == "verweij") {
                folds <- as.list(1:length(time))
            } else {
                folds <- as.list(which(status != 0))            
            }
        } else {
            while(TRUE) {
                folds <- split(sample(1:nrow(x)), rep(1:K, length = nrow(x)))
        
                #   make sure there is at least one event in training and test folds respectively
                #   Note: the Verweij approach actually could deal with folds that contain only
                #   censored observations, but this is expected to considerably increase variability
                #   and therefore alse prevented
                if (!any(unlist(lapply(folds,function(fold) sum(status[fold]))) == 0) &&
                    !any(unlist(lapply(folds,function(fold) sum(status[-fold]))) == 0)) 
                {
                    break
                }
            }
        }
    }
    
    criterion <- NULL
    
    eval.fold <- function(actual.fold,...) {
        if (trace) cat("cv fold ",actual.fold,": ",sep="")
        cv.fit <- CoxBoost(time=time[-folds[[actual.fold]]],status=status[-folds[[actual.fold]]],x=x[-folds[[actual.fold]],,drop=FALSE],
                           stepno=maxstepno,trace=trace,...)

        if (type == "verweij") {
            full.ploglik <- predict(cv.fit,newdata=x,newtime=time,newstatus=status,type="logplik",at.step=0:maxstepno)
            fold.ploglik <- predict(cv.fit,newdata=x[-folds[[actual.fold]],,drop=FALSE],newtime=time[-folds[[actual.fold]]],newstatus=status[-folds[[actual.fold]]],
                                    type="logplik",at.step=0:maxstepno)
                                    
            return(full.ploglik - fold.ploglik)
        } else {
            return(predict(cv.fit,newdata=x[folds[[actual.fold]],,drop=FALSE],newtime=time[folds[[actual.fold]]],
                                                 newstatus=status[folds[[actual.fold]]],
                                                 type="logplik",at.step=0:maxstepno))
        }
    }

    if (parallel) {
        if (!require(snowfall)) {
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
            criterion <- matrix(unlist(lapply(1:length(folds),eval.fold,...)),nrow=length(folds),byrow=TRUE)        
        } else {
            sfLibrary(CoxBoost)
            if (upload.x) {
                sfExport("time","status","x","maxstepno","trace","type","folds")
            } else {
                sfExport("time","status","maxstepno","trace","type","folds")
            }
            criterion <- matrix(unlist(sfClusterApplyLB(1:length(folds),eval.fold,...)),nrow=length(folds),byrow=TRUE)                
        }
    } else {
        criterion <- matrix(unlist(lapply(1:length(folds),eval.fold,...)),nrow=length(folds),byrow=TRUE)        
    }
    mean.criterion <- apply(criterion,2,mean)
    
    list(mean.logplik=mean.criterion,se.logplik=apply(criterion,2,sd)/sqrt(nrow(criterion)),optimal.step=which.max(mean.criterion)-1,
         folds=folds)
}


optimCoxBoostPenalty <- function(time,status,x,minstepno=50,maxstepno=200,start.penalty=9*sum(status==1),
                                 iter.max=10,upper.margin=0.05,parallel=FALSE,trace=FALSE,...)
{
    if (parallel) {
        if (!require(snowfall)) {
            parallel <- FALSE
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            sfExport("x")            
        }
    }
    
    actual.penalty <- start.penalty
    
    #   default: start from a large penalty and go down, when gone to far use small steps up
    step.up <- 1.2
    step.down <- 0.5
    
    actual.res <- NULL
    
    for (i in 1:iter.max) {
        if (trace) cat("iteration",i,": evaluating penalty",actual.penalty,"\n")
        
        actual.res <- cv.CoxBoost(time=time,status=status,x=x,maxstepno=maxstepno,penalty=actual.penalty,
                                  parallel=parallel,upload.x=FALSE,trace=trace,...)
        actual.max <- actual.res$optimal.step
        
        if (trace) cat("maximum partial log-likelihood at boosting step",actual.max,"\n")
        if (actual.max >= minstepno && actual.max < maxstepno*(1-upper.margin)) break

        #   check whether we are in a scenario where penalty is far to low to start with
        if (i == 1 && actual.max < minstepno) {
            step.up <- 2
            step.down <- 0.8
        }

        if (actual.max < minstepno) {
            actual.penalty <- actual.penalty * step.up
        } else {
            actual.penalty <- actual.penalty * step.down
        }
        
        if (i == iter.max) warning("Exceeded iter.max in search for penalty parameter")
    }

    list(penalty=actual.penalty,cv.res=actual.res)
}

