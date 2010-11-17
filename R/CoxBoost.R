
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

CoxBoost <- function(time,status,x,unpen.index=NULL,standardize=TRUE,subset=1:length(time),
                     stepno=100,penalty=9*sum(status[subset]==1),
                     stepsize.factor=1,sf.scheme=c("sigmoid","linear"),pendistmat=NULL,connected.index=NULL,
                     x.is.01=FALSE,return.score=TRUE,trace=FALSE)
{
    sf.scheme <- match.arg(sf.scheme)
    
    object <- list()

    #   reduce response to subset
    time <- time[subset]
    status <- status[subset]
    
    #   reorder observations according to time to speed up computations
    object$time <- time
    object$status <- status    
    time.order <- order(time,decreasing=TRUE)
    subset.time.order <- (1:nrow(x))[subset][time.order]
    status <- status[time.order]
    time <- time[time.order]
    
    object$stepno <- stepno
    object$unpen.index <- unpen.index
    pen.index <- 1:ncol(x)
    if (!is.null(unpen.index)) pen.index <- pen.index[-unpen.index]
    if (length(penalty) < length(pen.index)) penalty <- rep(penalty[1],length(pen.index))
    
    if (any(stepsize.factor != 1)) {
        object$penalty <- matrix(NA,stepno,length(penalty))     
    } else {
        object$penalty <- penalty
    }

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
        unpen.x <- x[subset.time.order,unpen.index,drop=FALSE]
    }

    uncens <- which(status == 1)
    n <- length(status)
    n.uncens <- length(uncens)
    p <- length(pen.index)

    penpos <- match(1:p,connected.index)

    object$n <- n
    object$p <- p
    object$event.times <- sort(unique(time[uncens]))

    object$coefficients <- Matrix(0,stepno+1,p)
    if (!is.null(unpen.index)) {
        unpen.coefficients <- matrix(NA,stepno+1,ncol(unpen.x)) 
    } else {
        unpen.coefficients <- NULL
    }
    object$linear.predictor <- matrix(NA,stepno+1,n)

    object$meanx <- rep(0,length(object$xnames))
    object$sdx <- rep(1,length(object$xnames))
    if (standardize) {
        pen.sdx <- apply(x[subset,pen.index],2,sd)
        pen.sdx <- ifelse(pen.sdx == 0,1,pen.sdx)
        pen.meanx <- apply(x[subset,pen.index],2,mean)        
        x[subset,pen.index] <- scale(x[subset,pen.index],center=pen.meanx,scale=pen.sdx)

        object$meanx[pen.index] <- pen.meanx
        object$sdx[pen.index] <- pen.sdx        
        object$standardize <- TRUE
    } else {
        object$standardize <- FALSE
    }
    
    object$Lambda <- matrix(NA,stepno+1,length(object$event.times))
    if (return.score) object$scoremat <- matrix(NA,max(1,stepno),object$p)

    #   Efron handling of ties

    weightmat <- efron.weightmat(time,status)

    actual.beta <- rep(0,p)
    if (!is.null(unpen.index)) actual.unpen.beta <- rep(0,ncol(unpen.x))
    actual.linear.predictor <- rep(0,n)
    actual.risk.score <- rep(1,n)

    ml.fraction <- rep(0,p)

    #   boosting iterations

    x.double.vec <- as.double(x[subset.time.order,pen.index])
    weight.double.vec <- as.double(weightmat)
    max.nz.vec <- as.integer(apply(weightmat,2,function(arg) max(which(arg != 0))))
    max.1.vec <- as.integer(c(0,rev(cummin(rev(apply(weightmat,2,function(arg) ifelse(!any(arg != 1),length(arg),min(which(arg != 1)-1))))))))
    uncens.C <- as.integer(uncens - 1)
    
    warnstep <- NULL

    for (actual.step in 0:stepno) {
        if (actual.step > 0 && any(stepsize.factor != 1)) {
            object$penalty[stepno,] <- penalty            
        }         
        
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
                
        if (x.is.01) {
            res <- .C("find_best01",
                      x.double.vec,
                      as.integer(n),
                      as.integer(p),
                      uncens.C,
                      as.integer(length(uncens)),
                      as.double(actual.beta),
                      as.double(actual.risk.score),
                      as.double(actual.linear.predictor),
                      weight.double.vec,
                      max.nz.vec,
                      max.1.vec,
                      as.double(weightmat.times.risk),
                      as.double(weightmat.times.risk.sum),
                      as.double(penalty),
                      warncount=integer(1),
                      min.index=integer(1),
                      min.deviance=double(1),
                      min.beta.delta=double(1),
                      score.vec=double(p),
                      DUP=FALSE
                      )                
        } else {
            res <- .C("find_best",
                      x.double.vec,
                      as.integer(n),
                      as.integer(p),
                      uncens.C,
                      as.integer(length(uncens)),
                      as.double(actual.beta),
                      as.double(actual.risk.score),
                      as.double(actual.linear.predictor),
                      weight.double.vec,
                      max.nz.vec,
                      max.1.vec,
                      as.double(weightmat.times.risk),
                      as.double(weightmat.times.risk.sum),
                      as.double(penalty),
                      warncount=integer(1),
                      min.index=integer(1),
                      min.deviance=double(1),
                      min.beta.delta=double(1),
                      score.vec=double(p),
                      DUP=FALSE
                      )                
        }

        if (is.null(warnstep) && res$warncount > 0) warnstep <- actual.step

        min.index <- res$min.index
        min.deviance <- res$min.deviance
        min.beta.delta <- res$min.beta.delta

        if (return.score) object$scoremat[actual.step,] <- res$score.vec

        #cat("selected:",min.index,"(",min.deviance,")\n")
        if (trace) cat(object$xnames[pen.index][min.index]," ",sep="")

        #   update the maximum likelihood fractions needed for penalty distribution
        if (!is.null(pendistmat)) { 
            actual.x.bar <- apply(weightmat*actual.risk.score*x[subset.time.order,pen.index[min.index]],2,sum)/apply(weightmat*actual.risk.score,2,sum)
            I <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[subset.time.order,pen.index[min.index]],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))
            nu <- I / (I + penalty[min.index])
            
            ml.fraction[min.index] <- ml.fraction[min.index] + (1-ml.fraction[min.index])*nu
        }

        actual.beta[min.index] <- actual.beta[min.index] + min.beta.delta
        #print(actual.beta)
        actual.linear.predictor <- actual.linear.predictor + x[subset.time.order,pen.index[min.index]]*min.beta.delta

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
                I.index <- min.index

                actual.x.bar <- apply(weightmat*actual.risk.score*x[subset.time.order,pen.index[min.index]],2,sum)/apply(weightmat*actual.risk.score,2,sum)
                I <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[subset.time.order,pen.index[min.index]],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))

                if (!is.null(pendistmat)) {
                    connected <- connected.index[which(pendistmat[penpos[min.index],] != 0)]
                    if (length(connected) > 0) {
                        change.index <- connected[ml.fraction[connected] < 1]
                        I.index <- c(I.index,change.index)
                    }
                }

                I.vec <- .C("get_I_vec",
                          as.double(x[subset.time.order,pen.index[I.index]]),
                          as.integer(n),
                          as.integer(length(I.index)),
                          as.integer(length(uncens)),
                          as.double(weightmat.times.risk),
                          as.double(weightmat.times.risk.sum),
                          I.vec=double(length(I.index)),
                          DUP=FALSE
                          )$I.vec

                old.penalty <- penalty[min.index]
                if (sf.scheme == "sigmoid") {
                    new.nu <- max(1 - (1-(I.vec[1]/(I.vec[1]+penalty[min.index])))^actual.stepsize.factor,0.00001)  # prevent penalty -> Inf
                    penalty[min.index] <- (1/new.nu - 1)*I
                } else {
                    penalty[min.index] <- (1/actual.stepsize.factor - 1)*I.vec[1] + penalty[min.index]/actual.stepsize.factor
                }
                if (penalty[min.index] < 0) penalty[min.index] <- 0
                
                if (length(I.vec) > 1) {
                    if (trace) {
                        cat("\npenalty update for ",object$xnames[pen.index][min.index]," (mlf: ",round(ml.fraction[min.index],3),"): ",old.penalty," -> ",penalty[min.index],"\n",sep="")
                    }

                    change.I <- I.vec[2:length(I.vec)]
                    if (trace) cat("adjusting penalty for covariates",paste(object$xnames[pen.index][change.index],collapse=", ",sep=""),"\n")

                    new.target.penalty <- pendistmat[penpos[min.index],penpos[change.index]]*
                                          (1 - ml.fraction[change.index])*change.I/
                                          ((1-actual.stepsize.factor)*pendistmat[penpos[min.index],penpos[change.index]]*
                                                                      (1-ml.fraction[min.index])*I.vec[1]/(I.vec[1]+old.penalty) +
                                           (1-ml.fraction[change.index])*change.I/(change.I+penalty[change.index])) -
                                          change.I
                                          
                    penalty[change.index] <- ifelse(new.target.penalty > 0,new.target.penalty,penalty[change.index])

                    #   loop version
                    # for (actual.target in connected) {
                    #     if (ml.fraction[actual.target] < 1) {
                    #         if (trace) cat(object$xnames[pen.index][actual.target]," (mlf: ",round(ml.fraction[actual.target],3),"): ",penalty[actual.target]," -> ",sep="")
                    #         actual.x.bar <- apply(weightmat*actual.risk.score*x[subset.time.order,pen.index[actual.target]],2,sum)/apply(weightmat*actual.risk.score,2,sum)
                    #         I.target <- sum(apply((weightmat*actual.risk.score)*t(t(matrix(rep(x[subset.time.order,pen.index[actual.target]],n.uncens),nrow(weightmat),ncol(weightmat))) - actual.x.bar)^2,2,sum)/apply(weightmat*actual.risk.score,2,sum))                            
                    #         new.target.penalty <- pendistmat[penpos[min.index],penpos[actual.target]]*(1 - ml.fraction[actual.target])*I.target/
                    #                                   ((1-actual.stepsize.factor)*pendistmat[penpos[min.index],penpos[actual.target]]*(1-ml.fraction[min.index])*I/(I+old.penalty) +
                    #                                    (1-ml.fraction[actual.target])*I.target/(I.target+penalty[actual.target])) - 
                    #                                   I.target
                    #         if (new.target.penalty > 0) penalty[actual.target] <- new.target.penalty
                    #         if (trace) cat(penalty[actual.target],"\n")
                    #     }
                    # }                    
                }
            }
        }
    }
    if (trace) cat("\n")

    if (!is.null(warnstep)) warning(paste("potentially attempted to move towards a nonexisting maximum likelihood solution around step",warnstep))
    
    #   combine penalized and unpenalized covariates
    if (!is.null(object$unpen.index)) {
        object$p <- object$p + length(object$unpen.index)
        combined.coefficients <- Matrix(0,nrow(object$coefficients),object$p)
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
        sum(x$coefficients[x$stepno+1,] != 0),
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

predict.CoxBoost <- function(object,newdata=NULL,newtime=NULL,newstatus=NULL,subset=NULL,at.step=NULL,times=NULL,type=c("lp","logplik","risk","CIF"),...) {
    if (is.null(at.step)) at.step <- object$stepno

    if (is.null(newdata)) {
        linear.predictor <- object$linear.predictor[at.step+1,,drop=FALSE]
    } else {
        if (is.null(subset)) {
            subset.index <- 1:nrow(newdata)
        } else {
            subset.index <- (1:nrow(newdata))[subset]
            if (!is.null(newtime)) newtime <- newtime[subset]
            if (!is.null(newstatus)) newstatus <- newstatus[subset]
        }

        nz.index <- which(Matrix::colSums(abs(object$coefficients[at.step+1,,drop=FALSE])) > 0)
        if (length(nz.index) > 0) {
            if (object$standardize) {
                linear.predictor <- as.matrix(Matrix::tcrossprod(object$coefficients[at.step+1,nz.index,drop=FALSE],scale(newdata[subset.index,nz.index,drop=FALSE],center=object$meanx[nz.index],scale=object$sdx[nz.index])))
            } else {
                # linear.predictor <- t(newdata[subset.index,] %*% t(object$coefficients[at.step+1,,drop=FALSE]))
                linear.predictor <- as.matrix(Matrix::tcrossprod(object$coefficients[at.step+1,nz.index,drop=FALSE],newdata[subset.index,nz.index,drop=FALSE]))
            }
        } else {
            linear.predcitor <- matrix(0,length(at.step),length(subset.index))
        }
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

cv.CoxBoost <- function(time,status,x,maxstepno=100,K=10,type=c("verweij","naive"),parallel=FALSE,upload.x=TRUE,
                        multicore=FALSE,folds=NULL,
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
        cv.fit <- CoxBoost(time=time,status=status,x=x,subset=-folds[[actual.fold]],
                           stepno=maxstepno,return.score=FALSE,trace=trace,...)

        if (type == "verweij") {
            full.ploglik <- predict(cv.fit,newdata=x,newtime=time,newstatus=status,type="logplik",at.step=0:maxstepno)
            fold.ploglik <- predict(cv.fit,newdata=x,newtime=time,newstatus=status,subset=-folds[[actual.fold]],
                                    type="logplik",at.step=0:maxstepno)
                                    
            return(full.ploglik - fold.ploglik)
        } else {
            return(predict(cv.fit,newdata=x,newtime=time,newstatus=status,subset=folds[[actual.fold]],
                                                 type="logplik",at.step=0:maxstepno))
        }
    }

    eval.success <- FALSE
    
    if (parallel) {
        if (!require(snowfall)) {
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed using this package")
        } else {
            sfLibrary(CoxBoost)
            if (upload.x) {
                sfExport("time","status","x","maxstepno","trace","type","folds")
            } else {
                sfExport("time","status","maxstepno","trace","type","folds")
            }
            criterion <- matrix(unlist(sfClusterApplyLB(1:length(folds),eval.fold,...)),nrow=length(folds),byrow=TRUE)                
            eval.success <- TRUE
        }
    } 

    if (!eval.success & multicore) {
        if (!require(multicore)) {
            warning("package 'multicore' not found, i.e., parallelization cannot be performed using this package")
        } else {
            if (multicore > 1) {
                criterion <- matrix(unlist(mclapply(1:length(folds),eval.fold,mc.preschedule=FALSE,mc.cores=multicore,...)),nrow=length(folds),byrow=TRUE)
            } else {
                criterion <- matrix(unlist(mclapply(1:length(folds),eval.fold,mc.preschedule=FALSE,...)),nrow=length(folds),byrow=TRUE)                
            }
            eval.success <- TRUE
        }        
    }

    if (!eval.success) {
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

