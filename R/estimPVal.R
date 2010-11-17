estimPVal <- function(object,x,permute.n=10,per.covariate=FALSE,parallel=FALSE,multicore=FALSE,trace=FALSE,...) {
    if (is.matrix(object$penalty)) {
        stop("Uncertainty cannot be estimated with penalty updates")
    }
    
    if (object$standardize) x <- scale(x,center=object$meanx,scale=object$sdx)

    permute.index <- matrix(NA,permute.n,nrow(x))    
    for (actual.permute in 1:permute.n) permute.index[actual.permute,] <- sample(nrow(x))

    time <- object$time
    status <- object$status
    unpen.index <- object$unpen.index
    if (length(object$unpen.index) > 0) {
        pen.index <- (1:ncol(x))[-object$unpen.index]
    } else {
        pen.index <- 1:ncol(x)
    }
    stepno <- object$stepno
    penalty <- rep(object$penalty,2)

    eval.permute <- function(actual.permute,...) {
        if (trace) cat("permutation",actual.permute,"\n")

        actual.x <- cbind(x,x[permute.index[actual.permute,],pen.index])
        
        permute.res <- CoxBoost(time=time,status=status,x=actual.x,unpen.index=unpen.index,
                                standardize=FALSE,stepno=stepno,penalty=penalty,trace=FALSE,...)
    
        actual.score <- colMeans(permute.res$scoremat[,1:length(pen.index)])
        null.score <- colMeans(permute.res$scoremat[,(length(pen.index)+1):ncol(permute.res$scoremat)])

        if (per.covariate) {
            return(actual.score <= null.score)
        } else {
            return(unlist(lapply(actual.score,function(arg) mean(arg <= null.score))))
        }
    }

    done.parallel <- FALSE

    if (parallel) {
        if (!require(snowfall)) {
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            sfLibrary(CoxBoost)
            sfExport("time","status","x","unpen.index","pen.index","stepno","penalty","trace","permute.index")
            permute.mat <- matrix(unlist(sfClusterApplyLB(1:permute.n,eval.permute,...)),length(pen.index),permute.n)
            done.parallel <- TRUE            
        }
    } 

    if (!done.parallel & multicore) {
        if (!require(multicore)) {
            warning("package 'multicore' not found, i.e., parallelization cannot be performed using this package")
        } else {
            if (multicore > 1) {
                permute.mat <- matrix(unlist(mclapply(1:permute.n,eval.permute,mc.preschedule=FALSE,mc.cores=multicore,...)),length(pen.index),permute.n)
            } else {
                permute.mat <- matrix(unlist(mclapply(1:permute.n,eval.permute,mc.preschedule=FALSE,...)),length(pen.index),permute.n)
            }
            done.parallel <- TRUE
        }        
    }

    if (!done.parallel) {
        permute.mat <- matrix(unlist(lapply(1:permute.n,eval.permute)),length(pen.index),permute.n)
    }

    rowMeans(permute.mat)
}
