estimPVal <- function(object,x,permute.n=10,parallel=FALSE,trace=FALSE) {
    if (nrow(object$penalty) > 1 && any(object$penalty[1,] != object$penalty[2,])) {
        stop("Uncertainty cannot be estimated with penalty updates")
    }
    
    x <- scale(x,center=object$meanx,scale=object$sdx)

    permute.index <- matrix(NA,permute.n,nrow(x))    
    for (actual.permute in 1:permute.n) permute.index[actual.permute,] <- sample(nrow(x))

    time <- object$time
    status <- object$status
    unpen.index <- object$unpen.index
    stepno <- object$stepno
    penalty <- object$penalty[1,]

    eval.permute <- function(actual.permute,...) {
        if (trace) cat("permutation",actual.permute,"\n")

        actual.x <- x[permute.index[actual.permute],]
        if (length(object$unpen.index) > 0) {
            actual.x[,unpen.index] <- x[,unpen.index]
        } 
        
        permute.res <- CoxBoost(time=time,status=status,x=x,unpen.index=unpen.index,
                                standardize=FALSE,stepno=stepno,penalty=penalty,trace=FALSE)
                               
        apply(apply(-abs(permute.res$scoremat),1,rank),1,median)
    }

    done.parallel <- FALSE

    if (parallel) {
        if (!require(snowfall)) {
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            sfLibrary(CoxBoost)
            sfExport("time","status","x","unpen.index","stepno","penalty","trace")
            permute.rank.vec <- unlist(sfClusterApplyLB(1:permute.n,eval.permute))
            done.parallel <- TRUE            
        }
    } 
    
    if (!done.parallel) {
        permute.rank.vec <- unlist(lapply(1:permute.n,eval.permute))
    }

    full.rank.median <- apply(apply(-abs(object$scoremat),1,rank),1,median)
    
    unlist(lapply(seq(along=full.rank.median),function(i) mean(permute.rank.vec <= full.rank.median[i])))
}
