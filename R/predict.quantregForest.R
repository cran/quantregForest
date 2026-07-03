## Vectorized row-wise quantiles (type 7, as in stats::quantile),
## sorting all rows in one pass instead of calling quantile() per row.
rowQuantilesType7 <- function(v, probs){
    m <- nrow(v)
    ntree <- ncol(v)
    ord <- order(row(v), v, na.last=TRUE)
    sorted <- matrix(v[ord], nrow=m, byrow=TRUE)
    k <- ntree - rowSums(is.na(sorted))
    out <- matrix(NA_real_, nrow=m, ncol=length(probs))
    valid <- k >= 1L
    rows <- which(valid)
    for (j in seq_along(probs)){
        h <- (k[valid] - 1) * probs[j] + 1
        lo <- floor(h)
        hi <- ceiling(h)
        g <- h - lo
        qlo <- sorted[cbind(rows, lo)]
        qhi <- sorted[cbind(rows, hi)]
        out[valid, j] <- (1 - g) * qlo + g * qhi
    }
    out
}

"predict.quantregForest" <- function(object, newdata=NULL, what=c(0.1,0.5,0.9), ...)
{
    class(object) <- "randomForest"
    if(is.null(newdata)){
        if(is.null(object[["valuesOOB"]])) stop("need to fit with option keep.inbag=TRUE \n if trying to get out-of-bag observations")
        valuesPredict <- object[["valuesOOB"]]
    }else{
        predictNodes <- attr(predict(object,newdata=newdata,nodes=TRUE),"nodes")
        rownames(predictNodes) <- NULL
        valuesPredict <- 0*predictNodes
        ntree <- ncol(object[["valuesNodes"]])
        for (tree in 1:ntree){
            valuesPredict[,tree] <- object[["valuesNodes"]][ predictNodes[,tree],tree]
        }
    }
    if(is.function(what)){
        if(is.function(what(1:4))){
            result <- apply(valuesPredict,1,what)
        }else{
            if(length(what(1:4))==1){
                result <- apply(valuesPredict,1,what)
            }else{
                result <- t(apply(valuesPredict,1,what))
            }
        }
    }else{
        if( !is.numeric(what)) stop(" `what' needs to be either a function or a vector with quantiles")
        if( min(what)<0) stop(" if `what' specifies quantiles, the minimal values needs to be non-negative")
        if( max(what)>1) stop(" if `what' specifies quantiles, the maximal values cannot exceed 1")
        result <- rowQuantilesType7(valuesPredict, what)
        if(length(what)==1){
            result <- as.vector(result)
        }else{
            colnames(result) <- paste("quantile=",what)
        }
    }
    return(result)


}
