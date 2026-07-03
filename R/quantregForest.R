"quantregForest" <-
function(x,y, nthreads = 1, keep.inbag=FALSE, ...){

  ## Some checks
  #if(! class(y) %in% c("numeric","integer") )
  #  stop(" y must be numeric ")

  if(is.null(nrow(x)) || is.null(ncol(x)))
    stop(" x contains no data ")


  if( nrow(x) != length(y) )
    stop(" predictor variables and response variable must contain the same number of samples ")

  if (any(is.na(x))) stop("NA not permitted in predictors")
  if (any(is.na(y))) stop("NA not permitted in response")



  ## Check for categorial predictors with too many categories (copied from randomForest package)
   if (is.data.frame(x)) {
        ncat <- sapply(x, function(x) if(is.factor(x) && !is.ordered(x))
                       length(levels(x)) else 1)
      } else {
        ncat <- 1
    }
    maxcat <- max(ncat)
    if (maxcat > 32)
        stop("Can not handle categorical predictors with more than 32 categories.")


  ## Note that crucial parts of the computation
  ## are only invoked by the predict method
  cl <- match.call()
  cl[[1]] <- as.name("quantregForest")

  qrf <- if(nthreads > 1){
    parallelRandomForest(x=x, y=y, nthreads = nthreads,keep.inbag=keep.inbag, ...)
  }else{
    randomForest( x=x,y=y ,keep.inbag=keep.inbag,...)
  }

  nodesX <- attr(predict(qrf,x,nodes=TRUE),"nodes")
  rownames(nodesX) <- NULL
  nnodes <- max(nodesX)
  ntree <- ncol(nodesX)
  n <- nrow(x)
  valuesNodes  <- matrix(nrow=nnodes,ncol=ntree)

  ## store one randomly chosen response value per terminal node and tree
  for (tree in 1:ntree){
      shuffled <- sample.int(n)
      nodesTree <- nodesX[shuffled,tree]
      keep <- !duplicated(nodesTree)
      valuesNodes[nodesTree[keep],tree] <- y[shuffled[keep]]
  }



  qrf[["call"]] <- cl
  qrf[["valuesNodes"]] <- valuesNodes

  if(keep.inbag){

    valuesPredict <- matrix(NA_real_,nrow=n,ncol=ntree)

    ## for each tree and out-of-bag observation, sample the response of
    ## another observation falling into the same terminal node
    for (tree in 1:ntree){

      oobIdx <- which(qrf$inbag[,tree] == 0)
      if(length(oobIdx)==0) next

      nodesTree <- nodesX[,tree]
      ord <- order(nodesTree)
      cnt <- tabulate(nodesTree)
      offset <- cumsum(cnt) - cnt

      k <- cnt[nodesTree[oobIdx]]
      pick <- rep(NA_integer_,length(oobIdx))

      ## draw uniformly among all node members, redrawing the few
      ## observations that sampled themselves
      todo <- which(k > 1L)
      while(length(todo) > 0L){
        idx <- oobIdx[todo]
        pos <- offset[nodesTree[idx]] + pmin.int(k[todo], 1L + floor(runif(length(todo)) * k[todo]))
        cand <- ord[pos]
        ok <- cand != idx
        pick[todo[ok]] <- cand[ok]
        todo <- todo[!ok]
      }
      valuesPredict[oobIdx,tree] <- y[pick]
    }

    minoob <- min( rowSums(!is.na(valuesPredict)))
    if(minoob<10) stop("need to increase number of trees for sufficiently many out-of-bag observations")
    valuesOOB <- t(apply( valuesPredict,1 , function(x) sample( x[!is.na(x)], minoob)))
    qrf[["valuesOOB"]] <- valuesOOB
  }
  class(qrf) <- c("quantregForest","randomForest")

  return(qrf)
}
