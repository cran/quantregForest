## Fit a random forest on several cores by growing ntree/nthreads trees
## per fork and combining the forests. Forking is unavailable on Windows,
## where a single forest is grown instead.
"parallelRandomForest" <-
  function(x, y, nthreads = 1, keep.inbag = FALSE, ntree = 500, ...) {

    if (nthreads <= 1 || get_os() == "windows")
      return(randomForest(x = x, y = y, keep.inbag = keep.inbag,
                          ntree = ntree, ...))

    ntreeEach <- ceiling(ntree / nthreads)
    rflist <- parallel::mclapply(seq_len(nthreads),
                                 function(i) randomForest(x = x, y = y,
                                                          keep.inbag = keep.inbag,
                                                          ntree = ntreeEach, ...),
                                 mc.cores = nthreads)
    failed <- !vapply(rflist, inherits, logical(1), what = "randomForest")
    if (any(failed))
      stop("parallel random forest fit failed: ",
           paste(unique(vapply(rflist[failed], function(r)
             conditionMessage(attr(r, "condition")), character(1))), collapse = "; "))

    do.call(combine, rflist)
  }
