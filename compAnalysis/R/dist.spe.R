# Overwrite spearman distance function from TreeExp2 package with metric pearson distance
    
  dist.spe = function (expMat = NULL) {

    object_n <- ncol(expMat)
    gene_n <- nrow(expMat)

    dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


    for (i in 1:(object_n-1)) {

      for (j in (i+1):object_n) {

        dis.mat[j,i] <- sqrt(1/2*(1 - cor(expMat[,i],expMat[,j],
                                  method = "spearman")))

      }

    }

    #browser()
    colnames(dis.mat) <- colnames(expMat)
    rownames(dis.mat) <- colnames(dis.mat)
    dis.mat  + t(dis.mat)

  }
