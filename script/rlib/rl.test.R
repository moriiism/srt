###
### rl.test.R
###
### 2018.01.27 M.Morii
###   simple R-L method
### 

GetFitsHeader <- function(file)
{
    printf("--- GetFitsHeader --- \n")
    library(FITSio)

    fits <- readFITS(in.file)
    bitpix = fits$hdr[which(fits$hdr=="BITPIX")+1]
    naxis  = fits$hdr[which(fits$hdr=="NAXIS")+1]
    naxis1 = fits$hdr[which(fits$hdr=="NAXIS1")+1]
    naxis2 = fits$hdr[which(fits$hdr=="NAXIS2")+1]
    
    print(bitpix)
}

LoadModel <- function(respdir)
{
    printf("--- LoadModel --- \n")
    library(FITSio)

    in.file = sprintf("%s/gimage_000_000.img", respdir)
    fits <- readFITS(in.file)
    bitpix = as.integer(fits$hdr[which(fits$hdr=="BITPIX")+1])
    naxis  = as.integer(fits$hdr[which(fits$hdr=="NAXIS")+1])
    naxis1 = as.integer(fits$hdr[which(fits$hdr=="NAXIS1")+1])
    naxis2 = as.integer(fits$hdr[which(fits$hdr=="NAXIS2")+1])

    printf("bitpix = %d\n", bitpix)
    printf("naxis  = %d\n", naxis)
    printf("naxis1  = %d\n", naxis1)
    printf("naxis2  = %d\n", naxis2)

    # row: detector
    # col: sky
    
    nrow = naxis1 * naxis2
    nibin = 60
    njbin = 60
    ncol = nibin * njbin
    printf("nrow = %d, ncol = %d\n", nrow, ncol)
    mat = matrix(0.0, nrow=nrow, ncol=ncol)
    # fits$imDat
    for(jbin in 0:(njbin - 1)){
        for(ibin in 0:(nibin -1) ){
            in.file = sprintf("%s/gimage_%3.3d_%3.3d.img", respdir, ibin, jbin)
            ## print(in.file)
            fits <- readFITS(in.file)
            kbin = nibin * jbin + ibin + 1
            mat[, kbin] = fits$imDat / sum(fits$imDat)
        }
    }
    return (mat)
}

LoadData <- function(file)
{
    printf("--- LoadData --- \n")
    library(FITSio)

    fits <- readFITS(file)
    bitpix = as.integer(fits$hdr[which(fits$hdr=="BITPIX")+1])
    naxis  = as.integer(fits$hdr[which(fits$hdr=="NAXIS")+1])
    naxis1 = as.integer(fits$hdr[which(fits$hdr=="NAXIS1")+1])
    naxis2 = as.integer(fits$hdr[which(fits$hdr=="NAXIS2")+1])

    printf("bitpix = %d\n", bitpix)
    printf("naxis  = %d\n", naxis)
    printf("naxis1  = %d\n", naxis1)
    printf("naxis2  = %d\n", naxis2)
    D.vec = as.vector(fits$imDat)
    return (D.vec)
}

SolveByRL <- function(x.vec, D.vec, R.mat, outfile, tol){
    k.max = 10000
    for(k in 1 : k.max){
        den.vec = R.mat %*% x.vec
        x.new.vec = t(R.mat) %*% (D.vec / den.vec) * x.vec

        lx.vec = log(x.vec)
        lx.new.vec = log(x.new.vec)
        delta.lx.vec = lx.new.vec - lx.vec
        nstep = 100

        logl.pre = FuncL(x.vec, R.mat, D.vec)
        x.pre.vec = x.vec
        for( istep in 0 : nstep){
            lx.this.vec = lx.vec + delta.lx.vec * (1.0 + 0.2 * istep)
            x.this.vec = exp(lx.this.vec)
            logl = FuncL(x.this.vec, R.mat, D.vec)
            if(logl.pre > logl){
                x.new.vec = x.pre.vec
                break
            }
            logl.pre = logl
            x.pre.vec = x.this.vec
        }

        kldiv = KLDiv(x.vec, x.new.vec, R.mat)
        logl = FuncL(x.new.vec, R.mat, D.vec)
        printf("k = %d, KLdiv = %e, min(den.vec) = %e, max(den.vec) = %e\n",
               k, kldiv, min(den.vec), max(den.vec))
        fprintf(outfile, "%d  %e  %e\n", k, kldiv, logl)
        if(kldiv < tol){
            break
        }

        ## for next
        x.vec = x.new.vec
    }
    return (x.vec)
}

KLDiv <- function(x.vec, x.new.vec, R.mat)
{
    q.vec = R.mat %*% x.vec
    q.new.vec = R.mat %*% x.new.vec
    q.vec = q.vec / sum(q.vec)
    q.new.vec = q.new.vec / sum(q.new.vec)

    ### avoid zero component
    ans = 0.0
    for( index in 1 : length(q.new.vec) ){
        if(q.new.vec[index] > 0.0){
            ans = ans + q.new.vec[index] * log( q.new.vec[index] / q.vec[index] )
        }
    }
    return (ans)
}

FuncL <- function(y.vec, R.mat, D.vec)
{
    num.vec = R.mat %*% y.vec
    ans = sum( D.vec * log( num.vec ) ) - sum(y.vec)
    return(ans)
}

FuncLSupp <- function(y.vec, R.mat, D.vec, supp.vec)
{
    y.supp.vec = y.vec * supp.vec
    num.vec = R.mat %*% y.supp.vec
    ans = sum( D.vec * log( num.vec ) ) - sum(y.supp.vec)
    return(ans)
}

GetSuppVec <- function(vec){
    supp.vec = rep(0, length(vec))
    for(index in 1 : length(vec)){
        if(vec[index] > 0.0){
            supp.vec[index] = 1
        }
    }
    return(supp.vec)
}
