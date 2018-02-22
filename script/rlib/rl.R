###
### rl.R
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


SolveByRL <- function(x.vec, D.vec, R.mat, nrow, ncol){
    eta = 1.2
    y.vec = x.vec
    x.pre.vec = x.vec
    t = 1
    cost = 0.0
    cost.pre = 0.0
    tol = 1e-10
    k.max = 100
   

    for(k in 1 : k.max){
        printf("SolveByProxMap: k = %d\n", k)

        array = array(y.vec, dim=c(ncol, nrow))

        ProxMap(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)

        kldiv = KLDiv(y.vec, y.new.vec, R.mat)
        printf("SolveByProxMap: KL divergence = %e\n", kldiv)

        if(k > 1 && kldiv < tol){
            break
        }
        ##if(k > 1 && abs(cost.pre - cost) / cost < tol){
        ##    break
        ##}
        t = t.new
        y.vec = y.new.vec
        x.pre.vec = x.vec
        cost.pre = cost
    }
    x.vec = y.vec
    return (x.vec)
}


# y.vec ---> y.new.vec
ProxMap <- function(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
{
##    printf("ProxMap: max(y.vec), min(y.vec) = %e, %e\n", max(y.vec), min(y.vec))
    
    y.new.vec = y.vec
    y.pre.vec = y.vec
    nem.step = 500
    tol = 1.0e-8
    for(iem.step in 1:nem.step){
        sigma.vec = y.new.vec - 1.0 / L * DiffF(y.new.vec, mu, nrow, ncol, lin.or.log)
        num.vec = R.mat %*% y.new.vec
        term1 = (t(R.mat) %*% (D.vec / num.vec)) * y.new.vec
        term2 = (1.0 - beta) * nrow * ncol * y.new.vec / sum(y.new.vec) - (1.0 - beta)
        m.vec = term1 + term2
        if(min(num.vec) < 1.0e-10){
            break
        }
        y.new.vec = mapply(Mfunc, m.vec, sigma.vec, L)
        kldiv = KLDiv(y.pre.vec, y.new.vec, R.mat)
##        printf("KL divergence = %e\n", kldiv)
##        printf("iem.step = %d\n", iem.step)
        if(kldiv < tol){
            printf("ProxMap: KL divergence = %e\n", kldiv)
            printf("ProxMap: iem.step = %d\n", iem.step)
            break
        }
        y.pre.vec = y.new.vec
    }
    return(y.new.vec)
}


FilterEpsilon <- function(xval){
    ans = xval
    epsilon = 1.0e-10
    if(xval < epsilon){
        ans = epsilon
    }
    return(ans)
}


KLDiv <- function(y.vec, y.new.vec, R.mat)
{
    q.vec = R.mat %*% y.vec
    q.new.vec = R.mat %*% y.new.vec

    q.vec = q.vec / sum(q.vec)
    q.new.vec = q.new.vec / sum(q.new.vec)

### debug
##    printf("KLDiv: max(y.vec), max(y.new.vec) = %e, %e\n", max(y.vec), max(y.new.vec))
##    printf("KLDiv: max(q.vec), max(q.new.vec) = %e, %e\n", max(q.vec), max(q.new.vec))
##    printf("KLDiv: min(y.vec), min(y.new.vec) = %e, %e\n", min(y.vec), min(y.new.vec))
##    printf("KLDiv: min(q.vec), min(q.new.vec) = %e, %e\n", min(q.vec), min(q.new.vec))

###    q.filt.vec = mapply(FilterEpsilon, q.vec)
###    q.new.filt.vec = mapply(FilterEpsilon, q.new.vec)


##    printf("KLDiv: min(q.filt.vec), min(q.new.filt.vec) = %e, %e\n", min(q.filt.vec), min(q.new.filt.vec))

##    printf("KLDiv: min: log( q.new.filt.vec / q.filt.vec ) = %e\n", min(log( q.new.filt.vec / q.filt.vec )))
##    printf("KLDiv: max: log( q.new.filt.vec / q.filt.vec ) = %e\n", max(log( q.new.filt.vec / q.filt.vec )))
    
###     ans = sum( q.new.filt.vec * log( q.new.filt.vec / q.filt.vec ) )
    ans = sum( q.new.vec * log( q.new.vec / q.vec ) )

    return (ans)
}

