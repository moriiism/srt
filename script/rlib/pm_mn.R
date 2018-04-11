###
### pm_mn.R
###
### 2018.03.07 M.Morii
###   prox-map, multinominal
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

###
###
###

SolveByProxMapMN <- function(rho.vec, Nph, D.vec, R.mat, beta, mu, L, tol, nstep, outdir, outfile.head, nrow, ncol){
    eta = 1.2
    rho.new.vec = rho.vec

    L = FindL(rho.vec, D.vec, R.mat, beta, mu, eta, nrow, ncol)
    
    for(istep in 1 : nstep){
        ik = FindIk(rho.new.vec, D.vec, R.mat, beta, mu, L, eta, nrow, ncol)
        L.pre = L
        L = eta**ik * L
        printf("SolveByProxMap: ik = %d, L = %e, L.pre = %e\n", ik, L, L.pre)
        rho.supp.vec = GetSuppVec(rho.new.vec)
        rho.new.vec = ProxMap(rho.new.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.supp.vec)

        ## dump
        array = array(rho.new.vec, dim=c(ncol, nrow))
        file.tmp = sprintf("temp%3.3d.fits", istep)
        writeFITSim(array, file=file.tmp)
        ## dump
    }
    return (rho.new.vec)
}

FindL <- function(rho.vec, D.vec, R.mat, beta, mu, eta, nrow, ncol){
    nstep = 1000
    L.org = 1e-10
    L = L.org
    for(istep in 1 : nstep){
        L = eta**istep * L.org
        rho.supp.vec = GetSuppVec(rho.vec)
        sigma.vec = FuncSigma(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.supp.vec)
        printf("FindL: L = %e, min(sigma.vec) = %e, max(sigma.vec) = %e\n", L, min(sigma.vec), max(sigma.vec))
        rho.new.vec = mapply(max, 0, sigma.vec - min(sigma.vec))

        printf("sum(rho.new.vec) = %e\n", sum(rho.new.vec))
        
        if( sum(rho.new.vec) < 1){
            break
        }
    }
    return(L)
}


FindIk <- function(rho.vec, D.vec, R.mat, beta, mu, L, eta, nrow, ncol){
    ik.max = 1000
    ik = 0
    L.org = L
    while(ik <= ik.max){
        L = eta**ik * L.org
        rho.supp.vec = GetSuppVec(rho.vec)
        printf("FindIK: GetNZeroVec(rho.supp.vec) = %d\n", GetNZeroVec(rho.supp.vec))
        
        pLy = ProxMap(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.supp.vec)

        printf("FindIK: GetNZeroVec(pLy) = %d\n", GetNZeroVec(pLy))
        
        rho.new.supp.vec = GetSuppVec(pLy)
        qminusf = QMinusF(pLy, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.new.supp.vec, rho.supp.vec)
        printf("FindIk: (ik, qminusf) = (%d, %e)\n", ik, qminusf)
        if(qminusf >= 0){
            break
        }
        ik = ik + 1
    }
    return(ik)
}


ProxMap <- function(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.supp.vec){
    
    sigma.vec = FuncSigma(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.supp.vec)
    printf("ProxMap: min(sigma.vec) = %e, max(sigma.vec) = %e\n", min(sigma.vec), max(sigma.vec))

    rho.new.vec = mapply(max, 0, sigma.vec - min(sigma.vec))
    printf("sum = %e\n", sum(rho.new.vec))
    if( sum(rho.new.vec) < 1 ) {
        num.nonzero = nrow * ncol - GetNZeroVec(rho.supp.vec)
        rho.new.vec = sigma.vec + (1 - sum(sigma.vec)) / num.nonzero
    }
    else {
        index.vec = order(sigma.vec)
        for( index in 1 : ( length(sigma.vec) - 1) ){
            ## printf("index = %d\n", index)
            
            rho.new.vec = mapply(max, 0, sigma.vec - sigma.vec[index.vec[index]])
            rho.new.p1.vec = mapply(max, 0, sigma.vec - sigma.vec[index.vec[index + 1]])

            ## printf("sum(rho.new.vec) = %e, sum(rho.new.p1.vec) = %e\n", sum(rho.new.vec), sum(rho.new.p1.vec) )
            
            if( sum(rho.new.vec) > 1 && sum(rho.new.p1.vec) < 1){
                ## sum = sum( sigma.vec[index.vec[(index + 1) : length(sigma.vec) ]] )
                ##num.nonzero = length(sigma.vec) - index
                ##rho.new.vec = sigma.vec + (1 - sum ) / num.nonzero

                
                break
            }
        }
    }
    printf("sum(rho.new.vec) = %e\n", sum(rho.new.vec))
    return(rho.new.vec)
}


KLDiv <- function(y.vec, y.new.vec, R.mat)
{
    q.vec = R.mat %*% y.vec
    q.new.vec = R.mat %*% y.new.vec
    q.vec = q.vec / sum(q.vec)
    q.new.vec = q.new.vec / sum(q.new.vec)

    ### avoid zero component
    ans = 0.0
    for( index in 1 : length(q.new.vec) ){
        ## printf("index = %d, q.new.vec[index] = %e\n", index, q.new.vec[index])
        
        if(q.new.vec[index] > 0.0){
            ans = ans + q.new.vec[index] * log( q.new.vec[index] / q.vec[index] )
        }
    }
    return (ans)
}

###
###
###

FuncLp <- function(rho.vec, sigma.vec, tau, rho.supp.vec){
    term1 = sum( (rho.vec - sigma.vec) * (rho.vec - sigma.vec) * rho.supp.vec )
    term2 = tau * (sum(rho.vec * rho.supp.vec) - 1)
    ans = term1 + term2
    return (ans)
}

QMinusF <- function(rho.new.vec, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.new.supp.vec, rho.supp.vec)
{
    rho.supp2.vec = rho.new.supp.vec * rho.supp.vec
    term1 =      FuncF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.supp2.vec)
    term2 = -1 * FuncF(rho.new.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.supp2.vec)
    term3 = sum( (rho.new.vec - rho.vec) * DiffF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.supp2.vec) )
    term4 = L / 2.0 * sum( (rho.new.vec - rho.vec) * (rho.new.vec - rho.vec) * rho.supp2.vec)
    ans = term1 + term2 + term3 + term4
    return (ans)
}

FuncSigma <- function(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.supp.vec){
    term1 = rho.vec * rho.supp.vec
    term2 = - 1.0 / L * DiffF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.supp.vec)
    ans.vec = term1 + term2

    printf("FuncSigma: min(term1) = %e, max(term1) = %e \n", min(term1), max(term1))
    printf("FuncSigma: min(term2) = %e, max(term2) = %e \n", min(term2), max(term2))
    
    return(ans.vec)
}

FuncF <- function(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.supp.vec)
{
    num.vec = R.mat %*% (rho.vec * rho.supp.vec)
    term1 = -1 * sum( D.vec * log( num.vec ) )
    term2 = (1.0 - beta) * SumLogPlus(rho.vec)
    term3 = mu * TermV(rho.vec, nrow, ncol)
    ans = term1 + term2 + term3
    return(ans)
}

DiffF <- function(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.supp.vec)
{
    num.vec = R.mat %*% (rho.vec * rho.supp.vec)
    term1 = -1 * (t(R.mat) %*% DivVect(D.vec, num.vec))
    term2 = (1 - beta) * DivVect(rho.supp.vec, rho.vec)
    term3 = mu * DiffTermV(rho.vec, nrow, ncol)
    ans = term1 + term2 + term3
    return (ans)
}

DiffTermV <- function(rho.vec, nrow, ncol)
{
    mat = matrix(rho.vec, nrow=nrow, ncol=ncol)
    mat.aug.p.tmp = rbind(mat, mat[nrow,])
    mat.aug.p = cbind(mat.aug.p.tmp, mat.aug.p.tmp[,ncol])
    mat.aug.m.tmp = rbind(mat[1,], mat)
    mat.aug.m = cbind(mat.aug.m.tmp[,1], mat.aug.m.tmp)
    
    mat.p.0 = mat.aug.p[2:(nrow+1),1:ncol]
    mat.0.p = mat.aug.p[1:nrow,2:(ncol+1)]
    mat.m.0 = mat.aug.m[1:nrow,2:(ncol+1)]
    mat.0.m = mat.aug.m[2:(nrow+1),1:ncol]

    mat.diff = 2.0 * ( 4 * mat - mat.m.0 - mat.p.0 - mat.0.m - mat.0.p )
    rho.diff.vec = as.vector(mat.diff)
    rho.diff.vec = rho.diff.vec / (nrow * ncol)
    return (rho.diff.vec)
}

TermV <- function(rho.vec, nrow, ncol)
{
    mat = matrix(rho.vec, nrow=nrow, ncol=ncol)
    mat.aug.tmp = rbind(mat, mat[nrow,])
    mat.aug = cbind(mat.aug.tmp, mat.aug.tmp[,ncol])
    mat.p1.p0 = mat.aug[2:(nrow+1),1:ncol]
    mat.p0.p1 = mat.aug[1:nrow,2:(ncol+1)]
    ans = sum( (mat - mat.p1.p0) * (mat - mat.p1.p0) ) + sum( (mat - mat.p0.p1) * (mat - mat.p0.p1) )
    ans = ans / (nrow * ncol)
    return (ans)
}


###
###
###

SumLogPlus <- function(vec){
    ans = 0.0
    for(index in 1 : length(vec)){
        if(vec[index] > 0.0){
            ans = ans + log(vec[index])
        }
    }
    return(ans)
}

DivVect <- function(num.vec, den.vec){
    epsilon = 1.0e-15
    len = length(num.vec)
    ans.vec = rep(0.0, len)
    for(index in 1 : len){
        if(abs(num.vec[index]) > epsilon){
            ans.vec[index] = num.vec[index] / den.vec[index]
        }
        else {
            ans.vec[index] = 0.0
        }
    }
    return(ans.vec)
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

GetNZeroVec <- function(vec){
    epsilon = 1.0e-15
    nzero = 0
    for(index in 1 : length(vec)){
        if(abs(vec[index]) < epsilon){
            nzero = nzero + 1
        }
    }
    return(nzero)
}

