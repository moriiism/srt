###
### pm_mn_e.R
###
### 2018.03.07 M.Morii
###   prox-map, multinominal, with eplison offset
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

LineSearch <- function(x.vec, x.new.vec, R.mat, D.vec, beta, mu, nrow, ncol){
    x1     = x.vec[1]
    x1.new = x.new.vec[1]
    x2.vec = x.vec[2:length(x.vec)]
    x2.new.vec = x.new.vec[2:length(x.new.vec)]
    theta.vec = log( x2.vec / (1 - x1))
    theta.new.vec = log( x2.new.vec / (1 - x1.new))
    nstep = 100
    logl.init = FuncF(x.vec, D.vec, R.mat, beta, mu, nrow, ncol)
    logl.pre  = logl.init
    x.pre.vec = x.vec

    printf("logl.init = %e\n", logl.init)
    

    outfile = sprintf("temp.dat")
    fopen(outfile)
    fprintf(outfile, "skip sing\n")
    fprintf(outfile, "read\n")
    
    eta = 1.2
    for( istep in 1 : nstep){
        factor = eta**istep
        lx1.this = factor * (log(x1.new) - log(x1)) + log(x1)
        x1.this = exp(lx1.this)
        theta.this.vec = factor * (theta.new.vec - theta.vec) + theta.vec
        alpha = sum(exp(theta.this.vec))
        x2.this.vec = (1 - x1.this) * exp(theta.this.vec) / alpha
        x.this.vec = c(x1.this, x2.this.vec)

        logl = FuncF(x.this.vec, D.vec, R.mat, beta, mu, nrow, ncol)

        fprintf(outfile, "%d  %e\n", istep, logl - logl.init)

        printf("logl = %e\n", logl)

        if(logl.pre < logl){
            printf("factor(istep) = %e (%d)\n", factor, istep)
            if(istep != 1){
                x.new.vec = x.pre.vec
            }
            break
        }
        logl.pre = logl
        x.pre.vec = x.this.vec
    }
    return(x.new.vec)
}

SolveByProxMapMN <- function(rho.vec, Nph, D.vec, R.mat, beta, mu, L, tol, nstep, outdir, outfile.head, nrow, ncol, epsilon){
    outfile = sprintf("%s/%s_moni.dat", outdir, outfile.head)
    fopen(outfile)
    fprintf(outfile, "# istep, kldiv, logl, logl.inc, delta.logl, tdiff\n")

    outfile.time.logl = sprintf("%s/%s_time_logl.dat", outdir, outfile.head)
    fopen(outfile.time.logl)
    fprintf(outfile.time.logl, "# tdiff logl.inc\n")

    time.st = Sys.time()
    logl.init = FuncF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol)
    logl.pre = logl.init
    fprintf(outfile, "0  0  %.10e  0.0  0.0  0.0\n", logl.init)
    
    eta = 1.2
    rho.new.vec = rho.vec
    rho.pre.vec = rho.vec
    for(istep in 1 : nstep){
        ik = FindIk(rho.new.vec, D.vec, R.mat, beta, mu, L, eta, nrow, ncol)
        L.pre = L
        L = eta**ik * L
        printf("SolveByProxMap: ik = %d, L = %e, L.pre = %e\n", ik, L, L.pre)
        rho.tmp.vec = ProxMap(rho.new.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon)
        rho.new.vec = LineSearch(rho.new.vec, rho.tmp.vec, R.mat, D.vec, beta, mu, nrow, ncol)

        kldiv = KLDiv(rho.pre.vec, rho.new.vec, R.mat)
        logl = FuncF(rho.new.vec, D.vec, R.mat, beta, mu, nrow, ncol)
       
        delta.logl = logl - logl.pre
        logl.inc   = logl - logl.init
        time = Sys.time()
        tdiff = difftime(time, time.st, units=c("secs"))
        printf("%d  %e  %.10e  %e  %e  %e\n",
               istep, kldiv, logl, logl.inc, delta.logl, tdiff)
        fprintf(outfile, "%d  %e  %.10e  %e  %e  %e\n",
                istep, kldiv, logl, logl.inc, delta.logl, tdiff)
        fprintf(outfile.time.logl, "%e  %e\n", tdiff, logl.inc)

        ## save
##        if(istep %% 100 == 0){
            outimg = sprintf("%s/outimg_%d.fits", outdir, istep)
            array = array(rho.new.vec * Nph, dim=c(ncol, nrow))
            writeFITSim(array, file=outimg)
##        }
        
        if(kldiv < tol){
            break
        }

        ## for next
        rho.pre.vec = rho.new.vec
        logl.pre = logl

    }
    return (rho.new.vec)
}


FindIk <- function(rho.vec, D.vec, R.mat, beta, mu, L, eta, nrow, ncol){
    ik.max = 1000
    ik = 0
    L.org = L

    logl.init = FuncF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol)
    while(ik <= ik.max){
        L = eta**ik * L.org
        pLy = ProxMap(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon)
        qminusf = QMinusF(pLy, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)
        printf("FindIk: (ik, L, qminusf) = (%d, %e, %e)\n", ik, L, qminusf)
        logl = FuncF(pLy, D.vec, R.mat, beta, mu, nrow, ncol)
        printf("logl, logl.init = %e, %e \n", logl, logl.init)
        
        if(logl > logl.init){
            printf("logl(%e) > logl.init(%e), then next \n", logl, logl.init)
            ik = ik + 1
            next
        }
        if(qminusf >= 0){
            outimg = sprintf("findik_%d.fits", ik)
            array = array(pLy, dim=c(ncol, nrow))
            writeFITSim(array, file=outimg)
            break
        }
        ik = ik + 1
    }
    return(ik)
}

FuncSum <- function(tau, sigma.vec, epsilon){
    rho.vec = mapply(max, epsilon, sigma.vec - tau / 2.0)
    ans = sum(rho.vec) - 1.0
    return(ans)
}

ProxMap <- function(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon){
    npix = nrow * ncol
    sigma.vec = FuncSigma(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)
    tau.min = 2 * ( min(sigma.vec) - epsilon )
    rho.new.vec = rho.vec
    if( FuncSum(tau.min, sigma.vec, epsilon) <= 0){
        rho.new.vec = sigma.vec + (1 - sum(sigma.vec)) / (nrow * ncol)
    }
    else {
        index.vec = order(sigma.vec)
        tau.ord.vec = 2 * ( sigma.vec[index.vec] - epsilon )
        ipoint = npix
        for( istep in 1 : 10000 ){
            delta.index = as.integer(npix / 2**istep)
            if(delta.index < 1){
                break
            }
            if(FuncSum(tau.ord.vec[ipoint], sigma.vec, epsilon) <= 0){
                ipoint = ipoint - delta.index
            }
            else {
                ipoint = ipoint + delta.index
            }
        }
        if(FuncSum(tau.ord.vec[ipoint], sigma.vec, epsilon) <= 0){
            for( index in ipoint : 1){
                if( FuncSum(tau.ord.vec[index], sigma.vec, epsilon) > 0 &&  FuncSum(tau.ord.vec[index + 1], sigma.vec, epsilon) <= 0 ){
                    rho.new.vec[index.vec[1 : index]] = epsilon
                    term = (sum(sigma.vec[ index.vec[(index + 1) : (nrow * ncol)]]) + index * epsilon - 1) / ( (nrow * ncol) - index )
                    rho.new.vec[index.vec[(index + 1) : (nrow * ncol)]] = sigma.vec[ index.vec[(index + 1) : (nrow * ncol)]] - term
                    break
                }
            }
        }
        else {
            for( index in ipoint : npix){
                if( FuncSum(tau.ord.vec[index], sigma.vec, epsilon) > 0 &&  FuncSum(tau.ord.vec[index + 1], sigma.vec, epsilon) <= 0 ){
                    rho.new.vec[index.vec[1 : index]] = epsilon
                    term = (sum(sigma.vec[ index.vec[(index + 1) : (nrow * ncol)]]) + index * epsilon - 1) / ( (nrow * ncol) - index )
                    rho.new.vec[index.vec[(index + 1) : (nrow * ncol)]] = sigma.vec[ index.vec[(index + 1) : (nrow * ncol)]] - term
                    break
                }
            }
        }
    }
    return(rho.new.vec)
}

KLDiv <- function(y.vec, y.new.vec, R.mat)
{
    q.vec = R.mat %*% y.vec
    q.new.vec = R.mat %*% y.new.vec
    q.vec = q.vec / sum(q.vec)
    q.new.vec = q.new.vec / sum(q.new.vec)
    ans = 0.0
    for( index in 1 : length(q.new.vec) ){
        if(q.new.vec[index] > 0.0){
            ans = ans + q.new.vec[index] * log( q.new.vec[index] / q.vec[index] )
        }
    }
    return (ans)
}

###
###
###

FuncLp <- function(rho.vec, sigma.vec, tau){
    term1 = sum( (rho.vec - sigma.vec) * (rho.vec - sigma.vec) )
    term2 = tau * (sum(rho.vec) - 1)
    ans = term1 + term2
    return (ans)
}

QMinusF <- function(rho.new.vec, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)
{
    term1 =      FuncF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol)
    term2 = -1 * FuncF(rho.new.vec, D.vec, R.mat, beta, mu, nrow, ncol)
    term3 = sum( (rho.new.vec - rho.vec) * DiffF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol) )
    term4 = L / 2.0 * sum( (rho.new.vec - rho.vec) * (rho.new.vec - rho.vec) )
    ans = term1 + term2 + term3 + term4
    return (ans)
}

FuncSigma <- function(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol){
    term1 = rho.vec
    term2 = -1.0 / L * DiffF(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol)
    ans.vec = term1 + term2
    return(ans.vec)
}

FuncF <- function(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol)
{
    num.vec = R.mat %*% rho.vec
    term1 = -1 * sum( D.vec * log( num.vec ) )
    term2 = (1.0 - beta) * sum(log(rho.vec))
    term3 = mu * TermV(rho.vec, nrow, ncol)
    ans = term1 + term2 + term3
    return(ans)
}

DiffF <- function(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol)
{
    num.vec = R.mat %*% rho.vec
    term1 = -1 * (t(R.mat) %*% (D.vec / num.vec))
    term2 = (1 - beta) / rho.vec
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


######
######
######
###
###SumLogPlus <- function(vec, epsilon){
###    ans = 0.0
###    for(index in 1 : length(vec)){
###        if(vec[index] > epsilon){
###            ans = ans + log(vec[index])
###        }
###    }
###    return(ans)
###}
###
###DivVect <- function(num.vec, den.vec, epsilon){
###    len = length(num.vec)
###    ans.vec = rep(0.0, len)
###    for(index in 1 : len){
###        if(abs(num.vec[index]) > epsilon){
###            ans.vec[index] = num.vec[index] / den.vec[index]
###        }
###        else {
###            ans.vec[index] = 0.0
###        }
###    }
###    return(ans.vec)
###}
###
###GetSuppVec <- function(vec, epsilon){
###    supp.vec = rep(0, length(vec))
###    for(index in 1 : length(vec)){
###        if(vec[index] > epsilon){
###            supp.vec[index] = 1
###        }
###    }
###    return(supp.vec)
###}
###
###GetNZeroVec <- function(vec, epsilon){
###    nzero = 0
###    for(index in 1 : length(vec)){
###        if(abs(vec[index]) < epsilon){
###            nzero = nzero + 1
###        }
###    }
###    return(nzero)
###}
###
