###
### pm_mn_case2.R
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

SolveByProxMapMN <- function(rho.vec, Nph, D.vec, R.mat, beta, mu, L, tol, nstep, outdir, outfile.head, nrow, ncol, epsilon){
    outfile = sprintf("%s/%s_moni.dat", outdir, outfile.head)
    fopen(outfile)
    fprintf(outfile, "# istep, kldiv, logl, logl.inc, delta.logl, tdiff\n")

    outfile.time.logl = sprintf("%s/%s_time_logl.dat", outdir, outfile.head)
    fopen(outfile.time.logl)
    fprintf(outfile.time.logl, "# tdiff logl.inc\n")

    time.st = Sys.time()
    logl.init = FuncL(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol)
    logl.pre = logl.init
    fprintf(outfile, "0  0  %.10e  0.0  0.0  0.0\n", logl.init)
    
    eta = 1.2
    rho.new.vec = rho.vec
    rho.pre.vec = rho.vec
    for(istep in 1 : nstep){
        
        ik = FindIk(rho.new.vec, D.vec, R.mat, beta, mu, L, eta, nrow, ncol, epsilon)
        L.pre = L
        L = eta**ik * L
        printf("SolveByProxMap: ik = %d, L = %e, L.pre = %e\n", ik, L, L.pre)
        ans.proxmap = ProxMap(rho.new.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon)
        rho.new.vec = ans.proxmap[[1]]
        
        kldiv = KLDiv(rho.pre.vec, rho.new.vec, R.mat)

        ## rho.pre.supp.vec = GetSuppVec(rho.pre.vec, epsilon)
        ## rho.new.supp.vec = GetSuppVec(rho.new.vec, epsilon)
        ## logl.pre = FuncL.supp(rho.pre.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.new.supp.vec)
        ## logl     = FuncL.supp(rho.new.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.new.supp.vec)

        logl.pre = FuncL(rho.pre.vec, D.vec, R.mat, beta, mu, nrow, ncol)
        logl     = FuncL(rho.new.vec, D.vec, R.mat, beta, mu, nrow, ncol)
        
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
            outimg = sprintf("%s/outimg_%4.4d.fits", outdir, istep)
            array = array(rho.new.vec * Nph, dim=c(ncol, nrow))
            writeFITSim(array, file=outimg)
##        }

        if( -1.0e-3 < delta.logl && delta.logl < 0.0){
            break
        }
        if(kldiv < tol){
            break
        }

        ## for next
        rho.pre.vec = rho.new.vec
        logl.pre = logl

    }
    return (rho.new.vec)
}


FindIk <- function(rho.vec, D.vec, R.mat, beta, mu, L, eta, nrow, ncol, epsilon){
    ik.max = 1000
    ik = 0
    L.org = L
    while(ik <= ik.max){
        L = eta**ik * L.org
        ans.proxmap = ProxMap(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon)
        pLy = ans.proxmap[[1]]
        ans.proxmap.flag.good = ans.proxmap[[2]]
        
        ## qminusf = QMinusF.supp(pLy, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon)
        qminusf = QMinusF(pLy, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)
        printf("FindIk: (ik, L, qminusf) = (%d, %e, %e)\n", ik, L, qminusf)
        if(qminusf >= 0 && ans.proxmap.flag.good == 1){
            break
        }
        ik = ik + 1
    }
    return(ik)
}

ProxMap <- function(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon){
    npix = nrow * ncol
    sigma.vec = FuncSigma(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)

    ## printf("proxmap: sigma.vec.min = %e, sigma.vec.max = %e\n", min(sigma.vec), max(sigma.vec))
    
    rho.new.vec = rho.vec
    nem = 1000
    flag.good = 1
    tau.pre = 1e-10
    for(iem in 1 : nem){
        lem = FuncLem(rho.new.vec, D.vec, R.mat, sigma.vec, L)
        m.vec = FuncM(rho.new.vec, D.vec, R.mat)

        ## printf("em: m.vec.min = %e, m.vec.max = %e\n", min(m.vec), max(m.vec))

        tau.thres.vec = FuncTauThres(sigma.vec, m.vec, L, epsilon)
        tau.thres.min = min(tau.thres.vec)
        tau.thres.max = max(tau.thres.vec)
        ##        printf("tau.thres.min = %e, tau.thres.max = %e\n", tau.thres.min, tau.thres.max)
        
        nnewton = 100
        tol.newton = 1.0e-3
        if(tau.pre > tau.thres.min){
            tau = tau.pre
        }
        else {
            tau = tau.thres.max
        }
        for(inewton in 1 : nnewton){
            tau = tau - FuncS(tau, sigma.vec, m.vec, L, epsilon) / FuncDiffS(tau, sigma.vec, m.vec, L, epsilon)
            if( abs(FuncS(tau, sigma.vec, m.vec, L, epsilon) ) < tol.newton){
                ##printf("inewton = %d, tau = %e\n", inewton, tau)
                break
            }
        }
        tau.pre = tau

        
        rho.tmp.vec = rho.new.vec
        rho.new.vec = FuncRho(tau, sigma.vec, m.vec, L, epsilon)
        rho.new.vec = LineSearch(rho.tmp.vec, rho.new.vec, D.vec, R.mat, sigma.vec, L, epsilon)
        
        lem.new = FuncLem(rho.new.vec, D.vec, R.mat, sigma.vec, L)
        diff.lem = lem.new - lem
        kldiv = KLDiv(rho.tmp.vec, rho.new.vec, R.mat)

        printf("iem = %d, kldiv = %e, diff.lem = %e, lem = %e, lem.new = %e \n", iem, kldiv, diff.lem, lem, lem.new)
        if(diff.lem > 0){
            flag.good = 0
            break
        }
        if(kldiv < 0.0){
            flag.good = 0
            break
        }
        ##if(abs(diff.lem) < 1.e-3){
        ##    break
        ##}

        ## line search on: 1.e-8
        ## line search off : 1.e-10
        if(kldiv < 1.e-10){
            break
        }
    }
    ans = list(rho.new.vec, flag.good)
    return(ans)
}


LineSearch <- function(x.vec, x.new.vec, D.vec, R.mat, sigma.vec, L, epsilon){
    x1     = x.vec[1]
    x1.new = x.new.vec[1]
    x2.vec = x.vec[2:length(x.vec)]
    x2.new.vec = x.new.vec[2:length(x.new.vec)]
    
    theta.vec = log( x2.vec / (1 - x1))
    theta.new.vec = log( x2.new.vec / (1 - x1.new))
    nstep = 100
    lem.init = FuncLem(x.vec, D.vec, R.mat, sigma.vec, L)
    lem.pre  = lem.init
    x.pre.vec = x.vec

    outfile = sprintf("temp.dat")
    fopen(outfile)
    fprintf(outfile, "skip sing\n")
    fprintf(outfile, "read\n")
    
    eta = 3.0
    for( istep in 1 : nstep){
        factor = eta**istep
        lx1.this = factor * (log(x1.new) - log(x1)) + log(x1)
        x1.this = exp(lx1.this)
        theta.this.vec = factor * (theta.new.vec - theta.vec) + theta.vec
        alpha = sum(exp(theta.this.vec))
        x2.this.vec = (1 - x1.this) * exp(theta.this.vec) / alpha
        x.this.vec = c(x1.this, x2.this.vec)

        if(min(x.this.vec) < epsilon){
            printf("min(x.this.vec) < epsilon: factor(istep) = %e (%d)\n", factor, istep)
            if(istep != 1){
                x.new.vec = x.pre.vec
            }
            break
        }

        lem = FuncLem(x.this.vec, D.vec, R.mat, sigma.vec, L)

        ##printf("linesearch: min(x.this.vec) = %e\n", min(x.this.vec))
        ##printf("linesearch: max(x.this.vec) = %e\n", max(x.this.vec))
        ##printf("linesearch: sum(x.this.vec) = %e\n", sum(x.this.vec))
        ##printf("linesearch: lem(x.this.vec) = %e\n", lem)

        fprintf(outfile, "%d  %e\n", istep, lem - lem.init)

        if(lem.pre < lem){
            ## printf("factor(istep) = %e (%d)\n", factor, istep)
            if(istep != 1){
                x.new.vec = x.pre.vec
            }
            break
        }
        lem.pre = lem
        x.pre.vec = x.this.vec
    }
    return(x.new.vec)
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

FuncM <- function(rho.vec, D.vec, R.mat){
    num.vec = R.mat %*% rho.vec
    ans.vec = t(R.mat) %*% (D.vec / num.vec) * rho.vec
    return(ans.vec)
}

FuncRho <- function(tau, sigma.vec, m.vec, L, epsilon){
    termb.vec = sigma.vec + tau / L
    ans.vec = ( termb.vec + sqrt( termb.vec * termb.vec + 4 * m.vec / L ) ) / 2.0
    ans.vec = mapply(ThresEpsilon, ans.vec, epsilon)
    return(ans.vec)
}

FuncDiffRhoVec <- function(tau, sigma.vec, m.vec, L, epsilon){
    tau.thres.vec = FuncTauThres(sigma.vec, m.vec, L, epsilon)
    ans.vec = mapply(FuncDiffRho, tau, sigma.vec, m.vec, L, epsilon, tau.thres.vec)
    return(ans.vec)
}

FuncDiffRho <- function(tau, sigma, m, L, epsilon, tau.thres){
    ans = 0.0
    if(tau <= tau.thres){
        ans = 0.0
    }
    else {
        termb = sigma + tau / L
        root = sqrt( termb * termb + 4 * m / L )
        ans = (termb + root) / (2 * L * root)
    }
    return(ans)
}

FuncTauThres <- function(sigma.vec, m.vec, L, epsilon){
    tau.thres.vec = L * (epsilon - sigma.vec) - m.vec / epsilon
    return(tau.thres.vec)
}

FuncS <- function(tau, sigma.vec, m.vec, L, epsilon){
    ans = sum( FuncRho(tau, sigma.vec, m.vec, L, epsilon) ) - 1
    return(ans)
}

FuncDiffS <- function(tau, sigma.vec, m.vec, L, epsilon){
    ans = sum(FuncDiffRhoVec(tau, sigma.vec, m.vec, L, epsilon))
    return(ans)
}


##################

FuncLem <- function(rho.vec, D.vec, R.mat, sigma.vec, L){
    term1 = L * sum( (rho.vec - sigma.vec) * (rho.vec - sigma.vec) ) / 2.0
    term2 = FuncG(rho.vec, D.vec, R.mat)
    ans = term1 + term2
    return(ans)
}

FuncL <- function(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol){
    term1 = FuncF(rho.vec, beta, mu, nrow, ncol)
    term2 = FuncG(rho.vec, D.vec, R.mat)
    ans = term1 + term2
    return(ans)
}

FuncL.supp <- function(rho.vec, D.vec, R.mat, beta, mu, nrow, ncol, rho.supp.vec){
    term1 = FuncF.supp(rho.vec, beta, mu, nrow, ncol, rho.supp.vec)
    term2 = FuncG(rho.vec, D.vec, R.mat)
    ans = term1 + term2
    return(ans)
}

QMinusF <- function(rho.new.vec, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)
{
    term1 =      FuncF(rho.vec, beta, mu, nrow, ncol)
    term2 = -1 * FuncF(rho.new.vec, beta, mu, nrow, ncol)
    term3 = sum( (rho.new.vec - rho.vec) * DiffF(rho.vec, beta, mu, nrow, ncol) )
    term4 = L / 2.0 * sum( (rho.new.vec - rho.vec) * (rho.new.vec - rho.vec) )
    ans = term1 + term2 + term3 + term4
    return (ans)
}

QMinusF.supp <- function(rho.new.vec, rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, epsilon)
{
    rho.supp.vec = GetSuppVec(rho.vec, epsilon)
    rho.new.supp.vec = GetSuppVec(rho.new.vec, epsilon)
    printf("Q-F: nzero = %d, nzero.new = %d\n", GetNZeroVec(rho.vec, epsilon), GetNZeroVec(rho.new.vec, epsilon))
    ## term1 =      FuncF(rho.vec, beta, mu, nrow, ncol)
    ## term2 = -1 * FuncF(rho.new.vec, beta, mu, nrow, ncol)
    term1 =      FuncF.supp(rho.vec, beta, mu, nrow, ncol, rho.new.supp.vec)
    term2 = -1 * FuncF.supp(rho.new.vec, beta, mu, nrow, ncol, rho.new.supp.vec)
    term3 = sum( (rho.new.vec - rho.vec) * DiffF.supp(rho.vec, beta, mu, nrow, ncol, rho.new.supp.vec) )
    term4 = L / 2.0 * sum( (rho.new.vec - rho.vec) * (rho.new.vec - rho.vec) )
    ans = term1 + term2 + term3 + term4
    return (ans)
}


FuncSigma <- function(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol){
    term1 = rho.vec
    term2 = -1.0 / L * DiffF(rho.vec, beta, mu, nrow, ncol)
    ans.vec = term1 + term2
    return(ans.vec)
}

FuncSigma.supp <- function(rho.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, rho.supp.vec){
    term1 = rho.vec
    term2 = -1.0 / L * DiffF.supp(rho.vec, beta, mu, nrow, ncol, rho.supp.vec)
    ans.vec = term1 + term2
    return(ans.vec)
}


FuncF <- function(rho.vec, beta, mu, nrow, ncol)
{
    term1 = (1.0 - beta) * sum(log(rho.vec))
    term2 = mu * TermV(rho.vec, nrow, ncol)
    ans = term1 + term2
    return(ans)
}

FuncF.supp <- function(rho.vec, beta, mu, nrow, ncol, rho.supp.vec)
{
    term1 = (1.0 - beta) * SumLogPlus(rho.vec, rho.supp.vec)
    term2 = mu * TermV(rho.vec, nrow, ncol)
    ans = term1 + term2
    return(ans)
}

DiffF <- function(rho.vec, beta, mu, nrow, ncol)
{
    term1 = (1 - beta) / rho.vec
    term2 = mu * DiffTermV(rho.vec, nrow, ncol)
    ans = term1 + term2
    return (ans)
}

DiffF.supp <- function(rho.vec, beta, mu, nrow, ncol, rho.supp.vec)
{
    term1 = (1 - beta) * InvVect(rho.vec, rho.supp.vec)
    term2 = mu * DiffTermV(rho.vec, nrow, ncol)
    ans = term1 + term2
    return (ans)
}


FuncG <- function(rho.vec, D.vec, R.mat)
{
    num.vec = R.mat %*% (rho.vec)
    ans = -1 * sum( D.vec * log( num.vec ) )
    return(ans)
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

ThresEpsilon <- function(xval, epsilon){
    ans = xval
    if(is.nan(xval)){
        ans = epsilon
    }
    else {
        if(xval < epsilon){
            ans = epsilon
        }
    }
    return(ans)
}

SumLogPlus <- function(vec, supp.vec){
    ans = 0.0
    for(index in 1 : length(vec)){
        if(supp.vec[index] > 0){
            ans = ans + log(vec[index])
        }
    }
    return(ans)
}

InvVect <- function(vec, supp.vec){
    len = length(vec)
    ans.vec = rep(0.0, len)
    for(index in 1 : len){
        if(supp.vec[index] > 0){
            ans.vec[index] = 1.0 / vec[index]
        }
        else {
            ans.vec[index] = 0.0
        }
    }
    return(ans.vec)
}

DivVect <- function(num.vec, den.vec, epsilon){
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

GetSuppVec <- function(vec, epsilon){
    supp.vec = rep(0, length(vec))
    for(index in 1 : length(vec)){
        if(vec[index] > epsilon){
            supp.vec[index] = 1
        }
    }
    return(supp.vec)
}

GetNZeroVec <- function(vec, epsilon){
    nzero = 0
    for(index in 1 : length(vec)){
        if(abs(vec[index]) <= epsilon){
            nzero = nzero + 1
        }
    }
    return(nzero)
}

