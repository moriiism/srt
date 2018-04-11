###
### empm.R
###
### 2018.01.26 M.Morii
###   EM + prox-map
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

LoadDataWithSize <- function(file, nx, ny)
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
    D.org.vec = as.vector(fits$imDat)
    D.vec = rep(0.0, nx * ny)
    
    for(jbin in 1:ny){
        for(ibin in 1:nx){
            index.org = (jbin - 1) * naxis1 + ibin
            index = (jbin - 1) * nx + ibin
            D.vec[index] = D.org.vec[index.org]
        }
    }
    return (D.vec)
}




SolveByProxMap <- function(x.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, lin.or.log){
    eta = 1.2
    y.vec = x.vec
    x.pre.vec = x.vec
    t = 1
    cost = 0.0
    cost.pre = 0.0
    tol = 1e-10
    k.max = 100
   
    cost = FuncFG(y.vec, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
    printf("SolveByProxMap:k = 0, cost = %e\n", cost)

###    ## find proper L
###    L1 = FindL1(x.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
###    printf("FindL1 out: L1= %e\n", L1)
###    if(L1 < 0.0){
###        L = FindL2(x.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
###    }
###    else{
###        L = L1
###    }
###    printf("SOlveByProxMap: find L = %e\n", L)

    L = FindLByEM(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
    printf("SolveByProxMap: findLByEM: L = %e\n", L)
    
    L.pre = L
    for(k in 1 : k.max){
        printf("SolveByProxMap: k = %d\n", k)

        array = array(y.vec, dim=c(ncol, nrow))

        ### dump
        #file.tmp = sprintf("temp%3.3d.fits", k)
        #writeFITSim(array, file=file.tmp)
        ### dump
        
        ik = FindIk(y.vec, L.pre, eta, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        L = eta**ik * L.pre
        printf("SolveByProxMap: ik = %d, L = %e\n", ik, L)
        # L.pre = eta**(ik - 1) * L.pre
        L.pre = L 
        x.vec = ProxMap(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        t.new = (1.0 + sqrt(1.0 + 4.0 * t**2)) / 2.0
	y.new.vec = x.vec + (t - 1.0) / t.new * (x.vec - x.pre.vec)
	if(1 > prod( (sign(y.new.vec) + 1)/2. )){
	     printf("***** sign minus ****\n")

             ## truncate
             y.new.trunc.vec = mapply(max, y.new.vec, 0.0)
             cost.trunc = FuncFG(y.new.trunc.vec, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
             cost.org = FuncFG(x.vec, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
             printf("cost.trunc = %e, cost.org = %e\n", cost.trunc, cost.org)
             if(cost.trunc < cost.org){
                 printf("use truncated\n")
                 y.new.vec = y.new.trunc.vec
                 t.new = t
             }
             else {
                 printf("use not truncated\n")
                 y.new.vec = x.vec
                 t.new = t
             }
 	}
        else{
             printf("#### sign plus ####\n")
        }

        ## y.new.vec = mapply(max, y.new.vec, 0.0)
        
        cost = FuncFG(y.new.vec, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        printf("SolveByProxMap:k = %d, L = %e, cost = %e\n", k, L, cost)

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


FindLByEM <- function(x.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log){
    ik.max = 100
    eta = 1.2
    L.out = L
    
    for(ik in 0 : ik.max){
        L.new = eta**ik * L
        printf("FindLByEM: L.new = %e\n", L.new)

        y.new.vec = x.vec
        y.pre.vec = x.vec
        nem.step = 50
        tol = 1.0e-10

        sum.y.vec = sum(x.vec)
        flag.good = 1
        lval.pre = -1e10
        for(iem.step in 1:nem.step){
            sigma.vec = y.new.vec - 1.0 / L.new * DiffF(y.new.vec, mu, nrow, ncol, lin.or.log)
            num.vec = R.mat %*% y.new.vec

            ## to avoid division of zero
            if(min(num.vec) < 1.0e-10){
                printf("FindLByEM: iem.step = %d: min(num.vec) < 1.e-10, then break\n", iem.step)
                flag.good = 0
                break
            }
            
            term1 = (t(R.mat) %*% (D.vec / num.vec)) * y.new.vec
            term2 = (1.0 - beta) * nrow * ncol * y.new.vec / sum(y.new.vec) - (1.0 - beta)
            m.vec = term1 + term2

            ## check the likelihood value
            lval1 = sum(D.vec * log(num.vec))
            lval2 = -1.0 * L.new / 2. * sum( (y.new.vec - sigma.vec) * (y.new.vec - sigma.vec) )
            lval3 = -1 * (1.0 - beta) * SumLogPlus(y.new.vec)
            lval4 = (1.0 - beta) * nrow * ncol * log(sum(y.new.vec))
            lval5 = -1 * sum(y.new.vec)
            lval = lval1 + lval2 + lval3 + lval4 + lval5
            printf("FindLByEM: iem.step = %d: lval = %e, (%e, %e, %e, %e, %e) \n",
                   iem.step, lval, lval1, lval2, lval3, lval4, lval5)

            if(lval < lval.pre){
                printf("FindLByEM: decreasing likelihood (lval(%e) < lval.pre(%e)\n", lval, lval.pre)
                flag.good = 0
                break
            }
            lval.pre = lval
            
            y.new.vec = mapply(Mfunc, m.vec, sigma.vec, L.new)

            kldiv = KLDiv(y.pre.vec, y.new.vec, R.mat)
            ##        printf("KL divergence = %e\n", kldiv)
            ##        printf("iem.step = %d\n", iem.step)

            if(sum.y.vec * 2 < sum(y.new.vec)){
                printf("FindLByEM: large sum(y.new.vec)\n")
                flag.good = 0
                break
            }
            y.pre.vec = y.new.vec
        }
        if(flag.good == 1){
            L.out = L.new * 1.2
            break
        }
    }
    return(L.out)
}



FindL1 <- function(x.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log){
    L1 = L
    pLy = ProxMap(x.vec, L1, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
    qminusf1 = QMinusF(pLy, x.vec, L1, mu, nrow, ncol, lin.or.log)
    printf("FindL1: L1, qminusf1 = %e, %e\n", L1, qminusf1)

    L2 = L1 * 2
    pLy = ProxMap(x.vec, L2, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
    qminusf2 = QMinusF(pLy, x.vec, L2, mu, nrow, ncol, lin.or.log)
    printf("FindL1: L2, qminusf2 = %e, %e\n", L2, qminusf2)
    L.out = L1 - (L2 - L1) / (qminusf2 - qminusf1) * qminusf1

    return(L.out)
}

FindL2 <- function(x.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log){
    pLy = ProxMap(x.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
    qminusf = QMinusF(pLy, x.vec, L, mu, nrow, ncol, lin.or.log)
    printf("FindL2: qminusf = %e\n", qminusf)
    
    sign   = 0.0
    factor = 0.0
    if(qminusf >= 0.0){
        factor = 0.5
        sign = 1.0
    }
    else{
        factor = 2.0
        sign = -1.0
    }
    isearch.max = 1000
    isearch = 0
    L.this = L
    L.out = 0.0
    while(isearch <= isearch.max){
        L.this = L.this * factor
        pLy = ProxMap(x.vec, L.this, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        qminusf = QMinusF(pLy, x.vec, L.this, mu, nrow, ncol, lin.or.log)
        printf("FindL2: (isearch, L.this, qminusf) = (%d, %e, %e)\n", isearch, L.this, qminusf)

        if(L.this < 1.0e-10){
            break
        }
        if(qminusf * sign < 0.0){
            break
        }
        isearch = isearch + 1
    }
    if(sign > 0.0){
        L.out = L.this
    }
    else {
        L.out = L.this / factor
    }
    return(L.out)
}

FuncFG <- function(y.vec, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
{
    ans = FuncF(y.vec, mu, nrow, ncol, lin.or.log) + FuncG(y.vec, R.mat, D.vec, beta, nrow, ncol)
    return(ans)
}

FuncF <- function(y.vec, mu, nrow, ncol, lin.or.log)
{
    ans = mu * TermV(y.vec, nrow, ncol, lin.or.log)
    printf("FuncF = %e\n", ans)
    return(ans)
}

FuncG <- function(y.vec, R.mat, D.vec, beta, nrow, ncol)
{
    term1 = -1 * sum( D.vec * log( R.mat %*% y.vec ) )
    term2 = (1.0 - beta) * SumLogPlus(y.vec)
    term3 = -1 * (1.0 - beta) * nrow * ncol * log(sum(y.vec))
    term4 = sum(y.vec)
    ans = term1 + term2 + term3 + term4
    printf("FuncG = %e (%e, %e, %e, %e)\n", ans, term1, term2, term3, term4)
    return(ans)
}

FindIk <- function(y.vec, L, eta, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log){
    ik.max = 1000
    ik = 0
    while(ik <= ik.max){
        L = eta**ik * L
        pLy = ProxMap(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        qminusf = QMinusF(pLy, y.vec, L, mu, nrow, ncol, lin.or.log)
        printf("FindIk: (ik, qminusf) = (%d, %e)\n", ik, qminusf)
        if(qminusf >= 0){
            break
        }
        ik = ik + 1
    }
    return(ik)
}


QMinusF <- function(y.new.vec, y.vec, L, mu, nrow, ncol, lin.or.log)
{
    term1 = sum(y.new.vec - y.vec)
    term2 = mu * (TermV(y.new.vec, nrow, ncol, lin.or.log) - TermV(y.vec, nrow, ncol, lin.or.log))
    term3 = sum( (y.new.vec - y.vec) * DiffF(y.vec, mu, nrow, ncol, lin.or.log) )
    term4 = L / 2.0 * sum( (y.new.vec - y.vec) * (y.new.vec - y.vec) )
    ans = term1 + term2 + term3 + term4
    return (ans)
}

Mfunc <- function(mval, sigma, L){
    ans = 0.0
    if(mval >= 0){
        ans = max( ( sigma - 1./L + sqrt( (sigma - 1./L)**2 + 4 * mval / L) ) / 2.0 , 0.0 )
    }
    else{
        ans = 0.0
    }
    return(ans)
}

SumLogPlus <- function(vec){
    ans = 0.0
    for(index in 1 : length(vec)){
        if(vec[index] > 0.0){
            ans = ans + log(vec[index])
        }
    }
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



# y.vec ---> y.new.vec
ProxMap <- function(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
{
    printf("ProxMap: max(y.vec), min(y.vec), sum(y.vec) = %e, %e, %e\n", max(y.vec), min(y.vec), sum(y.vec))

    sum.y.vec = sum(y.vec)
    
    y.new.vec = y.vec
    y.pre.vec = y.vec
    nem.step = 300
    tol = 1.0e-10
    for(iem.step in 1:nem.step){
       
        sigma.vec = y.new.vec - 1.0 / L * DiffF(y.new.vec, mu, nrow, ncol, lin.or.log)
        num.vec = R.mat %*% y.new.vec

        ### to avoid division of zero
        if(min(num.vec) < 1.0e-10){
            printf("ProxMap: iem.step = %d: min(num.vec) < 1.e-10, then break\n", iem.step)
            break
        }
        
        term1 = (t(R.mat) %*% (D.vec / num.vec)) * y.new.vec
        term2 = (1.0 - beta) * nrow * ncol * y.new.vec / sum(y.new.vec) - (1.0 - beta)
        m.vec = term1 + term2

        ### printf("ProxMap: max(m.vec), min(m.vec), sum(m.vec) = %e, %e, %e\n", max(m.vec), min(m.vec), sum(m.vec))

        ## check the likelihood value
        lval1 = sum(D.vec * log(num.vec))
        lval2 = -1.0 * L / 2. * sum( (y.new.vec - sigma.vec) * (y.new.vec - sigma.vec) )
        lval = lval1 + lval2
        ## printf("ProxMap: iem.step = %d: lval = %e, lval1 = %e, lval2 = %e\n", iem.step, lval, lval1, lval2)
        

##        tmp.vec = (D.vec / num.vec)
##        array = array(tmp.vec, dim=c(ncol, nrow))
##        file.tmp = sprintf("norm%3.3d.fits", iem.step)
##        writeFITSim(array, file=file.tmp)
##        
##        array = array(num.vec, dim=c(ncol, nrow))
##        file.tmp = sprintf("num%3.3d.fits", iem.step)
##        writeFITSim(array, file=file.tmp)
##
##        array = array(sigma.vec, dim=c(ncol, nrow))
##        file.tmp = sprintf("sigma%3.3d.fits", iem.step)
##        writeFITSim(array, file=file.tmp)
##
##        array = array(term1, dim=c(ncol, nrow))
##        file.tmp = sprintf("term1%3.3d.fits", iem.step)
##        writeFITSim(array, file=file.tmp)
##
##        array = array(term2, dim=c(ncol, nrow))
##        file.tmp = sprintf("term2%3.3d.fits", iem.step)
##        writeFITSim(array, file=file.tmp)
##        
##        array = array(m.vec, dim=c(ncol, nrow))
##        file.tmp = sprintf("mvec%3.3d.fits", iem.step)
##        writeFITSim(array, file=file.tmp)
##
##        tmp.vec = (sigma.vec - 1./L) * (sigma.vec - 1./L) + 4 * m.vec / L
##        array = array(tmp.vec, dim=c(ncol, nrow))
##        file.tmp = sprintf("sqrt%3.3d.fits", iem.step)
##        writeFITSim(array, file=file.tmp)

        
        y.new.vec = mapply(Mfunc, m.vec, sigma.vec, L)

##        printf("ProxMap: max(y.new.vec), min(y.new.vec), sum(y.new.vec) = %e, %e, %e\n", max(y.new.vec), min(y.new.vec), sum(y.new.vec))

        kldiv = KLDiv(y.pre.vec, y.new.vec, R.mat)
##        printf("KL divergence = %e\n", kldiv)
##        printf("iem.step = %d\n", iem.step)
        if(kldiv < tol){
            printf("ProxMap: KL divergence = %e\n", kldiv)
            printf("ProxMap: iem.step = %d\n", iem.step)
            break
        }

        if(sum.y.vec * 1.e3 < sum(y.new.vec)){
            printf("ProxMap: large sum(y.new.vec)\n")
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

DiffF <- function(y.vec, mu, nrow, ncol, lin.or.log)
{
    diffF.vec = mu * DiffTermV(y.vec, nrow, ncol, lin.or.log)
    return (diffF.vec)
}

DiffTermV <- function(y.vec, nrow, ncol, lin.or.log)
{
    ans = 0.0
    if(lin.or.log == "lin"){
        ans = DiffTermVlin(y.vec, nrow, ncol)
    }
    else if(lin.or.log == "log"){
        ans = DiffTermVlog(y.vec, nrow, ncol)          
    }
    return (ans)
}

DiffTermVlog <- function(y.vec, nrow, ncol)
{
    mat = matrix(y.vec, nrow=nrow, ncol=ncol)
    mat.aug.p.tmp = rbind(mat, mat[nrow,])
    mat.aug.p = cbind(mat.aug.p.tmp, mat.aug.p.tmp[,ncol])
    mat.aug.m.tmp = rbind(mat[1,], mat)
    mat.aug.m = cbind(mat.aug.m.tmp[,1], mat.aug.m.tmp)
    
    mat.p.0 = mat.aug.p[2:(nrow+1),1:ncol]
    mat.0.p = mat.aug.p[1:nrow,2:(ncol+1)]
    mat.m.0 = mat.aug.m[1:nrow,2:(ncol+1)]
    mat.0.m = mat.aug.m[2:(nrow+1),1:ncol]

    mat.diff = 2.0 / mat * ( 4 * log(mat) - log(mat.m.0) - log(mat.p.0) - log(mat.0.m) - log(mat.0.p) )
    y.diff.vec = as.vector(mat.diff)
    return (y.diff.vec)
}

DiffTermVlin <- function(y.vec, nrow, ncol)
{
    mat = matrix(y.vec, nrow=nrow, ncol=ncol)
    mat.aug.p.tmp = rbind(mat, mat[nrow,])
    mat.aug.p = cbind(mat.aug.p.tmp, mat.aug.p.tmp[,ncol])
    mat.aug.m.tmp = rbind(mat[1,], mat)
    mat.aug.m = cbind(mat.aug.m.tmp[,1], mat.aug.m.tmp)
    
    mat.p.0 = mat.aug.p[2:(nrow+1),1:ncol]
    mat.0.p = mat.aug.p[1:nrow,2:(ncol+1)]
    mat.m.0 = mat.aug.m[1:nrow,2:(ncol+1)]
    mat.0.m = mat.aug.m[2:(nrow+1),1:ncol]

    mat.diff = 2.0 * ( 4 * mat - mat.m.0 - mat.p.0 - mat.0.m - mat.0.p )
    y.diff.vec = as.vector(mat.diff)
    return (y.diff.vec)
}

TermV <- function(y.vec, nrow, ncol, lin.or.log)
{
    ans = 0.0
    if("lin" == lin.or.log){
        ans = TermVlin(y.vec, nrow, ncol)
    }
    else if ("log" == lin.or.log) {
        ans = TermVlog(y.vec, nrow, ncol)
    }
    return (ans)
}

TermVlog <- function(y.vec, nrow, ncol)
{
    mat = matrix(y.vec, nrow=nrow, ncol=ncol)
    mat.aug.tmp = rbind(mat, mat[nrow,])
    mat.aug = cbind(mat.aug.tmp, mat.aug.tmp[,ncol])
    mat.p1.p0 = mat.aug[2:(nrow+1),1:ncol]
    mat.p0.p1 = mat.aug[1:nrow,2:(ncol+1)]
    ans = sum( (log(mat) - log(mat.p1.p0)) * (log(mat) - log(mat.p1.p0)) ) + sum( (log(mat) - log(mat.p0.p1)) * (log(mat) - log(mat.p0.p1)) )
    return (ans)
}

TermVlin <- function(y.vec, nrow, ncol)
{
    mat = matrix(y.vec, nrow=nrow, ncol=ncol)
    mat.aug.tmp = rbind(mat, mat[nrow,])
    mat.aug = cbind(mat.aug.tmp, mat.aug.tmp[,ncol])
    mat.p1.p0 = mat.aug[2:(nrow+1),1:ncol]
    mat.p0.p1 = mat.aug[1:nrow,2:(ncol+1)]
    ans = sum( (mat - mat.p1.p0) * (mat - mat.p1.p0) ) + sum( (mat - mat.p0.p1) * (mat - mat.p0.p1) )
    return (ans)
}



