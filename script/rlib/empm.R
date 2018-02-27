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


SolveByProxMap <- function(x.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, lin.or.log){
    eta = 1.2
    y.vec = x.vec
    x.pre.vec = x.vec
    t = 1
    cost = 0.0
    cost.pre = 0.0
    tol = 1e-10
    k.max = 100
   
    ## cost = FuncFG(y.vec, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
    ## printf("SolveByProxMap:k = 0, cost = %e\n", cost)

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

    
    L.pre = L
    for(k in 1 : k.max){
        printf("SolveByProxMap: k = %d\n", k)

        array = array(y.vec, dim=c(ncol, nrow))
        # file.tmp = sprintf("temp%3.3d.fits", k)
        # writeFITSim(array, file=file.tmp)
        
        ik = FindIk(y.vec, L.pre, eta, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        printf("SolveByProxMap: ik = %d\n", ik)
        L = eta**ik * L.pre
        # L.pre = eta**(ik - 1) * L.pre
        L.pre = L 
        x.vec = ProxMap(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        t.new = (1.0 + sqrt(1.0 + 4.0 * t**2)) / 2.0
	y.new.vec = x.vec + (t - 1.0) / t.new * (x.vec - x.pre.vec)
	if(0 > prod(sign(y.new.vec))){
	     printf("***** sign minus ****\n")
	      y.new.vec = x.vec
	      t.new = t
 	}
        else{
             printf("#### sign plus ####\n")
        }

        ## y.new.vec = mapply(max, y.new.vec, 0.0)
        
        ## cost = FuncFG(y.new.vec, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
        printf("SolveByProxMap:k = %d, L = %e\n", k, L)

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
    ans = sum(y.vec) + mu * TermV(y.vec, nrow, ncol, lin.or.log)
    return(ans)
}

FuncG <- function(y.vec, R.mat, D.vec, beta, nrow, ncol)
{
    term1 = -1 * sum( D.vec * log( R.mat %*% y.vec ) )
    term2 = (1.0 - beta) * sum( log(y.vec) )
    term3 = (1.0 - beta) * nrow * ncol * log(sum(y.vec))
    ans = term1 + term2 + term3
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
        ans = max( ( sigma + sqrt( sigma * sigma + 4 * mval / L) ) / 2.0 , 0.0 )
    }
    else{
        ans = 0.0
    }
    return(ans)
}

# y.vec ---> y.new.vec
ProxMap <- function(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol, lin.or.log)
{
    printf("ProxMap: max(y.vec), min(y.vec) = %e, %e\n", max(y.vec), min(y.vec))
    
    y.new.vec = y.vec
    y.pre.vec = y.vec
    nem.step = 300
    tol = 1.0e-10
    for(iem.step in 1:nem.step){
        sigma.vec = y.new.vec - 1.0 / L * DiffF(y.new.vec, mu, nrow, ncol, lin.or.log)
        num.vec = R.mat %*% y.new.vec
        term1 = (t(R.mat) %*% (D.vec / num.vec)) * y.new.vec
        term2 = (1.0 - beta) * nrow * ncol * y.new.vec / sum(y.new.vec) - (1.0 - beta)
        m.vec = term1 + term2
        ### if(min(num.vec) < 1.0e-10){
        ###    break
        ###}
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
    diffF.vec = 1.0 + mu * DiffTermV(y.vec, nrow, ncol, lin.or.log)
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
    ans = sum( (log(mat) - log(mat.p1.p0)) * (log(mat) - log(mat.p1.p0)) )
    + sum( (log(mat) - log(mat.p0.p1)) * (log(mat) - log(mat.p0.p1)) )
    return (ans)
}

TermVlin <- function(y.vec, nrow, ncol)
{
    mat = matrix(y.vec, nrow=nrow, ncol=ncol)
    mat.aug.tmp = rbind(mat, mat[nrow,])
    mat.aug = cbind(mat.aug.tmp, mat.aug.tmp[,ncol])
    mat.p1.p0 = mat.aug[2:(nrow+1),1:ncol]
    mat.p0.p1 = mat.aug[1:nrow,2:(ncol+1)]
    ans = sum( (mat - mat.p1.p0) * (mat - mat.p1.p0) )
    + sum( (mat - mat.p0.p1) * (mat - mat.p0.p1) )
    return (ans)
}



