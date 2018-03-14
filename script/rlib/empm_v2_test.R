###
### empm_v2.R
###
### 2018.03.07 M.Morii
###   EM + prox-map (V_rho)
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

#    for(index in 1 : nrow){
#        if(sum(mat[index,]) < 1e-3){
#            printf("%e, %e, %e\n", min(mat[index,]), max(mat[index,]), sum(mat[index,]))
#        }
#    }
    
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


SolveByProxMap <- function(x.vec, D.vec, R.mat, beta, mu, L, nrow, ncol){
    eta = 1.2
    y.vec = x.vec
    x.pre.vec = x.vec
    t = 1
    cost = 0.0
    cost.pre = 0.0
    tol = 1e-20
    k.max = 100

    supp.vec = GetSuppVec(y.vec)
    cost = FuncFGSupp(y.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec)
    printf("SolveByProxMap:k = 0, cost = %e\n", cost)

    L.find = FindLByEM(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol)
    printf("SolveByProxMap: findLByEM: L.find = %e\n", L.find)
    L = L.find
    
    L.pre = L
    for(k in 1 : k.max){
        printf("SolveByProxMap: k = %d\n", k)

        ### dump
        array = array(y.vec, dim=c(ncol, nrow))
        file.tmp = sprintf("temp%3.3d.fits", k)
        writeFITSim(array, file=file.tmp)
        ### dump

        ik = FindIk(y.vec, L.pre, eta, R.mat, D.vec, beta, mu, nrow, ncol)
        L = eta**ik * L.pre
        printf("SolveByProxMap: ik = %d, L = %e, L.pre = %e\n", ik, L, L.pre)
        L.pre = L
        ans.proxmap = ProxMap(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol)
        x.vec = ans.proxmap[[1]]
        
        t.new = (1.0 + sqrt(1.0 + 4.0 * t**2)) / 2.0
	y.new.vec = x.vec + (t - 1.0) / t.new * (x.vec - x.pre.vec)
	if(1 > prod( (sign(y.new.vec) + 1)/2. )){
	     printf("***** sign minus ****\n")

             ## truncate
             y.new.trunc.vec = mapply(max, y.new.vec, 0.0)
             supp.vec = GetSuppVec(y.new.trunc.vec)
             cost.trunc = FuncFGSupp(y.new.trunc.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec)
             cost.org   = FuncFGSupp(x.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec)
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

        nzero.pre = GetNZeroVec(y.vec)
        nzero     = GetNZeroVec(y.new.vec)
        supp.vec = GetSuppVec(y.new.vec)
        cost.new = FuncFGSupp(y.new.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec)
        cost     = FuncFGSupp(y.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec)
        cost.diff = cost.new - cost
        cost.diff.ratio = cost.diff / cost
        kldiv = KLDiv(y.vec, y.new.vec, R.mat)
        printf("SolveByProxMap:k = %d, L = %e, cost.diff.ratio = %e, KL divergence = %e, nzero.pre, nzero = %d, %d\n", k, L, cost.diff.ratio, kldiv, nzero.pre, nzero)
        if(k > 1 && abs(cost.diff.ratio) < 1e-10){
            break
        }        
        if(k > 1 && kldiv < tol){
            break
        }
        
        t = t.new
        y.vec = y.new.vec
        x.pre.vec = x.vec
        cost.pre = cost
    }
    x.vec = y.vec
    return (x.vec)
}


FindLByEM <- function(x.vec, L, R.mat, D.vec, beta, mu, nrow, ncol){
    ik.max = 100
    eta = 1.2
    L.out = L
    
    for(ik in 0 : ik.max){
        L.new = eta**ik * L
        printf("FindLByEM: L.new = %e\n", L.new)

        y.new.vec = x.vec
        y.pre.vec = x.vec
        nem.step = 100
        tol = 1.0e-10

        sum.y.vec = sum(x.vec)
        flag.good = 1
        lem.pre = -1e10
        for(iem.step in 1:nem.step){
            sigma.vec = y.new.vec - 1.0 / L.new * DiffF(y.new.vec, mu, nrow, ncol)
            num.vec = R.mat %*% y.new.vec
            ## to avoid division of zero
            if(min(num.vec) < 1.0e-10){
                num.vec = mapply(AddEpsilon, num.vec)
                ## printf("FindLByEM: iem.step = %d: min(num.vec) < 1.e-10, then add epsilon\n", iem.step)
            }
            term1 = (t(R.mat) %*% (D.vec / num.vec)) * y.new.vec
            term2 = (1.0 - beta) * nrow * ncol * y.new.vec / sum(y.new.vec) - (1.0 - beta)
            m.vec = term1 + term2
            y.new.vec = mapply(Mfunc, m.vec, sigma.vec, L.new)


            ## check the likelihood value
            nzero.pre = GetNZeroVec(y.pre.vec)
            nzero     = GetNZeroVec(y.new.vec)
            supp.vec = GetSuppVec(y.new.vec)
            lem.pre = FuncLemSupp(y.pre.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec, L.new)
            lem     = FuncLemSupp(y.new.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec, L.new)
            diff.lem = lem - lem.pre
            kldiv = KLDiv(y.pre.vec, y.new.vec, R.mat)
            printf("FindLByEM: iem.step = %d, diff.lem = %e: lem.pre = %e, lem = %e, kldiv = %e, nzero.pre, nzero = %d, %d\n",
                   iem.step, diff.lem, lem.pre, lem, kldiv, nzero.pre, nzero)
            
            if(lem < lem.pre){
                printf("FindLByEM: iem.step = %d, diff.lem = %e: lem.pre = %e, lem = %e, nzero.pre, nzero = %d, %d\n",
                       iem.step, diff.lem, lem.pre, lem, nzero.pre, nzero)
                printf("FindLByEM: decreasing likelihood (lem(%e) < lem.pre(%e)\n", lem, lem.pre)
                flag.good = 0
                break
            }

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

### supp

FuncLemSupp <- function(y.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec, L)
{
    sigma.vec = y.vec - 1.0 / L * DiffF(y.vec, mu, nrow, ncol)
    term1 = -1.0 * L / 2. * sum( (y.vec - sigma.vec) * (y.vec - sigma.vec) * supp.vec)
    term2 = -1 * FuncGSupp(y.vec, R.mat, D.vec, beta, nrow, ncol, supp.vec)
    ans = term1 + term2
    return(ans)
}


FuncFGSupp <- function(y.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec)
{
    ans = FuncF(y.vec, mu, nrow, ncol) + FuncGSupp(y.vec, R.mat, D.vec, beta, nrow, ncol, supp.vec)
    return(ans)
}

FuncGSupp <- function(y.vec, R.mat, D.vec, beta, nrow, ncol, supp.vec)
{
    y.supp.vec = y.vec * supp.vec
    num.vec = R.mat %*% y.supp.vec
    num.vec = mapply(AddEpsilon, num.vec)
    
    term1 = -1 * sum( D.vec * log( num.vec ) )
    term2 = (1.0 - beta) * SumLogPlus(y.supp.vec)
    term3 = -1 * (1.0 - beta) * nrow * ncol * log(sum(y.supp.vec))
    term4 = sum(y.supp.vec)
    ans = term1 + term2 + term3 + term4
    ## printf("FuncGSupp = %e (%e, %e, %e, %e)\n", ans, term1, term2, term3, term4)
    return(ans)
}


### 

FuncFG <- function(y.vec, R.mat, D.vec, beta, mu, nrow, ncol)
{
    ans = FuncF(y.vec, mu, nrow, ncol) + FuncG(y.vec, R.mat, D.vec, beta, nrow, ncol)
    return(ans)
}

FuncF <- function(y.vec, mu, nrow, ncol)
{
    ans = mu * TermV(y.vec, nrow, ncol) / sum(y.vec)**2
    ## printf("FuncF = %e\n", ans)
    return(ans)
}

FuncG <- function(y.vec, R.mat, D.vec, beta, nrow, ncol)
{
    num.vec = R.mat %*% y.vec
    num.vec = mapply(AddEpsilon, num.vec)
    
    term1 = -1 * sum( D.vec * log( num.vec ) )
    term2 = (1.0 - beta) * SumLogPlus(y.vec)
    term3 = -1 * (1.0 - beta) * nrow * ncol * log(sum(y.vec))
    term4 = sum(y.vec)
    ans = term1 + term2 + term3 + term4
    printf("FuncG = %e (%e, %e, %e, %e)\n", ans, term1, term2, term3, term4)
    return(ans)
}

FindIk <- function(y.vec, L, eta, R.mat, D.vec, beta, mu, nrow, ncol){
    ik.max = 1000
    ik = 0
    L.org = L
    while(ik <= ik.max){
        L = eta**ik * L.org
        ans.proxmap = ProxMap(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol)
        pLy = ans.proxmap[[1]]
        ans.proxmap.flag.good = ans.proxmap[[2]]

        qminusf = QMinusF(pLy, y.vec, L, mu, nrow, ncol)
        printf("FindIk: (ik, qminusf) = (%d, %e)\n", ik, qminusf)
        if(qminusf >= 0  && ans.proxmap.flag.good == 1 ){
            break
        }
        ik = ik + 1
    }
    return(ik)
}


QMinusF <- function(y.new.vec, y.vec, L, mu, nrow, ncol)
{
    term1 = mu * ( TermV(y.new.vec, nrow, ncol) / sum(y.new.vec)**2 - TermV(y.vec, nrow, ncol) / sum(y.vec)**2 )
    term2 = sum( (y.new.vec - y.vec) * DiffF(y.vec, mu, nrow, ncol) )
    term3 = L / 2.0 * sum( (y.new.vec - y.vec) * (y.new.vec - y.vec) )
    ans = term1 + term2 + term3
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

GetNZeroVec <- function(vec){
    nzero = 0
    for(index in 1 : length(vec)){
        if(abs(vec[index]) < 1e-20){
            nzero = nzero + 1
        }
    }
    return(nzero)
}






# y.vec ---> y.new.vec
ProxMap <- function(y.vec, L, R.mat, D.vec, beta, mu, nrow, ncol)
{
    sum.y.vec = sum(y.vec)
    
    y.new.vec = y.vec
    y.pre.vec = y.vec
    nem.step = 200
    tol = 1.0e-20
    flag.good = 1
    for(iem.step in 1:nem.step){
       
        sigma.vec = y.new.vec - 1.0 / L * DiffF(y.new.vec, mu, nrow, ncol)
        num.vec = R.mat %*% y.new.vec
        ## to avoid division of zero
        if(min(num.vec) < 1.0e-10){
            num.vec = mapply(AddEpsilon, num.vec)
            ## printf("ProxMap: iem.step = %d: min(num.vec) < 1.e-10, then add epsilon\n", iem.step)
        }
        
        term1 = (t(R.mat) %*% (D.vec / num.vec)) * y.new.vec
        term2 = (1.0 - beta) * nrow * ncol * y.new.vec / sum(y.new.vec) - (1.0 - beta)
        m.vec = term1 + term2
        y.new.vec = mapply(Mfunc, m.vec, sigma.vec, L)


        nzero.pre = GetNZeroVec(y.pre.vec)
        nzero     = GetNZeroVec(y.new.vec)
        supp.vec = GetSuppVec(y.new.vec)
        lem.pre = FuncLemSupp(y.pre.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec, L)
        lem     = FuncLemSupp(y.new.vec, R.mat, D.vec, beta, mu, nrow, ncol, supp.vec, L)
        diff.lem = lem - lem.pre
        ##printf("ProxMap: iem.step = %d, diff.lem = %e: lem.pre = %e, lem = %e, nzero.pre, nzero = %d, %d\n",
        ##       iem.step, diff.lem, lem.pre, lem, nzero.pre, nzero)
        if(lem < lem.pre){
            printf("ProxMap: decreasing likelihood, iem.step = %d (lem(%e) < lem.pre(%e), diff.lem = %e)\n",
                   iem.step, lem, lem.pre, diff.lem)
            flag.good = 0
            break
        }
        
        kldiv = KLDiv(y.pre.vec, y.new.vec, R.mat)
        if(kldiv < tol){
            printf("ProxMap: iem.step = %d, KL divergence = %e\n", iem.step, kldiv)
            flag.good = 1
            break
        }

        if(sum.y.vec * 1.e3 < sum(y.new.vec)){
            printf("ProxMap: iem.step = %d, large sum(y.new.vec)\n", iem.step)
            flag.good = 0
            break
        }
        
        y.pre.vec = y.new.vec
    }
    ans = list(y.new.vec, flag.good)
    return(ans)
}

AddEpsilon <- function(xval){
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

DiffF <- function(y.vec, mu, nrow, ncol)
{
    term1 = DiffTermV(y.vec, nrow, ncol) / ( sum(y.vec)**2 )
    term2 = -2 * TermV(y.vec, nrow, ncol) / ( sum(y.vec)**3 ) * y.vec
    diffF.vec = mu * (term1 + term2)
    return (diffF.vec)
}

DiffTermV <- function(y.vec, nrow, ncol)
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

TermV <- function(y.vec, nrow, ncol)
{
    mat = matrix(y.vec, nrow=nrow, ncol=ncol)
    mat.aug.tmp = rbind(mat, mat[nrow,])
    mat.aug = cbind(mat.aug.tmp, mat.aug.tmp[,ncol])
    mat.p1.p0 = mat.aug[2:(nrow+1),1:ncol]
    mat.p0.p1 = mat.aug[1:nrow,2:(ncol+1)]
    ans = sum( (mat - mat.p1.p0) * (mat - mat.p1.p0) ) + sum( (mat - mat.p0.p1) * (mat - mat.p0.p1) )
    return (ans)
}
