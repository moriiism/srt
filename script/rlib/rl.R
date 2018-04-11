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



SolveByRL <- function(x.vec, D.vec, R.mat, outfile, tol, nstep){
    fopen(outfile)
    fprintf(outfile, "# k, kldiv, logl, logl.inc, delta.logl, tdiff\n") 
    k.max = nstep
    time.st = Sys.time()
    logl.init = FuncL(x.vec, R.mat, D.vec)
    logl.pre = logl.init

    Nph = sum(D.vec)
    for(k in 1 : k.max){
        den.vec = R.mat %*% x.vec
        x.new.vec = t(R.mat) %*% (D.vec / den.vec) * x.vec / Nph
        
        kldiv = KLDiv(x.vec, x.new.vec, R.mat)
        logl = FuncL(x.new.vec, R.mat, D.vec)
        delta.logl = logl - logl.pre
        logl.inc =  logl - logl.init
        time = Sys.time()
        tdiff = difftime(time, time.st, units=c("secs"))
        printf("k = %d, KLdiv = %e, min(den.vec) = %e, max(den.vec) = %e, tdiff = %e\n",
               k, kldiv, min(den.vec), max(den.vec), tdiff)
        fprintf(outfile, "%d  %e  %.10e  %e  %e  %e\n", k, kldiv, logl, logl.inc, delta.logl, tdiff)
        
        if(kldiv < tol){
            break
        }

        ## for next
        x.vec = x.new.vec
        logl.pre = logl
        
        ## save
        if(k %% 1000 == 0){
            ncol = 60
            nrow = 60
            outimg = sprintf("casa_1e+06_rl/outimg_%d.fits", k)
            array = array(x.vec, dim=c(ncol, nrow))
            writeFITSim(array, file=outimg)
        }
        
    }
    return (x.vec)
}


LineSearch <- function(x.vec, x.new.vec, R.mat, D.vec){
    x1     = x.vec[1]
    x1.new = x.new.vec[1]
    x2.vec = x.vec[2:length(x.vec)]
    x2.new.vec = x.new.vec[2:length(x.new.vec)]
    theta.vec = log( x2.vec / (1 - x1))
    theta.new.vec = log( x2.new.vec / (1 - x1.new))
    nstep = 100
    logl.init = FuncL(x.vec, R.mat, D.vec)
    logl.pre  = logl.init
    x.pre.vec = x.vec

    outfile = sprintf("temp.dat")
    fopen(outfile)
    fprintf(outfile, "skip sing\n")
    fprintf(outfile, "read\n")
    
    eta = 5.0
    for( istep in 1 : nstep){
        factor = eta**istep
        lx1.this = factor * (log(x1.new) - log(x1)) + log(x1)
        x1.this = exp(lx1.this)
        theta.this.vec = factor * (theta.new.vec - theta.vec) + theta.vec
        alpha = sum(exp(theta.this.vec))
        x2.this.vec = (1 - x1.this) * exp(theta.this.vec) / alpha
        x.this.vec = c(x1.this, x2.this.vec)

        logl = FuncL(x.this.vec, R.mat, D.vec)

        fprintf(outfile, "%d  %e\n", istep, logl - logl.init)
        
        if(logl.pre > logl){
            printf("factor(istep) = %e (%d)\n", factor, istep)
            x.new.vec = x.pre.vec
            break
        }
        logl.pre = logl
        x.pre.vec = x.this.vec
    }
    return(x.new.vec)
}

SolveByRLLS <- function(x.vec, Nph, D.vec, R.mat, tol, nstep, outdir, outfile.head){
    outfile = sprintf("%s/%s_moni.dat", outdir, outfile.head)
    fopen(outfile)
    fprintf(outfile, "# istep, kldiv, logl, logl.inc, delta.logl, tdiff\n")

    outfile.time.logl = sprintf("%s/%s_time_logl.dat", outdir, outfile.head)
    fopen(outfile.time.logl)
    fprintf(outfile.time.logl, "# tdiff logl.inc\n")

    time.st = Sys.time()
    logl.init = FuncL(x.vec, R.mat, D.vec)
    logl.pre = logl.init
    fprintf(outfile, "0  0  %.10e  0.0  0.0  0.0\n", logl.init)
    
    for(istep in 1 : nstep){
        den.vec = R.mat %*% x.vec
        x.new.vec = t(R.mat) %*% (D.vec / den.vec) * x.vec / Nph
        x.new.vec = LineSearch(x.vec, x.new.vec, R.mat, D.vec)

        kldiv = KLDiv(x.vec, x.new.vec, R.mat)
        logl = FuncL(x.new.vec, R.mat, D.vec)
        delta.logl = logl - logl.pre
        logl.inc =  logl - logl.init
        time = Sys.time()
        tdiff = difftime(time, time.st, units=c("secs"))
        printf("istep = %d, KLdiv = %e, delta.logl = %e, min(den.vec) = %e, max(den.vec) = %e, tdiff = %e\n",
               istep, kldiv, delta.logl, min(den.vec), max(den.vec), tdiff)
        fprintf(outfile, "%d  %e  %.10e  %e  %e  %e\n",
                istep, kldiv, logl, logl.inc, delta.logl, tdiff)
        fprintf(outfile.time.logl, "%e  %e\n", tdiff, logl.inc)
        
        if(kldiv < tol){
            break
        }

        ## for next
        x.vec = x.new.vec
        logl.pre = logl
        
        ## save
        if(istep %% 1000 == 0){
            ncol = 60
            nrow = 60
            outimg = sprintf("%s/outimg_%d.fits", outdir, istep)
            array = array(x.vec * Nph, dim=c(ncol, nrow))
            writeFITSim(array, file=outimg)
        }
        
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
    ans = sum( D.vec * log( num.vec ) )
    return(ans)
}
