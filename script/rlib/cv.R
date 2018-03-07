###
### cv.R
###
### 2018.03.06 M.Morii
###   CV
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

KLDiv <- function(re.vec, vl.vec, R.mat)
{
    q.re.vec = R.mat %*% re.vec
    q.re.vec = q.re.vec / sum(q.re.vec)
    q.vl.vec = vl.vec / sum(vl.vec)

    ## avoid zero component 
    ans = 0.0
    for( index in 1 : length(q.re.vec) ){
        if(q.re.vec[index] > 0.0){
            ans = ans + q.vl.vec[index] * log( q.vl.vec[index] / q.re.vec[index] )
        }
    }
    return (ans)
}


EvalSqErr <- function(cvfile.list, nfold, R.mat)
{
    df.cvfile = read.table(cvfile.list)
    sum = 0.0
    sum2 = 0.0
    for( ifold in 1 : nfold){
        printf("refile = %s\n", df.cvfile[ifold, 1])
        printf("vlfile = %s\n", df.cvfile[ifold, 2])
        re.img = as.character(df.cvfile[ifold, 1])
        vl.img = as.character(df.cvfile[ifold, 2])
    
        re.vec = LoadData(re.img)
        printf("length(re.vec) = %d\n", length(re.vec))
        vl.vec = LoadData(vl.img)
        printf("length(vl.vec) = %d\n", length(vl.vec))

        sqerr = SqErr(re.vec, vl.vec, R.mat)
        printf("sqerr = %e\n", sqerr)
        sum = sum + sqerr
        sum2 = sum2 + sqerr * sqerr
    }
    ave = sum / nfold
    sigma = (sum2 - nfold * ave*ave)/(nfold - 1)
    printf("ave = %e, sigma = %e\n", ave, sigma)
    ans = list(ave, sigma)
    return (ans)
}


SqErr <- function(re.vec, vl.vec, R.mat)
{
    q.re.vec = R.mat %*% re.vec
    q.re.vec = q.re.vec / sum(q.re.vec)
    q.vl.vec = vl.vec / sum(vl.vec)

    ans = sqrt( sum((q.vl.vec - q.re.vec) * (q.vl.vec - q.re.vec)) )
    return (ans)
}




