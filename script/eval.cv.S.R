#!/usr/bin/Rscript

###
### eval.cv.S.R
###
### 2018.03.06 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
srtdir    = "/home/morii/work/github/moriiism/srt"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(srtdir, "script/rlib/cv.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
respdir     = args[1]
mu.list     = args[2]
beta.list   = args[3]
nfold       = as.integer(args[4])

printf("respdir = %s\n", respdir)
printf("mu.list = %s\n", mu.list)
printf("beta.list = %s\n", beta.list)
printf("nfold = %d\n", nfold)

R.mat = LoadModel(respdir)
printf("ncol(R.mat) = %d\n", ncol(R.mat))
printf("nrow(R.mat) = %d\n", nrow(R.mat))

ncol = 60
nrow = 60
nx.vec = ncol * nrow
lin.or.log = "lin"

df.mulist = read.table(mu.list)
df.betalist = read.table(beta.list)
nmu   = nrow(df.mulist)
nbeta = nrow(df.betalist)

nrow = nmu * nbeta
ncol = 4
mat.ans = matrix(0.0, nrow=nrow, ncol=4)


irow = 1
for( imu in 1 : nmu){
    for( ibeta in 1 : nbeta){
        mu   = as.numeric(df.mulist[imu, 1])
        beta = as.numeric(df.betalist[ibeta, 1])
        cvfile.list = sprintf("setup/cvfile_mu%.1e_beta%.1e.list", mu, beta)
        ans = EvalSqErr(cvfile.list, nfold, R.mat)
        ave = as.numeric(ans[1])
        sigma = as.numeric(ans[2])
        mat.ans[irow, 1] = mu
        mat.ans[irow, 2] = beta
        mat.ans[irow, 3] = ave
        mat.ans[irow, 4] = sigma
        irow = irow + 1
    }
}

print(mat.ans)
