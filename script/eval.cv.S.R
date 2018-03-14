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
outdir      = args[5]
outfile.head = args[6]

printf("respdir = %s\n", respdir)
printf("mu.list = %s\n", mu.list)
printf("beta.list = %s\n", beta.list)
printf("nfold = %d\n", nfold)
printf("outdir = %s\n", outdir)
printf("outfile.head = %s\n", outfile.head)

dir.create(outdir, recursive = TRUE)

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
        cvfile.list = sprintf("setup/cvfile_mu%.1e_beta%.2e.list", mu, beta)
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

df.ans = data.frame(mat.ans)
outfile = sprintf("%s/%s.dat", outdir, outfile.head)
colnames(df.ans) = c("mu","beta","ave","sigma")
write.table(df.ans, file=outfile, row.names=F)

##
## for hist
##

nrow.mu.plus = 0
for( imu in 1 : nmu){
    for( ibeta in 1 : nbeta){
        mu   = as.numeric(df.mulist[imu, 1])
        if(mu < 1e-10){
            next
        }
        nrow.mu.plus = nrow.mu.plus + 1
    }
}

mat.ans2 = matrix(0.0, nrow=nrow.mu.plus, ncol=4)

irow = 1
for( imu in 1 : nmu){
    for( ibeta in 1 : nbeta){
        mu   = as.numeric(df.mulist[imu, 1])
        beta = as.numeric(df.betalist[ibeta, 1])
        cvfile.list = sprintf("setup/cvfile_mu%.1e_beta%.2e.list", mu, beta)
        ans = EvalSqErr(cvfile.list, nfold, R.mat)
        ave = as.numeric(ans[1])
        sigma = as.numeric(ans[2])
        if(mu < 1e-10){
            next
        }
        mat.ans2[irow, 1] = log10(mu)
        mat.ans2[irow, 2] = beta
        mat.ans2[irow, 3] = ave
        mat.ans2[irow, 4] = sigma
        irow = irow + 1
    }
}

df.ans2 = data.frame(mat.ans2)
outfile2 = sprintf("%s/%s_mu_plus.dat", outdir, outfile.head)
write.table(df.ans2, file=outfile2, row.names=F, col.names=F)
