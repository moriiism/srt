#!/usr/local/bin/Rscript

###
### poiss.S.R
###
### 2017.11.21 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
srtdir    = "/home/morii/work/github/moriiism/srt"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(srtdir, "script/rlib/em.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
respdir     = args[1]
datafile    = args[2]

R.mat = LoadModel(respdir)
printf("ncol(R.mat) = %d\n", ncol(R.mat))
printf("nrow(R.mat) = %d\n", nrow(R.mat))

D.vec = LoadData(datafile)
printf("length(D.vec) = %d\n", length(D.vec))


ncol = 60
nrow = 60
nx.vec = ncol * nrow
x.vec = rep(1.0 / nx.vec, nx.vec)
x.vec = SolveByProxMap(x.vec, D.vec, R.mat, beta, mu, nrow, ncol, lin.or.log)


