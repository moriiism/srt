#!/usr/bin/Rscript

###
### run.cv.S.R
###
### 2018.03.05 M.Morii
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
cvdir       = args[2]

outfile     = args[5]
initfile    = args[6]

printf("respdir = %s", respdir)
printf("dvdir   = %s\n", dvdir)
printf("outfile = %s\n", outfile)
printf("initfile = %s\n", initfile)

dirname = dirname(outfile)
dir.create(dirname, recursive = TRUE)

R.mat = LoadModel(respdir)
printf("ncol(R.mat) = %d\n", ncol(R.mat))
printf("nrow(R.mat) = %d\n", nrow(R.mat))

D.vec = LoadData(datafile)
printf("length(D.vec) = %d\n", length(D.vec))

ncol = 60
nrow = 60
nx.vec = ncol * nrow

lin.or.log = "lin"

if(initfile == "no"){
    printf("no load file\n")
    x.vec = rep( sum(D.vec) / nx.vec, nx.vec)
} else {
    printf("load initfile ...\n")
    x.vec = LoadData(initfile)
}






array = array(x.vec, dim=c(ncol, nrow))
writeFITSim(array, file=outfile)

