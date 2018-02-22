#!/usr/bin/Rscript

## #!/usr/local/bin/Rscript

###
### run.deconv.rl.S.R
###
### 2018.01.27 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
srtdir    = "/home/morii/work/github/moriiism/srt"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(srtdir, "script/rlib/rl.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
respdir     = args[1]
datafile    = args[2]
outfile     = args[3]

printf("respdir = %s", respdir)
printf("datafile = %s\n", datafile)
printf("outfile = %s\n", outfile)

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
x.vec = rep( sum(D.vec) / nx.vec, nx.vec)

x.vec = SolveByRL(x.vec, D.vec, R.mat, nrow, ncol)

array = array(x.vec, dim=c(ncol, nrow))
writeFITSim(array, file=outfile)

