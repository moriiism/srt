#!/usr/bin/Rscript

## #!/usr/local/bin/Rscript

###
### run.deconv.rl.ls.S.R
###
### 2018.03.29 M.Morii
###  R-L method with line search
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
respdir      = args[1]
datafile     = args[2]
outdir       = args[3]
outfile.head = args[4]
tol          = as.numeric(args[5])
nstep        = as.integer(args[6])

printf("respdir  = %s\n", respdir)
printf("datafile = %s\n", datafile)
printf("outdir   = %s\n", outdir)
printf("outfile.head  = %s\n", outfile.head)
printf("tol      = %e\n", tol)
printf("nstep    = %d\n", nstep)

dir.create(outdir, recursive = TRUE)

R.mat = LoadModel(respdir)
printf("ncol(R.mat) = %d\n", ncol(R.mat))
printf("nrow(R.mat) = %d\n", nrow(R.mat))

D.vec = LoadData(datafile)
printf("length(D.vec) = %d\n", length(D.vec))
Nph = sum(D.vec)
printf("Nph = %e\n", Nph)

ncol = 60
nrow = 60
nx.vec = ncol * nrow
x.vec = rep( 1.0 / nx.vec, nx.vec)
x.vec = SolveByRLLS(x.vec, Nph, D.vec, R.mat, tol, nstep, outdir, outfile.head)

outfile = sprintf("%s/%s.fits", outdir, outfile.head)
array = array(Nph * x.vec, dim=c(ncol, nrow))
writeFITSim(array, file=outfile)
