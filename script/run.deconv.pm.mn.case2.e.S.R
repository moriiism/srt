#!/usr/bin/Rscript

## #!/usr/local/bin/Rscript

###
### run.deconv.pm.mn.e.S.R
###
### 2018.04.09 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
srtdir    = "/home/morii/work/github/moriiism/srt"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(srtdir, "script/rlib/pm_mn_case2.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
respdir      = args[1]
datafile     = args[2]
mu           = as.numeric(args[3])
beta         = as.numeric(args[4])
outdir       = args[5]
outfile.head = args[6]
tol          = as.numeric(args[7])
nstep        = as.integer(args[8])
L            = as.numeric(args[9])
epsilon      = as.numeric(args[10])

printf("respdir = %s", respdir)
printf("datafile = %s\n", datafile)
printf("mu = %e\n", mu)
printf("beta = %e\n", beta)
printf("outdir = %s\n", outdir)
printf("outfile.head = %s\n", outfile.head)
printf("tol = %e\n", tol)
printf("nstep = %d\n", nstep)
printf("L = %e\n", L)
printf("epsilon = %e\n", epsilon)

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
x.vec = SolveByProxMapMN(x.vec, Nph, D.vec, R.mat, beta, mu, L, tol, nstep, outdir, outfile.head, nrow, ncol, epsilon)

outfile = sprintf("%s/%s.fits", outdir, outfile.head)
array = array(Nph * x.vec, dim=c(ncol, nrow))
writeFITSim(array, file=outfile)

warnings()
