#!/usr/bin/Rscript

## #!/usr/local/bin/Rscript

###
### run.deconv.v2.test.S.R
###
### 2018.03.07 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
srtdir    = "/home/morii/work/github/moriiism/srt"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(srtdir, "script/rlib/empm_v2_test.R", sep="/") )

###
###
###

args <- commandArgs(TRUE)
respdir     = args[1]
datafile    = args[2]
mu          = as.numeric(args[3])
beta        = as.numeric(args[4])
outfile     = args[5]
initfile    = args[6]
L           = as.numeric(args[7])

printf("respdir = %s", respdir)
printf("datafile = %s\n", datafile)
printf("mu = %e\n", mu)
printf("beta = %e\n", beta)
printf("outfile = %s\n", outfile)
printf("initfile = %s\n", initfile)
printf("L = %e\n", L)

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

if(initfile == "no"){
    printf("no load file\n")
    x.vec = rep( sum(D.vec) / nx.vec, nx.vec)
} else {
    printf("load initfile ...\n")
    x.vec = LoadData(initfile)
}

x.vec = SolveByProxMap(x.vec, D.vec, R.mat, beta, mu, L, nrow, ncol)

array = array(x.vec, dim=c(ncol, nrow))
writeFITSim(array, file=outfile)

warnings()
