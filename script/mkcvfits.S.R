#!/usr/bin/Rscript

## #!/usr/local/bin/Rscript

###
### run.deconv.S.R
###
### 2017.12.01 M.Morii
###
###

mitooldir = "/home/morii/work/github/moriiism/mitool"
srtdir    = "/home/morii/work/github/moriiism/srt"

###
source( paste(mitooldir, "script/mirlib/iolib.R", sep="/") )
source( paste(srtdir, "script/rlib/empm.R", sep="/") )

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

printf("respdir = %s", respdir)
printf("datafile = %s\n", datafile)
printf("mu = %e\n", mu)
printf("beta = %e\n", beta)
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


# L = 1e-10
L = 1e-3
# L = 1.0
# L = 1.0e1
x.vec = SolveByProxMap(x.vec, D.vec, R.mat, beta, mu, L, nrow, ncol, lin.or.log)

array = array(x.vec, dim=c(ncol, nrow))
writeFITSim(array, file=outfile)

