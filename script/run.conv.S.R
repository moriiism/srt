#!/usr/bin/Rscript

## #!/usr/local/bin/Rscript

###
### run.conv.S.R
###
### 2018.03.27 M.Morii
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
infile      = args[2]
outfile     = args[3]

printf("respdir = %s", respdir)
printf("infile = %s\n", infile)
printf("outfile = %s\n", outfile)

dirname = dirname(outfile)
dir.create(dirname, recursive = TRUE)

R.mat = LoadModel(respdir)
printf("ncol(R.mat) = %d\n", ncol(R.mat))
printf("nrow(R.mat) = %d\n", nrow(R.mat))


ncol = 60
nrow = 60
nx.vec = ncol * nrow
D.vec = LoadDataWithSize(infile, ncol, nrow)
printf("length(D.vec) = %d\n", length(D.vec))

x.vec = R.mat %*% D.vec

array = array(x.vec, dim=c(ncol, nrow))
writeFITSim(array, file=outfile)

D.norm.vec = D.vec / sum(D.vec)

logl = FuncL(D.norm.vec, R.mat, D.vec)
printf("logl = %.10e\n", logl)
