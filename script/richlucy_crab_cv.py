#!/usr/bin/env python3
#
# richlucy_crab_cv.py
#
'''
Overview:

Input:

Output:

Data:

Details:

'''
import os
import sys
import subprocess
import pandas as pd
import math

iarg = 1
orgfile        = sys.argv[iarg]; iarg += 1
rand_seed      = int(sys.argv[iarg]); iarg += 1
nfold          = int(sys.argv[iarg]); iarg += 1
outdir         = sys.argv[iarg]; iarg += 1
outfile_head   = sys.argv[iarg]; iarg += 1
nskyx          = int(sys.argv[iarg]); iarg += 1
nskyy          = int(sys.argv[iarg]); iarg += 1
ndetx          = int(sys.argv[iarg]); iarg += 1
ndety          = int(sys.argv[iarg]); iarg += 1
respdir        = sys.argv[iarg]; iarg += 1
posx_point_src = int(sys.argv[iarg]); iarg += 1
posy_point_src = int(sys.argv[iarg]); iarg += 1
nem            = int(sys.argv[iarg]); iarg += 1
tol_em         = float(sys.argv[iarg]); iarg += 1
npm            = int(sys.argv[iarg]); iarg += 1
tol_pm         = float(sys.argv[iarg]); iarg += 1
nnewton        = int(sys.argv[iarg]); iarg += 1
tol_newton     = float(sys.argv[iarg]); iarg += 1
mu_list_file   = sys.argv[iarg]; iarg += 1

print("orgfile = ", orgfile)
print("rand_seed = ", rand_seed)
print("nfold = ", nfold)
print("outdir = ", outdir)
print("outfile_head = ", outfile_head)
print("nskyx = ", nskyx)
print("nskyy = ", nskyy)
print("ndetx = ", ndetx)
print("ndety = ", ndety)
print("respdir = ", respdir)
print("posx_point_src = ", posx_point_src)
print("posy_point_src = ", posy_point_src)
print("nem = ", nem)
print("tol_em = ", tol_em)
print("npm = ", npm)
print("tol_pm = ", tol_pm)
print("nnewton = ", nnewton)
print("tol_newton = ", tol_newton)
print("mu_list_file = ", mu_list_file)

# make observation image for cross-validation
outdir_mkobs_cv = outdir + "/" + "obs_cv"
outfile_head_mkobs_cv = outfile_head
cmd = ["/home/morii/work/github/moriiism/srt/mkobs_cv/mkobs_cv",
       orgfile, str(rand_seed), str(nfold), outdir_mkobs_cv,
       outfile_head_mkobs_cv]
print(cmd)
subprocess.call(cmd)


# make response matrix, normed response matrix, 
# and efficiency matrix files
outdir_resp = outdir + "/" + "resp"
outfile_head_resp = "hxt"
nphoton_input_resp = 100
cmd = ["/home/morii/work/github/moriiism/srt/mkresp/mkresp", 
       respdir, outdir_resp, outfile_head_resp,
       str(nskyx), str(nskyy), str(nphoton_input_resp)]
print(cmd)
subprocess.call(cmd)


# mkimg_points
cmd = ["mkdir", outdir + "/" + "skysrc"]
print(cmd)
subprocess.call(cmd)
point_src_dat_file = outdir + "/" + "skysrc" + "/" + "point_src.dat"
point_src_dat_file_fptr = open(point_src_dat_file, "w")
print(f"{posx_point_src} {posy_point_src} 1.0", file=point_src_dat_file_fptr)
point_src_dat_file_fptr.close()

outdir_point_src = outdir + "/" + "skysrc"
outfile_head_point_src = "point_src"
cmd = ["/home/morii/work/github/moriiism/srt/mkimg_points/mkimg_points",
       point_src_dat_file, outdir_point_src, outfile_head_point_src,
       str(nskyx), str(nskyy)]
print(cmd)
subprocess.call(cmd)

# mu_list
mu_lst = []
mu_list_file_fptr = open(mu_list_file, "r")
for line in mu_list_file_fptr:
    mu = line.rstrip('\n')
    mu_lst.append(mu)
mu_list_file_fptr.close()

print(mu_lst)

for mu in mu_lst:
    mu = float(mu)
    for ifold in range(nfold):
        print("ifold = ", ifold)
        datafile = (f"{outdir}/obs_cv/{outfile_head}_obs_{rand_seed:04}_" +
                    f"{nfold:02}fold{ifold:02}_tr.fits")
        fixed_src_norm_file = (f"{outdir}/skysrc/point_src_norm.fits")
        skyfile = "none"
        resp_file = f"{outdir}/resp/hxt_resp.fits"
        eff_file = f"{outdir}/resp/hxt_eff.fits"
        outdir_ifold = f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}"
        cmd = ["/home/morii/work/github/moriiism/srt/richlucy_crab/richlucy_crab",
               datafile, fixed_src_norm_file, skyfile, resp_file, eff_file,
               str(nskyx), str(nskyy), str(ndetx), str(ndety),
               outdir_ifold, outfile_head,
               str(nem), str(tol_em), str(npm), str(tol_pm),
               str(nnewton), str(tol_newton), str(mu)]
        print(cmd)
        subprocess.call(cmd)


cmd = ["mkdir", outdir + "/" + "smr"]
print(cmd)
subprocess.call(cmd)
mu_rmse_file = f"{outdir}/smr/mu_rmse.qdp"
mu_rmse_file_fptr = open(mu_rmse_file, "w")
print("skip sing", file=mu_rmse_file_fptr)
print("read serr 2", file=mu_rmse_file_fptr)
print("! mu  rmse_ave  rmse_stddev", file=mu_rmse_file_fptr)

for mu in mu_lst:
    mu = float(mu)
    # eval
    for ifold in range(nfold):
        print("ifold = ", ifold)
        resp_file = f"{outdir}/resp/hxt_resp.fits"
        recfile = (f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}/{outfile_head}_pulsar+nebula.fits")
        valfile = (f"{outdir}/obs_cv/{outfile_head}_obs_{rand_seed:04}_" +
                   f"{nfold:02}fold{ifold:02}_vl.fits")
        outdir_eval = (f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}")
        outfile_head_eval = "eval"
        
        cmd = ["/home/morii/work/github/moriiism/srt/richlucy_crab/eval_val", 
               resp_file, recfile, valfile, str(nfold), outdir_eval, outfile_head_eval]
        print(cmd)
        subprocess.call(cmd)

        
    # average rmse
    rmse_sum = 0.0
    rmse_sum2 = 0.0
    nfold_exist = 0
    for ifold in range(nfold):
        rmse = 0.0
        rmse_file = f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}/eval_rmse.txt"
        if os.path.isfile(rmse_file) == False:
            print(f"warning: {rmse_file} does not exist.")
            continue
        else:
            nfold_exist += 1
        rmse_file_fptr = open(rmse_file, "r")
        for line in rmse_file_fptr:
            rmse = line.rstrip('\n')
        rmse_file_fptr.close()
        rmse_sum += float(rmse)
        rmse_sum2 += float(rmse) * float(rmse)

    if nfold_exist == 0:
        print(f"warning: nfold_exist == 0")
    else:
        if nfold_exist != nfold:
            print(f"warning: nfold_exist = {nfold_exist}, nfold = {nfold}")
        rmse_ave = rmse_sum / nfold_exist
        rmse_var = (rmse_sum2 - rmse_sum * rmse_sum / nfold_exist) / nfold_exist
        rmse_stddev = math.sqrt(rmse_var)
        print(f"{mu:.1e} {rmse_ave} {rmse_stddev}")
        print(f"{mu:.1e} {rmse_ave} {rmse_stddev}", file=mu_rmse_file_fptr)


print("log x on", file=mu_rmse_file_fptr)
print("lw 7", file=mu_rmse_file_fptr)
print("line on", file=mu_rmse_file_fptr)
print("ma 6 on", file=mu_rmse_file_fptr)
mu_rmse_file_fptr.close()


#
#if __name__ == "__main__":
#    main()

