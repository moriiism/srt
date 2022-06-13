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



cmd = ["mkdir", outdir + "/" + "smr"]
print(cmd)
subprocess.call(cmd)
mu_heldist_file = f"{outdir}/smr/mu_heldist.dat"
mu_heldist_file_fptr = open(mu_heldist_file, "w")
print("! mu heldist", file=mu_heldist_file_fptr)

mu_lst = [1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9]
# [1.0e3, 2e3, 6e3, 1.0e4, 2e4, 6e4, 1.0e5]
# [1.0e0, 1.0e1, 1.0e2, 1.0e3, 1.0e4,
# [1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6]
# [1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1]


for mu in mu_lst:
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
               resp_file, recfile, valfile, outdir_eval, outfile_head_eval]
        print(cmd)
        subprocess.call(cmd)

        
    # average helldist
    heldist_sum = 0.0
    for ifold in range(nfold):
        heldist = 0.0
        helldist_file = f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}/eval_heldist.txt"
        helldist_file_fptr = open(helldist_file, "r")
        for line in helldist_file_fptr:
            heldist = line.rstrip('\n')
        helldist_file_fptr.close()
        heldist_sum += float(heldist)
            
    heldist_ave = heldist_sum / nfold
    print(heldist_ave)
    print(f"{mu:.1e} {heldist_ave}", file=mu_heldist_file_fptr)

mu_heldist_file_fptr.close()


#
#if __name__ == "__main__":
#    main()

