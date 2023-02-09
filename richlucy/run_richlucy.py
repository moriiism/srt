#!/usr/bin/env python3
#
# run_richlucy.py
#
'''
Overview:
  run richlucy for many data
Input:

Output:

Data:

Details:

'''
import os
import sys
import subprocess
import math
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import time

srt_dir = "/home/morii/work/github/moriiism/srt"

def run_cmd(cmd):
    subprocess.call(cmd)

def LoadList1(list_file):
    line_lst = []
    list_fptr = open(list_file, "r")
    for line in list_fptr:
        line = line.rstrip()
        if(line == ""):
            continue
        if(line[0] == "#"):
            continue
        # print(line)
        line_lst.append(line)

    list_fptr.close()
    return(line_lst)

def LoadPhaseListFile(phase_list_file):
    phase_id_lst = []
    phase_st_lst = []
    phase_ed_lst = []
    phase_ratio_lst = []
    phase_tag_lst = []    
    phase_list_fptr = open(phase_list_file, "r")
    for line in phase_list_fptr:
        line = line.rstrip()
        if(line == ""):
            continue
        if(line[0] == "#"):
            continue
        print(line)
        (phase_id, phase_st, phase_ed,
         phase_ratio, phase_tag) = line.split()
        phase_id_lst.append(phase_id)
        phase_st_lst.append(phase_st)
        phase_ed_lst.append(phase_ed)
        phase_ratio_lst.append(phase_ratio)        
        phase_tag_lst.append(phase_tag)
    
    phase_list_fptr.close()
    print(phase_id_lst)
    print(phase_st_lst)
    print(phase_ed_lst)
    print(phase_ratio_lst)    
    print(phase_tag_lst)

    return (phase_id_lst, phase_st_lst, phase_ed_lst,
            phase_ratio_lst, phase_tag_lst)

#
# main
#

iarg = 1
work_dir       = sys.argv[iarg]; iarg += 1
target_list    = sys.argv[iarg]; iarg += 1
telescope_list = sys.argv[iarg]; iarg += 1
data_ver_list  = sys.argv[iarg]; iarg += 1
resp_ver_list  = sys.argv[iarg]; iarg += 1
ene_band_list  = sys.argv[iarg]; iarg += 1
rec_prog_list  = sys.argv[iarg]; iarg += 1
phase_list     = sys.argv[iarg]; iarg += 1
nskyx          = int(sys.argv[iarg]); iarg += 1
nskyy          = int(sys.argv[iarg]); iarg += 1
ndetx          = int(sys.argv[iarg]); iarg += 1
ndety          = int(sys.argv[iarg]); iarg += 1
cpu_num        = int(sys.argv[iarg]); iarg += 1

print("work_dir = ", work_dir)
print("target_list = ", target_list)
print("telescope_list = ", telescope_list)
print("data_ver_list = ", data_ver_list)
print("resp_ver_list = ", resp_ver_list)
print("ene_band_list = ", ene_band_list)
print("rec_prog_list = ", rec_prog_list)
print("phase_list = ", phase_list)
print("nskyx = ", nskyx)
print("nskyy = ", nskyy)
print("ndetx = ", ndetx)
print("ndety = ", ndety)
print("cpu_num = ", cpu_num)

# target object
# telescope name
# data version
# resp version
# energy band
# reconstruction program
target_lst = LoadList1(target_list)
telescope_lst = LoadList1(telescope_list)
data_ver_lst = LoadList1(data_ver_list)
resp_ver_lst = LoadList1(resp_ver_list)
ene_band_lst = LoadList1(ene_band_list)
rec_prog_lst = LoadList1(rec_prog_list)
(phase_id_lst, phase_st_lst, phase_ed_lst,
 phase_ratio_lst, phase_tag_lst) = LoadPhaseListFile(phase_list)

out_dir_lst = []
target_joint_lst = []
telescope_joint_lst = []
data_ver_joint_lst = []
resp_ver_joint_lst = []
ene_band_joint_lst = []
rec_prog_joint_lst = []
for target in target_lst:
    for telescope in telescope_lst:
        for data_ver in data_ver_lst:
            for resp_ver in resp_ver_lst:        
                for ene_band in ene_band_lst:
                    for rec_prog in rec_prog_lst:
                        out_dir = (work_dir + "/"
                                   + target + "/"
                                   + telescope + "/"
                                   + "data_" + data_ver + "/"
                                   + "resp_" + resp_ver + "/"
                                   + ene_band + "/"
                                   + rec_prog)
                        out_dir_lst.append(out_dir)
                        target_joint_lst.append(target)
                        telescope_joint_lst.append(telescope)
                        data_ver_joint_lst.append(data_ver)
                        resp_ver_joint_lst.append(resp_ver)
                        ene_band_joint_lst.append(ene_band)
                        rec_prog_joint_lst.append(rec_prog)

for out_dir in out_dir_lst:
    cmd = ["mkdir", "-p", out_dir]
    print(cmd)
    subprocess.call(cmd)

# ymaeda share data
ymaeda_share_dir = "/home/ymaeda/work/share/crab_hxt"

cmd_lst = []
for index in range(len(out_dir_lst)):
    print(out_dir_lst[index])
    data_ver = data_ver_joint_lst[index]
    resp_ver = resp_ver_joint_lst[index]    
    telescope = telescope_joint_lst[index]
    ene_band = ene_band_joint_lst[index]    

    # mkresp
    respdir = (ymaeda_share_dir + "/" + "resp" + "/"
               + resp_ver + "/" + telescope + "/"
               + ene_band + "/" + "gimage_100")
    outdir_mkresp = (out_dir_lst[index] + "/" + "resp")
    outfile_head = "rl"
    nphoton_input = 100
    cmd = [srt_dir + "/" + "mkresp" + "/" + "mkresp",
           respdir, outdir_mkresp, outfile_head,
           str(nskyx), str(nskyy), str(nphoton_input)]
    print(cmd)
    subprocess.call(cmd)

    for index_phase in range(len(phase_id_lst)):
        # richlucy
        datafile = (ymaeda_share_dir + "/"
                    + "data" + "/"
                    + data_ver + "/"
                    + telescope + "/"
                    + telescope + "_"
                    + "phase" + "_"
                    + phase_id_lst[index_phase] + "_"
                    + ene_band + "_"
                    + "x1178y1175_w80.fits")
        bg_file = "none"
        resp_norm_file = (out_dir_lst[index] + "/" + "resp" +
                          "/" + "rl_resp_norm.fits")
        eff_file = (out_dir_lst[index] + "/" + "resp" +
                    "/" + "rl_eff.fits")
        outdir_rl = (out_dir_lst[index] + "/" + "rec")
        outfile_head = "rl"
        nem = 1000
        tol_em = 1.0e-8
        acc_method = "squarem"

        cmd = [srt_dir + "/" + "richlucy" + "/" + "richlucy",
               datafile, bg_file, resp_norm_file, eff_file,
               str(nskyx), str(nskyy), str(ndetx), str(ndety),
               outdir_rl, outfile_head, str(nem), str(tol_em),
               acc_method]
        print(cmd)
        cmd_lst.append(cmd)

print("len of cmd_lst = ", len(cmd_lst))
for cmd in cmd_lst:
    print(cmd)

with ProcessPoolExecutor(max_workers=cpu_num) as executor:
    executor.map(run_cmd, cmd_lst)

    
#
#if __name__ == "__main__":
#    main()

