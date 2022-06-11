#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# richlucy_crab_cv.py
#
'''
Overview:

Input:
    $KUZE_ANA_DIR/input/record_path_ndraw.txt
Output:
    $KUZE_ANA_DIR/record_path_ndraw/
Data:
      $KUZE_ANA_DIR/data/edge_weight_add_bench_jig.csv
Details:

'''
import os
import sys
import subprocess
import pandas as pd

## check environment variable
#flag_set_env = 1
#if not "KUZE_PROG_DIR" in os.environ:
#    flag_set_env = 0
#if not "KUZE_DATA_DIR" in os.environ:
#    flag_set_env = 0
#if not "KUZE_ANA_DIR" in os.environ:
#    flag_set_env = 0
#if flag_set_env == 0:
#    print("Error: Set environment variables by source env.sh.")
#    exit()

#sys.path.append(os.environ["KUZE_PROG_DIR"])
#from mkpath.script.func.file import isOnlyLatinFile
#from mkpath.script.func.spec import Spec
#from mkpath.script.func.load_edge import loadEdgeRecordDF
#from mkpath.script.func.mk_path_fig import mkPathFigCSV_Record
#from mkpath.script.func.mk_path_df_record import (
#    getPathDFFromOptPath_Record)
#from mkpath.script.task.record_path.func import getRecordPath_NumBA

iarg = 1
orgfile       = sys.argv[iarg]; iarg += 1
rand_seed     = int(sys.argv[iarg]); iarg += 1
nfold         = int(sys.argv[iarg]); iarg += 1
outdir        = sys.argv[iarg]; iarg += 1
outfile_head  = sys.argv[iarg]; iarg += 1
nskyx         = int(sys.argv[iarg]); iarg += 1
nskyy         = int(sys.argv[iarg]); iarg += 1
ndetx         = int(sys.argv[iarg]); iarg += 1
ndety         = int(sys.argv[iarg]); iarg += 1
respdir       = sys.argv[iarg]; iarg += 1



print("orgfile = ", orgfile)
print("rand_seed = ", rand_seed)
print("nfold = ", nfold)
print("outdir = ", outdir)
print("outfile_head = ", outfile_head)


# 1. make response matrix, normed response matrix, 
#    and efficiency matrix files

outdir_resp = "resp"
outfile_head_resp = "arb"
nphoton_input_resp = 100
cmd = ["/home/morii/work/github/moriiism/srt/mkresp/mkresp", 
       respdir, outdir_resp, outfile_head_resp,
       str(nskyx), str(nskyy), str(nphoton_input_resp)]
print(cmd)
subprocess.call(cmd)


cmd = ["/home/morii/work/github/moriiism/srt/mkobs_cv/mkobs_cv",
       orgfile, str(rand_seed), str(nfold), outdir, outfile_head]
print(cmd)
subprocess.call(cmd)


for ifold in range(nfold):
    print(ifold)
    datafile = (f"{outdir}/{outfile_head}_obs_{rand_seed:04}_" +
                f"{nfold:02}fold{ifold:02}_tr.fits")
    cmd = ["/home/morii/work/github/moriiism/srt/richlucy_crab/richlucy_crab",
           datafile, fixed_src_norm_file, skyfile, resp_file, eff_file,
           nskyx, nskyy, ndetx, ndety, outdir, outfile_head,
           nem, tol_em, npm, tol_pm, nnewton, tol_newton, mu]
    print(cmd)
    subprocess.call(cmd)


# hxt_301_700_obs_0001_05fold00_tr.fits

#
#def exec_record_path_ndraw(exec_type: str,
#                           inputs):
#
#    task_name = "record_path_ndraw"
#    output_dir = None
#    input_df = None
#    if (exec_type == "term"):
#        # output_dir
#        output_dir = os.environ["KUZE_ANA_DIR"] + "/" + task_name
#        if not os.path.exists(output_dir):
#            os.makedirs(output_dir)
#        
#        # input_dir
#        input_dir = os.environ["KUZE_ANA_DIR"] + "/" + "input"
#        if not os.path.exists(input_dir):
#            os.makedirs(input_dir)
#            
#        input_file = task_name + ".txt"
#        input_file_full = input_dir + '/' + input_file
#        if not os.path.exists(input_file_full):
#            print(f"Error: Set input file at {input_dir}")
#            exit()
#        if not isOnlyLatinFile(input_file_full):
#            print(f"Error: Contains invalid characters at " +
#                  f"{input_file_full}")
#            exit()
#    
#        # read input file
#        input_df = pd.read_table(input_file_full, 
#                                 header=None,
#                                 sep=' *= *', 
#                                 engine='python', 
#                                 dtype=str, 
#                                 comment='#')
#        input_df.columns = ['name', 'val']
#        input_df.set_index("name", inplace=True)
#    
#    elif (exec_type == "api"):
#        from io import BytesIO
#        import base64
#
#        path_data = []
#        image_data = ""
#        message = ""
#
#        input_array = []
#        params = inputs.split('\n')
#        for param in params:
#            input_array.append(param.replace(' ', '').split('='))
#        input_df = pd.DataFrame(data=input_array)
#        input_df.columns = ['name', 'val']
#        input_df.set_index("name", inplace=True)
#
#    else:
#        print(f"Error: bad exec_type(={exec_type})")
#        exit()
#
#    # 工程数
#    ndraw = int(input_df.loc["ndraw", "val"])
#    # 材料外径
#    diam_org = float(input_df.loc["diam_org", "val"])
#    # 材料肉厚
#    thick_org = float(input_df.loc["thick_org", "val"])
#    # 引上外径
#    diam_prod = float(input_df.loc["diam_prod", "val"])
#    # 引上肉厚
#    thick_prod = float(input_df.loc["thick_prod", "val"])
#    # 鋼種CD
#    steel_cd = int(input_df.loc["steel_cd", "val"])
#    # 規格CD
#    spec_cd = int(input_df.loc["spec_cd", "val"])
#    # BA引抜回数の最小値
#    num_ba_min = int(input_df.loc["num_ba_min", "val"])
#    # BA引抜回数の最大値
#    num_ba_max = int(input_df.loc["num_ba_max", "val"])
#    # 全引抜をBAにする
#    flag_all_ba = int(input_df.loc["flag_all_ba", "val"])
#    # 治具制約 (1-6)
#    flag_jig = int(input_df.loc["flag_jig", "val"])
#    # SKD材の治具
#    flag_skd = int(input_df.loc["flag_skd", "val"])
#    # 中間工程として算出
#    flag_as_mid = int(input_df.loc["flag_as_mid", "val"])
#    
#    
#    # print param
#    print("ndraw = ", ndraw)
#    print("diam_org = ", diam_org)
#    print("thick_org = ", thick_org)
#    print("diam_prod = ", diam_prod)
#    print("thick_prod = ", thick_prod)
#    print("steel_cd = ", steel_cd)
#    print("spec_cd = ", spec_cd)
#    print("num_ba_min = ", num_ba_min)
#    print("num_ba_max = ", num_ba_max)
#    print("flag_all_ba = ", flag_all_ba)
#    print("flag_jig = ", flag_jig)
#    print("flag_skd = ", flag_skd)
#    print("flag_as_mid = ", flag_as_mid)
#    print("------")
#    
#    # check ndraw
#    if ndraw < 1:
#        message = f"Error: bad ndraw(={ndraw})."
#        if (exec_type == "term"):
#            print(message)
#            exit()
#        elif (exec_type == "api"):
#            return path_data, image_data, message
#    else:
#        pass
#    
#    if flag_all_ba == 0:
#        # check spec_cd
#        if spec_cd <= 0:
#            if ((num_ba_min < 0) or (num_ba_max < 0)):
#                message = "Error: spec_cd <= 0, " + \
#                          "so set 0<= num_ba_min <= num_ba_max."
#                if (exec_type == "term"):
#                    print(message)
#                    exit()
#                elif (exec_type == "api"):
#                    return path_data, image_data, message
#            else:
#                pass
#            if (num_ba_min > num_ba_max):
#                message = "Error: spec_cd <= 0, " + \
#                          "so set 0<= num_ba_min <= num_ba_max."
#                if (exec_type == "term"):
#                    print(message)
#                    exit()
#                elif (exec_type == "api"):
#                    return path_data, image_data, message
#            else:
#                pass
#        else:
#            if ((num_ba_min >= 0) or (num_ba_max >= 0)):
#                message = f"Error: spec_cd is positive(={spec_cd}), " + \
#                          "so set num_ba_min, num_ba_max < 0."
#                if (exec_type == "term"):
#                    print(message)
#                    exit()
#                elif (exec_type == "api"):
#                    return path_data, image_data, message
#            else:
#                pass
#            
#            spec = Spec()
#            num_ba_min = spec.getNumBA(spec_cd,
#                                       diam_prod)
#            num_ba_max = num_ba_min + 1
#    elif flag_all_ba == 1:
#        if ((spec_cd > 0) or
#            (num_ba_min >= 0) or
#            (num_ba_max >= 0)):
#            message = f"Error: flag_all_ba == 1, " + \
#                      "so set spec_cd <= 0, " + \
#                      "num_ba_min < 0 and num_ba_max < 0."
#            if (exec_type == "term"):
#                print(message)
#                exit()
#            elif (exec_type == "api"):
#                return path_data, image_data, message
#        else:
#            num_ba_min = ndraw
#            num_ba_max = ndraw
#    else:
#        message = f"Error: bad flag_all_ba(={flag_all_ba})"
#        if (exec_type == "term"):
#            print(message)
#            exit()
#        elif (exec_type == "api"):
#            return path_data, image_data, message
#            
#    print(f"num_ba_min = {num_ba_min}")
#    print(f"num_ba_max = {num_ba_max}")
#    
#    # check num_ba_max, ndraw
#    if num_ba_max > ndraw:
#        message = f"Error: num_ba_max(={num_ba_max}) must be " + \
#                  f"<= ndraw(={ndraw})."
#        if (exec_type == "term"):
#            print(message)
#            exit()
#        elif (exec_type == "api"):
#            return path_data, image_data, message
#    else:
#        pass
#    
#    # check num_ba_min, ndraw
#    if num_ba_min > ndraw:
#        message = f"Error: num_ba_min(={num_ba_min}) must be " + \
#                  f"<= ndraw(={ndraw})."
#        if (exec_type == "term"):
#            print(message)
#            exit()
#        elif (exec_type == "api"):
#            return path_data, image_data, message
#    else:
#        pass
#    
#    (edge_mid_ap_df,
#     edge_mid_ba_df,
#     edge_last_ap_df,
#     edge_last_ba_df) = loadEdgeRecordDF(steel_cd,
#                                         flag_jig,
#                                         flag_skd,
#                                         flag_as_mid)
#    (flag_ok,
#     opt_path_lst,
#     path_graph,
#     path_df) = getRecordPath_NumBA(edge_mid_ap_df,
#                                    edge_mid_ba_df,
#                                    edge_last_ap_df,
#                                    edge_last_ba_df,
#                                    diam_org,
#                                    thick_org,
#                                    diam_prod,
#                                    thick_prod,
#                                    num_ba_min,
#                                    num_ba_max,
#                                    flag_jig,
#                                    flag_skd,
#                                    ndraw)
#    if (flag_ok == 1):
#        (path_data,
#         image_data) = mkPathFigCSV_Record(exec_type,
#                                           opt_path_lst,
#                                           path_graph,
#                                           edge_mid_ap_df,
#                                           edge_mid_ba_df,
#                                           edge_last_ap_df,
#                                           edge_last_ba_df,
#                                           task_name,
#                                           output_dir,
#                                           flag_jig,
#                                           flag_skd,
#                                           flag_as_mid,
#                                           steel_cd)
#        message = "record_path_ndraw.py: successfully done."
#    else:
#        message = f"Failed to get optimal path."
#
#    if (exec_type == "term"):
#        print(message)
#        exit()
#    elif (exec_type == "api"):
#        return path_data, image_data, message
#
#
#if __name__ == "__main__":
#    main()
#
#
#

