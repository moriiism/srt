#!/usr/bin/env python3
#
# logbin.py
#
'''
Overview:
  print steps in log scale
Input:

Output:

Data:

Details:

'''
import os
import sys
import subprocess
#import pandas as pd
import math

iarg = 1
lo = float(sys.argv[iarg]); iarg += 1
up = float(sys.argv[iarg]); iarg += 1
nbin = int(sys.argv[iarg]); iarg += 1

print("lo = ", lo)
print("up = ", up)
print("nbin = ", nbin)

lo_log = math.log10(lo)
up_log = math.log10(up)
delta_log = (up_log - lo_log) / nbin

print("! list of edge of bin")
for ibin in range(nbin + 1):
    val_log = lo_log + delta_log * ibin
    val = pow(10, val_log)
    print(f"{val:.4e}")

print("! list of center of bin")
for ibin in range(nbin):
    val_log = lo_log + delta_log * (ibin + 0.5)
    val = pow(10, val_log)
    print(f"{val:.4e}")

#print(cmd)
#subprocess.call(cmd)

