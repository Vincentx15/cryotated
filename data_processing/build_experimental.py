import os
import sys
import subprocess

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '..'))

from data_processing import mrc_utils

datadir_name = "../data/phenix"
files_list = os.listdir(datadir_name)

fail_list = []
for i, dirname in enumerate(files_list):
    if not i % 10:
        print("Done {}/{} files".format(i, len(files_list)))

    pdb_name, mrc = dirname.split("_")

    dir_path = os.path.join(datadir_name, dirname)
    pdb_path = os.path.join(datadir_name, dirname, "{}.pdb".format(pdb_name))
    mrcgz_path = os.path.join(datadir_name, dirname, f"emd_{mrc}.map.gz")
    mrc_path = os.path.join(datadir_name, dirname, f"emd_{mrc}.map")
    carved_name = os.path.join(datadir_name, dirname, f"{mrc}_carved.mrc")
    subsampled_name = os.path.join(datadir_name, dirname, f"{mrc}_subsampled.mrc")
    try:
        if not os.path.exists(mrc_path):
            cmd1 = f"wget -P {dir_path} ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-{mrc}/map/emd_{mrc}.map.gz"
            subprocess.call(cmd1, shell=True)
            cmd2 = f"gunzip {mrcgz_path}"
            subprocess.call(cmd2, shell=True)

        if not os.path.exists(carved_name):
            zoned = mrc_utils.carve(mrc=mrc_path, pdb_name=pdb_path, out_name=carved_name, filter_cutoff=6)
        if not os.path.exists(subsampled_name):
            subsampled = mrc_utils.subsample(mrc=carved_name, out_name=subsampled_name)
    except:
        fail_list.append(dirname)
# print(fail_list)
# print(len(fail_list))
# fail_list = ['6b23_7035', '3j9o_6266', '5adx_2857', '3j9x_6310', '3j9g_2699', '5flu_3222', '6b0x_7030']
# fail_list = [1, 1, 1, 1, 1, 1, 1]
# The failed examples consist in huge system for most. They have the pdb going out of the mrc, yet usually aligned.
# I choose to discard those as there are only 7/500 and because it would imply trimming the pdb for
# the further preprocessing.

