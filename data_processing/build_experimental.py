import os
import sys
import subprocess

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '..'))

from data_processing import preprocess_experimental

datadir_name = "../data/phenix"
files_list = os.listdir(datadir_name)
for i, dirname in enumerate(files_list):
    if not i % 10:
        print("Done {}/{} files".format(i, len(files_list)))

    dirname = "3j1q_5415"
    pdb_name, mrc = dirname.split("_")

    dir_path = os.path.join(datadir_name, dirname)
    pdb_path = os.path.join(datadir_name, dirname, "{}.pdb".format(pdb_name))
    mrcgz_path = os.path.join(datadir_name, dirname, f"emd_{mrc}.map.gz")
    mrc_path = os.path.join(datadir_name, dirname, f"emd_{mrc}.map")
    carved_name = os.path.join(datadir_name, dirname, f"{mrc}_carved.mrc")
    subsampled_name = os.path.join(datadir_name, dirname, f"{mrc}_subsampled.mrc")

    if not os.path.exists(mrc_path):
        cmd1 = f"wget -P {dir_path} ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-{mrc}/map/emd_{mrc}.map.gz"
        subprocess.call(cmd1, shell=True)
        cmd2 = f"gunzip {mrcgz_path}"
        subprocess.call(cmd2, shell=True)

    if not os.path.exists(carved_name):
        zoned = preprocess_experimental.carve(mrc=mrc_path, pdb_name=pdb_path, out_name=carved_name, filter_cutoff=6)
    if not os.path.exists(subsampled_name):
        subsampled = preprocess_experimental.subsample(mrc=carved_name, out_name=subsampled_name)
