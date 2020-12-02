import os
import numpy as np

from chimera import runCommand as rc  # use 'rc' as shorthand for runCommand
from chimera import replyobj  # for emitting status messages

# Get a resolution sample. Instead of reading the info in the file, we make a stochastic approximation
# The resolution mean is 3.736+/-0.439
# An alternative would be to parse the .Dat
resolution = np.random.normal(3.736, 0.439)

# loop through the files, opening, processing, and closing each in turn
datadir_name = "../data/phenix"
files_list = os.listdir(datadir_name)
for i, dirname in enumerate(files_list):
    if not i % 10:
        print("Done {}/{} files".format(i, len(files_list)))
    pdb_name, mrc = dirname.split("_")
    pdb_path = os.path.join(datadir_name, dirname, "{}.pdb".format(pdb_name))
    mrc_path = os.path.join(datadir_name, dirname, "{}_simulated.mrc".format(mrc))
    if os.path.exists(mrc_path):
        continue

    # replyobj.status("Processing " + fn)  # show what file we're working on
    rc("open {}".format(pdb_path))
    rc("molmap #0 {} gridSpacing 1 cutoffRange 5".format(resolution))
    rc("volume #0.1 save {}".format(mrc_path))
    rc("close all")

rc("stop now")
