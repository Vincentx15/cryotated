# Data generation

## Getting raw data

There are two avenues for obtaining data :

- The first one involves using simulated maps.
The method we chose for doing so is the molmap function
used in Chimera. This method relies on a chimera python
script that iterates through a set of pdb files and creates
a set of corresponding simulated maps.
Run ```chimera --nogui build_simulated.py```

- The second method involves experimental maps. 
Once we download the experimental maps and their corresponding pdb,
We first 'carve' the mrc to get a box around the pdb to have lighter mrc files.
We then optionally filter the values far away from the PDB
Finally we need to resample the experimental maps to get a fixed voxel_size value of 1.
This is done in the build_experimental.py

The output of those two steps are mrc files centered around
a protein and their corresponding pdbs.


## Processing data

Then there is a common processing pipeline to refine the
maps and normalize them. The goal is to normalize the values
and to threshold the maps to zero out some parts.
Then we might want to additionally chunk the data.
This is done in the refine_data.py script

Once we have those maps and their corresponding pdb, we
need to also create output maps, that use the same grid sizes
than the input maps and other annotations, such as the CA
position or the vector going from one CA to the next.
This is done in the Complex.py script.


## Loading data

There are two avenues that need to be explored :
- creating an hdf5 from those processed files and using mostly
the same loading process as the one we had before
- Simply using the loading process of the files.

The efficiency of the first could avoid overloading the I/O
but the data loading v_0 should at least try the second way
and assess if it limiting in terms of time.
A guess is that since the 3D conv models are really expensive,
the loading time might end up being negligible.