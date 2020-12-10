import os
import sys

import pymol.cmd as cmd
import mrcfile
import numpy as np
import scipy.ndimage
import scipy.interpolate

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '..'))


def pymol_parse(pdbname):
    cmd.load(pdbname, "temp_parsing")
    xyz = cmd.get_coords('temp_parsing', 1)
    cmd.delete("temp_parsing")
    return xyz


'''
# Just to benchmark the parsing time against pymol.
# To use coordinates, pymol is 10 times faster.

from Bio.PDB import MMCIFParser
def bio_python(pdbname, parser):
    """
    Just to benchmark against
    """
    structure = parser.get_structure("poulet", pdbname)
    coords = [atom.get_vector() for atom in structure.get_atoms()]
    return coords



# 3.16s vs 0.27 for pymol
a = time.perf_counter()
for i in range(10):
    bio_python(pdbname, parser=parser)
print(f'time1 : {time.perf_counter() - a}')
'''


def load_mrc(mrc, mode='r+'):
    """
    returns an mrc from either a mrc or a mrc filename
    :param mrc:
    :param mode:
    :return:
    """
    if isinstance(mrc, str):
        mrc = mrcfile.open_async(mrc, mode=mode)
        return mrc.result()
    elif isinstance(mrc, mrcfile.mrcfile.MrcFile):
        return mrc
    else:
        raise ValueError("Wrong input to the MRC loading function")


"""
Once we download the experimental maps and their corresponding pdb,
We first 'carve' the mrc to get a box around the pdb to have lighter mrc files.
We then optionally filter the values far away from the PDB
Finally we need to resample the experimental maps to get a fixed voxel_size value of 1.
"""


def carve(mrc, pdb_name, out_name='carved.mrc', padding=4, filter_cutoff=-1):
    """
        This goes from full size to a reduced size, centered around a pdb.
        The main steps are :
            - Getting the coordinates to get a box around the PDB
            - Selecting the right voxels in this box
            - Updating the origin and the size of the 'cell'
            - Optionally filter out the values further away from filter_cutoff
    :param mrc: Either name or mrc file
    :param pdb_name: path to the pdb
    :param out_name: path to the output mrc
    :param padding: does not need to be an integer
    :param filter_cutoff: negative value will skip the filtering step. Otherwise it's a cutoff in Angstroms
    """
    # Get the bounds from the pdb
    coords = pymol_parse(pdbname=pdb_name)
    xyz_min = coords.min(axis=0)
    xyz_max = coords.max(axis=0)

    # Load the mrc data
    # the data and the 'x,y,z' annotation might not match.
    # axis_mapping tells us which 'xyz' dimension match which axis of the data : {first_axis : X, Y or Z}
    mrc = load_mrc(mrc, mode='r+')
    axis_mapping = {0: int(mrc.header.maps) - 1,
                    1: int(mrc.header.mapr) - 1,
                    2: int(mrc.header.mapc) - 1}
    reverse_axis_mapping = {value: key for key, value in axis_mapping.items()}
    voxel_size = mrc.voxel_size[..., None].view(dtype=np.float32)
    voxel_size = np.array(voxel_size)
    origin = mrc.header.origin[..., None].view(dtype=np.float32)
    origin = np.array(origin)
    # The shift array convention is as crappy as the s,r,c order.
    shift_array = np.array((mrc.header.nzstart,
                            mrc.header.nystart,
                            mrc.header.nxstart))
    shift_array_xyz = np.array([shift_array[reverse_axis_mapping[i]] for i in range(3)])
    data_shape_xyz = np.array([mrc.data.shape[reverse_axis_mapping[i]] for i in range(3)])

    # Using the origin, find the corresponding cells in the mrc
    mins_array = ((xyz_min - origin) / voxel_size - shift_array_xyz).astype(np.int, casting='unsafe')
    maxs_array = ((xyz_max - origin) / voxel_size - shift_array_xyz).astype(np.int, casting='unsafe')
    mins_array = np.max((mins_array, np.zeros_like(mins_array)), axis=0)
    maxs_array = np.min((maxs_array, data_shape_xyz - 1), axis=0)
    x_min, y_min, z_min = mins_array
    x_max, y_max, z_max = maxs_array
    grouped_bounds = [(i, j) for i, j in zip(mins_array, maxs_array)]
    shifted_origin = origin + (mins_array + shift_array_xyz - (padding,) * 3) * voxel_size

    # Extract those cells, one must be careful because x is not always the columns index
    data = mrc.data.copy()
    for array_axis, data_axis in axis_mapping.items():
        data = np.take(data, axis=array_axis, indices=range(*grouped_bounds[data_axis]))

    data = np.pad(data, pad_width=padding)

    # Optionnaly select only the cells that are at a certain distance to the pdbs
    if filter_cutoff > 0:
        filter_array = np.zeros_like(data)
        for coord in coords:
            new_coord = ((coord - xyz_min) / voxel_size + (padding,) * 3).astype(np.int, casting='unsafe')
            new_coord_axis = tuple(new_coord[axis_mapping[i]] for i in range(3))
            filter_array[new_coord_axis] += 1

        filter_array = np.float_(1. - (filter_array > 0.))
        filter_array = scipy.ndimage.distance_transform_edt(filter_array, sampling=voxel_size)
        filter_array = np.array([filter_array < filter_cutoff]).astype(np.float32)
        filter_array = np.reshape(filter_array, data.shape)
        data = filter_array * data

    # Update meta-data
    with mrcfile.new(out_name) as mrc:
        mrc.header.cella.x = (x_max - x_min + 2 * padding) * voxel_size[0]
        mrc.header.cella.y = (y_max - y_min + 2 * padding) * voxel_size[1]
        mrc.header.cella.z = (z_max - z_min + 2 * padding) * voxel_size[2]
        mrc.header.origin.x = shifted_origin[0]
        mrc.header.origin.y = shifted_origin[1]
        mrc.header.origin.z = shifted_origin[2]
        mrc.set_data(data)
        mrc.update_header_from_data()
        mrc.update_header_stats()


def subsample(mrc, padding=0, out_name='subsample.mrc'):
    """
        A script to change the voxel size of a mrc to 1A.
        The main operation is building a linear interpolation model and doing inference over it.
    :param mrc: Either string that will be attempted to be loaded or an MRC object
    :param out_name: Name of the mrc to dump
    :return: None
    """
    mrc = load_mrc(mrc, mode='r+')
    voxel_size = mrc.voxel_size[..., None].view(dtype=np.float32)
    voxel_size = np.array(voxel_size)
    axis_mapping = {0: int(mrc.header.maps) - 1,
                    1: int(mrc.header.mapr) - 1,
                    2: int(mrc.header.mapc) - 1}
    data = mrc.data.copy()
    data_axes = tuple(np.arange(0, data.shape[i]) * voxel_size[axis_mapping[i]] for i in range(3))
    interpolator = scipy.interpolate.RegularGridInterpolator(data_axes,
                                                             data,
                                                             method='linear',
                                                             bounds_error=False,
                                                             fill_value=0)
    cella = [mrc.header.cella.x, mrc.header.cella.y, mrc.header.cella.z]
    new_axes = tuple(np.arange(0, int(cella[axis_mapping[i]])) for i in range(3))
    # new_axes = tuple(np.arange(0, data.shape[i] * voxel_size[axis_mapping[i]]) for i in range(3))
    x, y, z = np.meshgrid(*new_axes, indexing='ij')
    flat = x.flatten(), y.flatten(), z.flatten()
    new_grid = np.vstack(flat).T
    new_data_grid = interpolator(new_grid).reshape(x.shape).astype(np.float32)
    new_data_grid = np.pad(new_data_grid, pad_width=padding)

    with mrcfile.new(out_name) as mrc2:
        mrc2.header.cella.x = int(mrc.header.cella.x) + 2 * padding
        mrc2.header.cella.y = int(mrc.header.cella.y) + 2 * padding
        mrc2.header.cella.z = int(mrc.header.cella.z) + 2 * padding
        mrc2.header.origin.x = mrc.header.origin.x - padding
        mrc2.header.origin.y = mrc.header.origin.y - padding
        mrc2.header.origin.z = mrc.header.origin.z - padding
        mrc2.set_data(new_data_grid)
        mrc2.update_header_from_data()
        mrc2.update_header_stats()


class MRC_grid():
    """
    This class is just used to factor out the interfacing with MRC files easier
    We loose the utilities developped in the mrcfile package and assume a unit voxel_size for simplicity.
    We always have consistent x,y,z and easy to access class member
    The convention for the axis is not the dominant one for MRC, we set it to be (X, Y, Z)
    so that the shape and origin are aligned
    """

    def __init__(self, MRC_file):
        self.mrc_obj = load_mrc(MRC_file)

        self.voxel_size = np.array(self.mrc_obj.voxel_size[..., None].view(dtype=np.float32))
        assert np.allclose(self.voxel_size, np.ones_like(self.voxel_size), atol=0.01)
        self.origin = np.array(self.mrc_obj.header.origin[..., None].view(dtype=np.float32))
        axis_mapping = (int(self.mrc_obj.header.maps) - 1,
                        int(self.mrc_obj.header.mapr) - 1,
                        int(self.mrc_obj.header.mapc) - 1)

        data = self.mrc_obj.data.copy()
        self.data = np.transpose(data, axes=axis_mapping)


def save_coords(coords, topology, outfilename, selection=None):
    """
    Save the coordinates to a pdb file
    • coords: coordinates
    • topology: topology
    • outfilename: name of the oupyt pdb
    • selection: Boolean array to select atoms
    """
    object_name = 'struct_save_coords'
    cmd.delete(object_name)
    if selection is None:
        selection = np.ones(len(topology['resids']), dtype=bool)
    for i, coords_ in enumerate(coords):
        if selection[i]:
            name = topology['names'][i]
            resn = topology['resnames'][i]
            resi = topology['resids'][i]
            chain = topology['chains'][i]
            elem = name[0]
            cmd.pseudoatom(object_name,
                           name=name,
                           resn=resn,
                           resi=resi,
                           chain=chain,
                           elem=elem,
                           hetatm=0,
                           segi=chain,
                           pos=list(coords_))
    cmd.save(outfilename, selection=object_name)
    cmd.delete(object_name)


def save_density(density, outfilename, origin, spacing=1, padding=0):
    """
    Save the density file as mrc for the given atomname
    """
    density = density.astype('float32')
    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(density.T)
        mrc.voxel_size = spacing
        mrc.header['origin']['x'] = origin[0] - padding + .5 * spacing
        mrc.header['origin']['y'] = origin[1] - padding + .5 * spacing
        mrc.header['origin']['z'] = origin[2] - padding + .5 * spacing
        mrc.update_header_from_data()
        mrc.update_header_stats()


if __name__ == '__main__':
    # import time

    pdb_name = "../data/4ci0/4ci0.cif"
    mrc_name = "../data/4ci0/emd_2513.mrc"
    carved_name = "../data/4ci0/emd_2513_carved.mrc"
    subsampled_name = "../data/4ci0/emd_2513_subsampled_2.mrc"

    # pdb_name = "data/5a33/5a33.cif"
    # mrc_name = "data/5a33/emd_3014.mrc"
    # carved_name = "data/5a33/emd_3014_carved.mrc"
    # subsampled_name = "data/5a33/emd_3014_subsampled.mrc"

    # carve(mrc=mrc_name, pdb_name=pdb_name, out_name=carved_name, filter_cutoff=6)
    subsample(mrc=carved_name, out_name=subsampled_name, padding=0)
