"""
This script takes as input a pdb and an mrc and some extra selection tools and outputs a grid aligned with the mrc.
It also introduces the 'Complex' class that is fulling the Database object
"""

import os
import sys

from scipy import ndimage
import pymol.cmd as cmd
import numpy as np

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '..'))

from data_processing import mrc_utils


def get_coords_lig(pdb_name, selection=None):
    """
    Selects the CA and put them in a grid

    :param pdb_name: PDB file name for the protein
    :param selection:
    :return:
    """
    a = cmd.get_object_list('all')
    safe_sel = ['None', 'None'] + list({'prot'}.intersection(set(a)))
    safe_sel = ' or '.join(safe_sel)
    cmd.delete(safe_sel)
    cmd.load(pdb_name, 'prot')
    if selection is None:
        selection = 'prot'
    coords = cmd.get_coords(selection=f'{selection} and prot and name CA')
    return coords


def get_grid(coords, mrc, blur=True):
    """
    Generate a grid without channels from the coordinates
    """

    def gaussian_blur(data, sigma=1.):
        """
        Apply a gaussian blur to a grid object
        """
        if data.sum() > 0:
            data = np.float_(1. - (data > 0.))
            data = ndimage.distance_transform_edt(data)
            return np.exp(-data ** 2 / (2. * sigma ** 2))
        else:
            return data

    ranger = list(zip(mrc.origin, mrc.origin + mrc.data.shape * mrc.voxel_size))
    out, _ = np.histogramdd(coords, range=ranger, bins=mrc.data.shape)

    if blur:
        out = gaussian_blur(out)
    return out


# def get_grid_channels(coords, spacing, padding, xyz_min, xyz_max,
#                       grid_shape=None, blur=True):
#     """
#     Compute the 3D grids per channel from the coordinates
#     - coords: coordinates in the format [x, y, z, channel_id]
#     Remarks: - channel_id is zero-based numbering of the channels
#              - if channel_id is -1 only one channel is present (HETATM ligand
#                and not a protein partner)
#     """
#
#     def get_grid_shape(xyz_min, xyz_max, spacing, padding):
#         xm, ym, zm = xyz_min - (padding,) * 3
#         xM, yM, zM = xyz_max + (padding,) * 3
#         X = np.arange(xm, xM, spacing)
#         Y = np.arange(ym, yM, spacing)
#         Z = np.arange(zm, zM, spacing)
#         nx, ny, nz = len(X) - 1, len(Y) - 1, len(Z) - 1
#         return nx, ny, nz
#
#     channel_ids = coords[:, -1]
#     channel_ids_u = np.unique(channel_ids)
#     if grid_shape is None:
#         grid_shape = get_grid_shape(xyz_min, xyz_max, spacing, padding)
#     if len(channel_ids_u) == 1 and channel_ids_u[0] == -1:
#         # Its an organic HETATM ligand (from PL)
#         grid = get_grid(coords[:, :3], spacing, padding,
#                         blur=blur, xyz_max=xyz_max, xyz_min=xyz_min)
#         return grid[None, ...]
#     else:
#         grid = []
#         for channel_id, _ in enumerate(ATOMTYPES):
#             sel = (channel_ids == channel_id)
#             if sel.sum() > 0:
#                 coords_channel = coords[sel]
#                 grid_ = get_grid(coords_channel[:, :3], spacing, padding,
#                                  xyz_min=xyz_min, xyz_max=xyz_max,
#                                  blur=blur)
#                 grid.append(grid_)
#             else:  # emtpy channel
#                 grid.append(np.zeros(grid_shape))
#         grid = np.asarray(grid)
#         return grid


class Complex(object):
    """
    Object containing a protein-ligand system
    The main difficulty arises from the creation of the grid for the output,
    because we need those to align with the input mrc
    """

    def __init__(self, mrc, pdb_name, selection=None):
        """
        Code for this class is inspired by the Density class (PDB -> n,4)
        and the Complex class (n,4 + MRC -> Grid)
        """

        self.mrc = mrc_utils.MRC_grid(mrc)
        self.pdb_name = pdb_name
        self.n4 = get_coords_lig(pdb_name=pdb_name, selection=selection)

        self.out_grid = get_grid(coords=self.n4, mrc=self.mrc, blur=True)
        # self.save_mrc_lig()
        self.mrc.data
        self.out_grid
        pass

    def save_mrc_lig(self):
        """
        Save all the channels of the ligand in separate mrc files
        """
        outbasename = os.path.dirname(self.pdb_name)
        mrc_utils.save_density(density=self.out_grid,
                               outfilename=os.path.join(outbasename, 'out_grid.mrc'),
                               origin=self.mrc.origin)


if __name__ == '__main__':
    pdb_name = "../data/4ci0/4ci0.cif"
    # mrc_name = "../data/4ci0/emd_2513_subsampled.mrc"
    mrc_name = "../data/4ci0/simulated_cut.mrc"
    # pdb_name = "../data/phenix/5lzp_4128/5lzp.pdb"
    # mrc_name = "../data/phenix/5lzp_4128/4128_subsampled.mrc"
    # mrc_name = "../data/phenix/5lzp_4128/4128_simulated.mrc"
    # mrc_name = "../data/4ci0/simulated_cut.mrc"

    # Both were checked on this example and worked fine

    cp = Complex(mrc=mrc_name, pdb_name=pdb_name)
    print("mrc shape : ", cp.mrc.data.shape)
    print("lig shape : ", cp.out_grid.shape)
