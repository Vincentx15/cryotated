import os
import sys

import numpy as np
import torch
from torch.utils.data.dataset import Dataset
from torch.utils.data import Subset, DataLoader

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '..'))

from data_processing.Complex import Complex

"""
This script takes an hdf5 and organize it into a Database object that yields Complexes
Then this Database is fed to Pytorch Dataloaders
"""


class EMDataset(Dataset):
    """
    Encapsulates the Complex building
    """

    def __init__(self, data_path, simulated=True):
        self.data_path = data_path

        # For pytorch loading, we need the file reading in the get_item
        self.keys = os.listdir(data_path)
        self.simulated = simulated

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, item):
        """
        Returns the desired complex.
        :param item:
        :return:
        """

        dirname = self.keys[item]
        # print(dirname)
        pdb_name, mrc = dirname.split("_")
        pdb_name = os.path.join(self.data_path, dirname, f'{pdb_name}.pdb')
        mrc_file = f'{mrc}_simulated.mrc' if self.simulated else f'{mrc}_subsampled.mrc'
        mrc_path = os.path.join(self.data_path, dirname, mrc_file)

        complex = Complex(pdb_name=pdb_name, mrc=mrc_path)
        in_grid = complex.mrc.data
        out_grid = complex.out_grid

        return torch.from_numpy(in_grid), torch.from_numpy(out_grid)


class Loader:
    def __init__(self,
                 data_path,
                 batch_size=1,
                 num_workers=10,
                 simulated=True):
        """

        :param data_path: hdf5 file to load
        :param batch_size:
        :param num_workers:
        :param rotate:
        """
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.simulated = simulated
        self.dataset = EMDataset(data_path=data_path, simulated=simulated)

    def get_data(self):
        n = len(self.dataset)

        np.random.seed(0)
        torch.manual_seed(0)

        split_train, split_valid = 0.7, 0.85
        train_index, valid_index = int(split_train * n), int(split_valid * n)
        indices = list(range(n))

        train_indices = indices[:train_index]
        valid_indices = indices[train_index:valid_index]
        test_indices = indices[valid_index:]

        train_set = Subset(self.dataset, train_indices)
        valid_set = Subset(self.dataset, valid_indices)
        test_set = Subset(self.dataset, test_indices)

        train_loader = DataLoader(dataset=train_set, shuffle=True, batch_size=self.batch_size,
                                  num_workers=self.num_workers, worker_init_fn=np.random.seed)
        valid_loader = DataLoader(dataset=valid_set, shuffle=True, batch_size=self.batch_size,
                                  num_workers=self.num_workers)
        test_loader = DataLoader(dataset=test_set, shuffle=True, batch_size=self.batch_size,
                                 num_workers=self.num_workers)

        return train_loader, valid_loader, test_loader


if __name__ == '__main__':
    loader = Loader(data_path="../data/phenix", num_workers=0)
    tl, _, _ = loader.get_data()
    print(f'{len(tl)} systems in our data')
    print()
    for mrc, lig in tl:
        pass
        print("mrc shape : ", mrc.shape)
        print("lig shape : ", lig.shape)
