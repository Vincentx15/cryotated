import os
import torch
from torch.utils.data.dataset import Dataset
from torch.utils.data import Subset, DataLoader

import numpy as np

"""
This script takes an hdf5 and organize it into a Database object that yields Complexes
Then this Database is fed to Pytorch Dataloaders
"""

class EMDataset(Dataset):
    """
        Uses a HDF5 file as defined above and turn it into a Pytorch data set
        """

    def __init__(self, data_path, simulated=True):
        self.data_path = data_path

        # For pytorch loading, we need the file reading in the get_item
        self.keys = os.listdir(data_path)

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, item):
        """
        Returns the desired complex.
        :param item:
        :return:
        """

        dirname = self.keys[item]
        pdb_name, mrc = dirname.split("_")

        complex = 0
        cp = self.database.get_complex(prot, lig, rotate=self.rotate)

        gprot = cp.grid_prot.astype(np.float32)
        glig = cp.grid_lig.astype(np.float32)

        # cp.save_mrc_lig()
        # cp.save_mrc_prot()

        return torch.from_numpy(gprot), torch.from_numpy(glig), cp.is_PL


class InferenceDataset(Dataset):
    """
        Almost the same with less options and different returns, with the name of the ligand.
        """

    def __init__(self, data_file):
        self.data_file = data_file

        # For pytorch loading, we need the file reading in the get_item
        self.database = Database(data_file)
        self.keys = self.database.get_protlig_keys()
        self.database = None

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, item):
        """
        Returns the desired complex.
        :param item:
        :return:
        """

        if self.database is None:
            self.database = Database(self.data_file)

        prot, lig = self.keys[item]
        cp = self.database.get_complex(prot, lig, rotate=False)
        gprot = cp.grid_prot.astype(np.float32)

        return torch.from_numpy(gprot), prot, lig, cp.is_PL


class Loader:
    def __init__(self, df,
                 batch_size=1,
                 num_workers=10,
                 rotate=True):
        """

        :param df: hdf5 file to load
        :param batch_size:
        :param num_workers:
        :param rotate:
        """
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.dataset = HDPLDataset(data_file=df, rotate=rotate)

    def get_data(self):
        n = len(self.dataset)

        np.random.seed(0)
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
