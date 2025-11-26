import torch
import pyarrow.parquet as pq

from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from pytorch_lightning import LightningDataModule


AA = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {a: i for i, a in enumerate(AA)}

def one_hot_encode(seq):
    idxs = torch.tensor([AA_TO_IDX[x] for x in seq], dtype=torch.long)
    return torch.nn.functional.one_hot(idxs, num_classes=len(AA)).float()


class RosettaDataModule(LightningDataModule):
    def __init__(self,
                 parquet_path=None,
                 batch_size=16,
                 set_key="set",
                 num_workers=1,
                 persistent_workers=True,
                 sample_from_groups=True,
                 group_sample_n=1000):
        super().__init__()
        self.parquet_path = parquet_path
        self.sample_from_groups = sample_from_groups
        self.group_sample_n = group_sample_n
        self.batch_size = batch_size
        self.set_key = set_key
        self.num_workers = num_workers
        self.persistent_workers = persistent_workers

    def setup(self, stage=None):
        if stage == 'fit' or stage is None:
            # read data
            self.df = pq.read_table(self.parquet_path).to_pandas()

            if self.sample_from_groups is not None:
                grouped = self.df.groupby("pdb_fn")
                # adapt for replace and n
                self.df = grouped.sample(n=self.group_sample_n, replace=True, random_state=42).reset_index(drop=True)

            # Masking
            train_mask = self.df[self.set_key] == "train"
            val_mask = self.df[self.set_key] == "val"
            test_mask = self.df[self.set_key] == "test"

            seqs = self.df['mutated Sequence'].tolist()
            rosetta_scores = torch.Tensor(self.df['total_Score'])  # [1,0.9,2, ...]

            self.train_dataset = RosettaDataset(seqs[train_mask],
                                                rosetta_scores[train_mask])
            self.val_dataset = RosettaDataset(seqs[val_mask],
                                              rosetta_scores[val_mask])
            self.test_dataset = RosettaDataset(seqs[test_mask],
                                               rosetta_scores[test_mask])
            print(
                f"Train/Val/Test split: {len(self.train_dataset)}/{len(self.val_dataset)}/{len(self.test_dataset)} samples")

    def train_dataloader(self):
        return DataLoader(self.train_dataset,
                          batch_size=self.batch_size,
                          shuffle=True,
                          num_workers=self.num_workers,
                          persistent_workers=self.persistent_workers)

    def val_dataloader(self):
        return DataLoader(self.val_dataset,
                          batch_size=self.batch_size,
                          shuffle=False,
                          num_workers=self.num_workers,
                          persistent_workers=self.persistent_workers)

    def test_dataloader(self):
        return DataLoader(self.test_dataset,
                          batch_size=self.batch_size,
                          shuffle=False,
                          num_workers=self.num_workers,
                          persistent_workers=self.persistent_workers)


class RosettaDataset(Dataset):
    def __init__(self, seqs, rosetta_scores):
        self.seqs = seqs
        self.rosetta_scores = rosetta_scores

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        seq = self.seqs[idx]
        x = one_hot_encode(seq)
        y = self.rosetta_scores[idx]
        return x, y