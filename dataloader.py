import torch
import torch.utils.data as data
import random
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
import anndata

from config import Config


def read_mtx(mtx_file_path: str):
    """
    read sparse matrix file, such as .mtx file

    :param mtx_File_path:str, the path of the file need to be read

    return compressed sparse row matrix: scipy.sparse.csr_matrix
    """
    sparse_matrix = mmread(mtx_file_path)
    expression_matrix = sparse_matrix.A
    # (cells x genes) or (cells x peaks)
    expression_matrix = expression_matrix.transpose((1, 0))
    csr_expression_matrix = csr_matrix(expression_matrix, dtype=np.float64)

    # return a compressed sparse row matrix
    return csr_expression_matrix

def initAnndata(mtx_path:str, features_path:str, barcodes_path:str):
    """
    Inititlize an anndata object from matrix.mtx, barcodes.tsv, and genes.tsv/peaks.tsv.

    return anndata obejct
    """
    matrix = read_mtx(mtx_path)
    features = pd.read_csv(features_path, sep='\t', header=None)
    barcodes = pd.read_csv(barcodes_path, sep='\t', header=None)
    # matrix = matrix.transpose()

    adata = anndata.AnnData(X=matrix, 
                            obs=pd.DataFrame(index=barcodes[0]), 
                            var=pd.DataFrame(index=features[0])
                            )
    return adata

class Dataloader(data.Dataset):
    def __init__(self, train = True, data_reader = None, labels = None, protein_reader = None):
        self.train = train        
        self.data_reader, self.labels, self.protein_reader = data_reader, labels, protein_reader
        self.input_size = self.data_reader.shape[1]
        self.sample_num = self.data_reader.shape[0]
        
        self.input_size_protein = None
        if protein_reader is not None:
            self.input_size_protein = self.protein_reader.shape[1]

    def __getitem__(self, index):
        if self.train:
            # get atac data            
            rand_idx = random.randint(0, self.sample_num - 1)
            # sample = np.array(self.data_reader[rand_idx].todense())
            sample = np.array(self.data_reader.todense())
            print(sample.shape)
            sample = sample.reshape((-1, self.input_size))
            print(sample.shape)
            # 大于0的数设置为1，其他的设置为0，并转换为float类型
            in_data = (sample>0).astype(np.float64)  # binarize data
            
            if self.input_size_protein is not None:
                sample_protein = np.array(self.protein_reader[rand_idx].todense())
                sample_protein = sample_protein.reshape((1, self.input_size_protein))
                in_data = np.concatenate((in_data, sample_protein), 1)
                
            in_label = self.labels[rand_idx]
 
            return in_data, in_label

        else:
            sample = np.array(self.data_reader[index].todense())
            sample = sample.reshape((1, self.input_size))
            in_data = (sample>0).astype(np.float64)  # binarize data

            if self.input_size_protein is not None:
                sample_protein = np.array(self.protein_reader[index].todense())
                sample_protein = sample_protein.reshape((1, self.input_size_protein))
                in_data = np.concatenate((in_data, sample_protein), 1)
                
            #in_data = in_data.reshape((1, self.input_size))
            in_label = self.labels[index]
 
            return in_data, in_label

    def __len__(self):
        return self.sample_num
    

class PrepareDataloader:
    def __init__(self, config:Config) -> None:
        self.config = config

        self.joint_profiles = {}
        # load RNA data
        self.joint_profiles['gene_expression'] = read_mtx(self.config.dataset_path_dict['gene_expression'])
        self.joint_profiles['gene_names'] = np.loadtxt(config.dataset_path_dict['gene_names'], dtype=str)
        self.joint_profiles['gene_barcodes'] = np.loadtxt(config.dataset_path_dict['gene_barcodes'], dtype=str)

        # load ATAC data
        self.joint_profiles['atac_expression'] = read_mtx(self.config.dataset_path_dict['atac_expression'])
        self.joint_profiles['atac_names'] = np.loadtxt(config.dataset_path_dict['atac_names'], dtype=str)
        self.joint_profiles['atac_barcodes'] = np.loadtxt(config.dataset_path_dict['atac_barcodes'], dtype=str)        
        
        pass
