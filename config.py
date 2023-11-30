import torch


class Config:
    def __init__(self) -> None:
        self.dataset = "SNARE_seq"
        self.use_cuda = True
        if not self.use_cuda:
            self.device = torch.device('cpu')
        else:
            self.device = torch.device('cuda')

        if self.dataset == "SNARE_seq":
            self.dataset_path_dict = {
                "gene_names": 'SNARE_seq/cDNA_genes.txt',
                "gene_expression": 'SNARE_seq/cDNA_sparse_matrix.mtx',
                "gene_barcodes": 'SNARE_seq/cDNA_barcodes.txt',
                "atac_names": 'SNARE_seq/ATAC_peaks.txt',
                "atac_expression": 'SNARE_seq/ATAC_sparse_matrix.mtx',
                "atac_barcodes": 'SNARE_seq/ATAC_barcodes.txt'
            }
            self.batch_size = 256


        