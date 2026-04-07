#!/usr/bin/env python
# coding: utf-8

# # tesnorQTL analysis for all chromosomes

# Import packages
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans, post, pgen
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas {pd.__version__}")
import os
import os.path

print(tensorqtl.__version__)

import sys
print(sys.version)

### Run the analysis for all chromosomes
n = int(os.getenv("n"))
num = n + 1

for i in range(1,num):

    list_of_df_trans = []
    full_table_trans = pd.DataFrame()
    
    list_of_df_cis = []
    full_table_cis = pd.DataFrame()
    
    for x in range(1,23):

        # Set paths
        plink_prefix_path = f'brca_chr{x}'
        methylation_bed = 'methylation_residuals_3_{}.bed.gz'.format(i)
        
        prefix ='/significant_meqtls_lists_1000_ancestry_control_cis_trans/{}_meQTL'.format(i)
        
        prefix_final = '/significant_meqtls_lists_1000_ancestry_control_cis_trans/{}'.format(i)
        
        # Load phenotypes (methylation residuals after controlling for covariates of interest)
        phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(methylation_bed)
        
        # Load vcf files
        pr = genotypeio.PlinkReader(f'/vcfs/brca_chr{x}')
        
        # Load genotypes and variants into dataframes
        genotype_df = pr.load_genotypes()
        variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

        ####### cis-meQTL: obtaining nominal p-values for all variant-phenotype pairs per chrosomome
        cis.map_nominal(genotype_df, variant_df,
                    phenotype_df.loc[phenotype_pos_df['chr'] == f'chr{x}'],
                    phenotype_pos_df.loc[phenotype_pos_df['chr'] == f'chr{x}'],
                    prefix, window=1000000) 
        
        # If cis variants were detected
        if os.path.exists(f'{prefix}.cis_qtl_pairs.chr{x}.parquet'):
            # Load results
            pairs_df = pd.read_parquet(f'{prefix}.cis_qtl_pairs.chr{x}.parquet')
            # Add significant cis results to common list (pval < 2e-5 & effect size >= 0.25)
            list_of_df_cis.append(pairs_df[(pairs_df['pval_nominal'] < 2e-5) & (abs(pairs_df['slope']) >= 0.25)])

        ####### trans-meQTLs #############
        trans_df = trans.map_trans(genotype_df, phenotype_df,
                               return_sparse=True, pval_threshold=2e-5, maf_threshold=0,
                                batch_size=20000) 

        # Remove cis-associations from list with trans results
        trans_df = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=1000000)

        # Add significant trans results to significant trans result list
        list_of_df_trans.append(trans_df[abs(trans_df['b']) >= 0.25]) #trans_df
        
        # Set directory to delete individual chromosome mapping stats per permutation
        dir_name = "significant_meqtls_lists_1000_ancestry_control_cis_trans/"
        
        # Now delete .parquet files to save space
        test = os.listdir(dir_name)

        for item in test:
            if item.endswith(".parquet"):
                os.remove(os.path.join(dir_name, item))
                
    # Save list with significant results across all chromosomes
    full_df_cis = pd.concat(list_of_df_cis,ignore_index=True)
    full_df_trans = pd.concat(list_of_df_trans,ignore_index=True)

    full_df_cis.to_csv(f'{prefix_final}_cis_sig_all_3.csv')
    full_df_trans.to_csv(f'{prefix_final}_trans_sig_all_3.csv')
    
    
