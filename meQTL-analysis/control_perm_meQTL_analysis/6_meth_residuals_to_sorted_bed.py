import pandas as pd
import qtl
from qtl import io
import os

# Function to sort the bed file. It used to be integrated in an older version of the qtl package.
def sort_bed(bed_df, inplace=True):
    """Sort BED DataFrame"""
    sorted_df = bed_df.sort_values(['chr', 'start', 'end'], key=lambda x:
                    x.str.replace('chr','').str.replace('X','23').astype(int) if x.dtype == object else x,
                    inplace=inplace)
    if inplace:
        bed_df.reset_index(drop=True, inplace=True)
    else:
        sorted_df.reset_index(drop=True, inplace=True)
        return sorted_df

# Now run analysis in a loop
n = int(os.getenv("n"))
num = n + 1

for i in range(1, num):
    file = 'methylation_residuals_3_v2_{}.csv'.format(i)
    
    df = pd.read_csv(file)
    
    sort_bed(df, inplace=True)
    df

    result = '/methylation_residuals_3_{}.bed.gz'.format(i)
    
    io.write_bed(df, result,header=True, float_format=None)
   
