# Load basic modules
import os
import pandas as pd
import hail as hl
hl.init(default_reference='GRCh38')

## Load and process meSNP location file 
data = pd.read_csv("locus_allele_rs_id_meSNPs_autosomal.csv", sep="\t") # or "locus_allele_rs_id_meSNPs_autosomal_control_perm.csv"" for control perm meSNPs

ht = hl.Table.from_pandas(data)

ht = ht.annotate(locus=hl.parse_locus(ht.locus))
ht = ht.key_by('locus')

sorted_table = ht.order_by(ht.locus)
sorted_table = sorted_table.key_by('locus')

## Load gnomad VAT table
vat_table =hl.read_table("gnomad.genomes.v4.1.sites.ht") 

# Change VAT table key
vat2 = vat_table.key_by("locus")

## Keep only variants that exist in both table
vat_annotations = vat2.semi_join(sorted_table)
vat_annotations = vat_annotations.key_by('locus','alleles')

# Save global frequency and frequency entry info (globals for frequency)
freq=vat_annotations.freq.AF
globals_freq=vat_annotations.globals.freq_meta

## Now save results
freq.export("freq_gnomad_meSNPs.tsv")
globals_freq.export("freq_globals.tsv")
