import pandas as pd
gene_list = pd.read_csv(
    '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_ksk.csv')
bed_list = pd.read_csv(
    "/Volumes/lamb/shared/MedGen/ACGT/nf_KSK_gene_panel/project/xsvato01/KSK_ACGT/utils/genes_noscaffold.txt", sep="\t")
merged = gene_list.merge(bed_list, left_on='genes',
                         right_on='Name', how='left').sort_values(by=["Chromosome", "Start"])
genesByChrom = merged.groupby(by=["Chromosome"]).agg(
    {"genes": list}).reset_index()
genesByChrom.to_csv("GenesPerChr.csv")
