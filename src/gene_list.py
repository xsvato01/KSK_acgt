import pandas as pd
import h5py
import argparse
import numpy as np
HDF5_USE_FILE_LOCKING = 'FALSE'

#path = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_ksk.csv'
# bedpath = '/mnt/shared/MedGen/ACGTdatabase/data/hdf5/genes_noscaffold.txt'
#result_filename = "/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/result.csv"
#hdf_path = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_10k/'
#hdf_path = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_KSK/chrY_toAPP.hdf5'
annot_desc = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE"


def DeDups(atr, key, num_samples):
    altColumn = atr["alt"].reshape(-1, atr["alt"].shape[-1])
    gtColumn = atr["gt"].reshape(-1, atr["gt"].shape[-1])
    pdAtr = pd.concat([
        pd.DataFrame(altColumn),
        pd.DataFrame(gtColumn),
        pd.DataFrame(atr["samples"]),
        pd.DataFrame(atr["sex"]),
        pd.DataFrame(atr["region"]),
        pd.DataFrame(atr["annot"]),
    ], axis=1)
    pdAtr.columns = ["Alt1", "Alt2", "Alt3", "Gt1",
                     "Gt2", "Sample", "Sex", "Reg", "Annot"]
    pdAtr['is_homo'] = pdAtr.apply(
        lambda x: 1 if x.iloc[3] == x.iloc[4] else 0, axis=1)
    pdAtrMelt = pdAtr.melt(id_vars=[
                           'Alt1', 'Alt2', 'Alt3', 'Sample', 'Sex', 'Reg', 'Annot', 'is_homo'], value_vars=['Gt1', 'Gt2'])  # make genotypes in one column
    # not interested in reference alt
    pdAtrMelt_filt = pdAtrMelt[pdAtrMelt["value"] > 0].reset_index()
    pdAtrMelt_filt["Alt"] = pdAtrMelt_filt.apply(  # Genotype number is taken as an index to select correct Alt
        lambda x: x.iloc[x["value"]], axis=1)

    data = []
    for _, row in pdAtrMelt_filt.iterrows():
        data.append({
            "Chr": key.split(':')[0],
            "Pos": key.split(':')[1],
            "Ref": atr['ref'][0],
            "Alt": row['Alt'],
            "Sample": row["Sample"],
            "Sex": row["Sex"],
            "Reg": row["Reg"],
            "Homo_count": row["is_homo"],
            "Annot": row["Annot"],
        })

    data = pd.DataFrame(data)
    data_grouped = data.groupby(by=["Alt"]).agg({  # make Alt key and reorder by that
        "Chr": lambda x: x.iloc[0],
        "Pos": lambda x: x.iloc[0],
        "Ref": lambda x: x.iloc[0],
        "Sample": list,
        "Sex": list,
        "Reg": list,
        "Homo_count": np.sum,
        "Annot": lambda x: x.iloc[0]

    }).reset_index()
    data_grouped["Count"] = data_grouped["Sex"].str.len()
    data_grouped["ACGT_freq"] = data_grouped["Count"] / (2 * int(num_samples))
    data_grouped.insert(3, 'Alt', data_grouped.pop('Alt')
                        )  # reorder Alt after Ref
    return (data_grouped)


def main(hdf5_dir, hdf5_filenames, genesByChrom, Chr):
    result = pd.DataFrame()
    hdf_path = f"{hdf5_dir}chr{Chr}_{hdf5_filenames}.hdf5"
    print("opening: " + hdf_path)
    f = h5py.File(hdf_path, 'r')
    num_samples = f.attrs['Total_sample']
    chr_genes = genesByChrom[genesByChrom["Chromosome"]
                             == Chr]["genes"].values[0]
    print("Genes on Chr" + Chr + ": " + ','.join(chr_genes))
    print("Loading keys for Chr"+Chr)
    fkeys = f.keys()
    poc = 0
    fkeyslen = len(fkeys)
    for key in fkeys:
        poc += 1
        gene = f"{(f[key].attrs['annot'][0][0])}".split('|')[3]
        if gene in chr_genes:
            print(f"Chr{Chr} {100*poc/fkeyslen} %")
            out = DeDups(f[key].attrs, key, num_samples)
            result = pd.concat([result, out])

    if not len(result):
        return (pd.DataFrame())

    annotation = result.Annot.str.split("|", expand=True)
    annotation.columns = annot_desc.split('|')
    result = pd.concat([result, annotation], axis=1)
    result = result.drop(columns='Annot')
    f.close()
    return (result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create variant table out of HDF files")
    parser.add_argument(
        '--save_dir', help='Path to save result', default='/mnt/shared/MedGen/ACGTdatabase/data/hdf5/', required=False)
    parser.add_argument(
        '--save_filename', help='Result filename', default="KSK_results", required=False)
    parser.add_argument(
        '--hdf5_dir', help='Path to HDF directory', default='/mnt/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_KSK/', required=False)
    parser.add_argument(
        '--hdf5_filenames', help='One HDF per chromosome, files must start with ../chr[NUMBER]_xxx.hdf. Enter the xxx variable ending of your files', default='toAPP', required=False)
    parser.add_argument(
        '--gene_list', help='Full path to csv file containing csv with genes to analyze in the first column', default='/mnt/shared/MedGen/ACGTdatabase/data/hdf5/genes_ksk.csv', required=False)
    parser.add_argument(
        '--bed_genome', help='Enter path of annotated bed genome', default='', required=True)
    parser.add_argument(
        '--chr', help='Enter chromosome to to run analysis on', default='22', required=False)
    # parser.add_argument(
    #     '--num_samples', help='Enter number of total samples', default='673', required=False)
    args = parser.parse_args()

    print(args)
    # args = pd.DataFrame()
    # args.bed_genome = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_noscaffold.txt'
    # args.gene_list = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_ksk.csv'
    # args.hdf5_filenames = 'toAPP'
    # args.hdf5_dir = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_10k/'
    # args.save_dir = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/'
    # args.chr = '10'
    # #args.num_samples = 48

    if not (isinstance(args.chr, str)):
        args.chr = str(args.chr)

    gene_list = pd.read_csv(args.gene_list)
    bed_list = pd.read_csv(args.bed_genome, sep="\t")
    merged = gene_list.merge(bed_list, left_on='genes',
                             right_on='Name', how='left').sort_values(by=["Chromosome", "Start"])
    genesByChrom = merged.groupby(by=["Chromosome"]).agg(
        {"genes": list}).reset_index()

    if not args.chr in list(genesByChrom["Chromosome"]):
        print(f"Chromosmome {args.chr} is not in list of analysed genes")
        exit()

    result = main(hdf5_dir=args.hdf5_dir, hdf5_filenames=args.hdf5_filenames,
                  genesByChrom=genesByChrom, Chr=args.chr)
    if not result.empty:
        result.to_csv(f"{args.chr}.tsv", index=False, sep="\t")
