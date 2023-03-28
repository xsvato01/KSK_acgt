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


def DeDups(atr, key):
    atrList = atr["alt"].reshape(-1, atr["alt"].shape[-1]).tolist()
    atrUniques = np.unique(atrList, axis=0)
    data = []
    for uniq in atrUniques:
        boolmask = (uniq == atr["alt"]).all(axis=2).squeeze()
        data.append({
            "chr": key.split(':')[0],
            "pos": key.split(':')[1],
            "ref": atr['ref'][0],
            "alt": atr['alt'][boolmask][0].squeeze().tolist(),
            "gt": atr['gt'][boolmask].squeeze().tolist(),
            "annot": atr['annot'][boolmask][0].squeeze().tolist(),
            # "sex": atr['sex'][boolmask],
            "samples": atr['samples'][boolmask].squeeze().tolist(),
            "ACGT_percentage": len(atr['samples'][boolmask]) / int(args.num_samples) * 100,
            "region": atr['region'][boolmask].squeeze().tolist()
        })
    return (data)


def main(hdf5_dir, hdf5_filenames, genesByChrom, Chr):
    result = []
    poc = 0
    hdf_path = f"{hdf5_dir}chr{Chr}_{hdf5_filenames}.hdf5"
    print("opening: " + hdf_path)
    f = h5py.File(hdf_path, 'r')
    chr_genes = genesByChrom[genesByChrom["Chromosome"]
                             == Chr]["genes"].values[0]
    # for j in chrom_genes[chrom_genes["Chromosome"] == i]["genes"].values[0]:
    print("Genes on Chr" + Chr + ": " + ','.join(chr_genes))
    print("Loading keys for Chr"+Chr)
    fkeys = f.keys()
    poc = 0
    fkeyslen = len(fkeys)
    for key in fkeys:
        poc += 1
        gene = f"{(f[key].attrs['annot'][0][0])}".split('|')[3]
        # if np.unique(f[key].attrs["alt"].reshape(-1, f[key].attrs["alt"].shape[-1]).tolist(), axis=0).shape[0] > 1:
        if gene in chr_genes:
            print(f"Chr{Chr} {100*poc/fkeyslen} %")
            if gene in chr_genes:
                out = DeDups(f[key].attrs, key)
                result.append(out)
    if not len(result):
        return (pd.DataFrame())
    flat_list = pd.DataFrame([y for x in result for y in x])
    annotation = flat_list.annot.str.split("|", expand=True)
    annotation.columns = annot_desc.split('|')
    flat_list = pd.concat(
        [flat_list, annotation], axis=1)
    flat_list = flat_list.drop(columns='annot')
    f.close()
    return (flat_list)


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
    parser.add_argument(
        '--num_samples', help='Enter number of total samples', default='673', required=False)
    args = parser.parse_args()

    print(args)
    # args.bed_genome = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_noscaffold.txt'
    # args.gene_list = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_ksk.csv'
    # args.hdf5_dir = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_10k/'
    # args.save_dir = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/'
    # args.num_samples = 673

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

    flat_list = main(hdf5_dir=args.hdf5_dir, hdf5_filenames=args.hdf5_filenames,
                     genesByChrom=genesByChrom, Chr=args.chr)
    if not flat_list.empty:
        flat_list.to_csv(f"{args.chr}.csv")
