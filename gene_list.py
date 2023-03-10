import pandas as pd
import h5py
import argparse
import numpy as np

path = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_ksk.csv'
bedpath = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_noscaffold.txt'
result_filename = "/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/result.csv"
#hdf_path = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_10k/'
#hdf_path = '/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_KSK/chrY_toAPP.hdf5'


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
            "region": atr['region'][boolmask].squeeze().tolist()
        })
    return (data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create variant table out of HDF files")
    parser.add_argument(
        '--save_dir', help='Path to save result', default='/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/', required=False)
    parser.add_argument(
        '--save_filename', help='Result filename', default="results", required=False)
    parser.add_argument(
        '--hdf5_dir', help='Path to HDF directory', default='/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_10K/', required=False)
    parser.add_argument(
        '--hdf5_filenames', help='One HDF per chromosome, files must start with ../chr[NUMBER]_xxx.hdf. Enter the xxx variable ending of your files', default='toAPP', required=False)
    parser.add_argument(
        '--gene_list', help='Full path to csv file containing csv with genes to analyze in the first column', default='/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/genes_ksk.csv', required=False)
    args = parser.parse_args()

    print(args)

    gene_list = pd.read_csv(args.gene_list)
    bed_list = pd.read_csv(bedpath, sep="\t")
    merged = gene_list.merge(bed_list, left_on='genes',
                             right_on='Name', how='left').sort_values(by=["Chromosome", "Start"])
    # merged.to_csv('/Volumes/lamb/shared/MedGen/ACGTdatabase/data/hdf5/merged.csv')
    genesByChrom = merged.groupby(by=["Chromosome"]).agg(
        {"genes": list}).reset_index()
    result = []
    # for i in genesByChrom["Chromosome"]:
    for i in list('X'):
        poc = 0
        hdf_path = f"{args.hdf5_dir}chr{i}_{args.hdf5_filenames}.hdf5"
        f = h5py.File(hdf_path, 'r')
        chr_genes = genesByChrom[genesByChrom["Chromosome"]
                                 == i]["genes"].values[0]
        # for j in chrom_genes[chrom_genes["Chromosome"] == i]["genes"].values[0]:
        fkeys = f.keys()
        for key in fkeys:
            gene = f"{(f[key].attrs['annot'][0][0])}".split('|')[3]
            if np.unique(f[key].attrs["alt"].reshape(-1, f[key].attrs["alt"].shape[-1]).tolist(), axis=0).shape[0] > 1:
                poc = poc+1
                print(f"Chr{i} {poc}")
            if gene in chr_genes:
                out = DeDups(f[key].attrs, key)
                result.append(out)
    flat_list = pd.DataFrame([y for x in result for y in x])
    flat_list = pd.concat(
        [flat_list, flat_list.annot.str.split("|", expand=True)], axis=1)
    flat_list.to_csv(f"{args.save_dir}{args.save_filename}.csv")
