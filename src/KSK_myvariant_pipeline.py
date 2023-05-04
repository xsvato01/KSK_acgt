import myvariant
from liftover import get_lifter
import omim
from omim.db import Manager, OMIM_DATA
import pandas as pd
import argparse
import os

manager = Manager(dbfile=omim.DEFAULT_DB)
mv = myvariant.MyVariantInfo()


def main(filename):
    chrom = filename
    csv = pd.read_csv(filename, sep="\t")
    converter = get_lifter('hg38', 'hg19')
    final_list = []
    for index, row in csv.iterrows():
        hg19_cors = converter[chrom][row.Pos]
        row['clinical_significance'] = []
        row['NM_HGVSc'] = list()
        HGVSc = ""
        var_dbs = False  # getvariant does not return False/empty object
        var_dbs = mv.getvariant(
            f'{hg19_cors[0][0]}:g.{hg19_cors[0][1]}{row.Ref}>{row.Alt}')

        if var_dbs and var_dbs.get('clinvar'):
            row['NM_HGVSc'] = [x for x in var_dbs['clinvar']['hgvs']["coding"]
                               if len(x) > 1]  # sometimes HVSc is parsed incorectly in DB
            row['NM_HGVSc'].extend(''.join([x for x in var_dbs['clinvar']['hgvs']["coding"] if len(
                x) == 1]))  # sometimes HVSc is parsed incorectly in DB

            if isinstance(var_dbs['clinvar']['rcv'], list):  # if one entry it is a dict
                row['clinical_significance'] = [x['clinical_significance']
                                                for x in var_dbs['clinvar']['rcv']]
                row['clinical_significance'] = ', '.join(
                    row['clinical_significance'])
            else:
                row['clinical_significance'] = var_dbs['clinvar']['rcv']['clinical_significance']

        if var_dbs and var_dbs.get('snpeff'):
            # append if already created form clinvar
            if isinstance(var_dbs['snpeff']['ann'], list):  # if one entry it is a dict
                HGVSc = [
                    f"{x['feature_id']}:{x['hgvs_c']}" for x in var_dbs['snpeff']['ann']]
                row['NM_HGVSc'].extend(HGVSc)

            else:
                HGVSc = f"{var_dbs['snpeff']['ann']['feature_id']}:{var_dbs['snpeff']['ann']['hgvs_c']}"
                row['NM_HGVSc'].append(HGVSc)
            row['NM_HGVSc'] = set(row.NM_HGVSc)

        if index == 0 or csv.loc[index-1].Gene != row.Gene:
            omim_vars = manager.query(
                OMIM_DATA, 'ensembl_gene_id', row.Gene).all()
            if omim_vars:
                omim_parsed = eval(omim_vars[0].geneMap)
                gene_phenotype = ', '.join(
                    [x['Phenotype'] for x in omim_parsed])
                gene_inheritance = ', '.join(
                    [x['Inheritance'] for x in omim_parsed])
            else:
                gene_phenotype = '?'
                gene_inheritance = '?'

        row['NM_HGVSc'] = ', '.join(row['NM_HGVSc'])
        row['gene_phenotype'] = gene_phenotype
        row['gene_inheritance'] = gene_inheritance
        final_list.append(row)

    pd_final_list = pd.DataFrame(final_list)
    pd_final_list.insert(4, 'ACGT_percentage',
                         pd_final_list.pop('ACGT_percentage'))
    pd_final_list.insert(5, 'NM_HGVSc', pd_final_list.pop('NM_HGVSc'))
    pd_final_list.insert(6, 'clinical_significance',
                         pd_final_list.pop('clinical_significance'))
    pd_final_list.insert(7, 'gene_phenotype',
                         pd_final_list.pop('gene_phenotype'))
    pd_final_list.insert(8, 'gene_inheritance',
                         pd_final_list.pop('gene_inheritance'))
    # pd_final_list.to_csv(f"{chrom}_appended.tsv", index=False, sep="\t")
    return (pd_final_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extend HDF with custom settings for KSK")
    parser.add_argument(
        '--filename', help='Input filename no extension aka Chrom', default='2', required=False)
    args = parser.parse_args()

    if os.stat(args.filename).st_size == 0:
        print("OS stat 0")
        exit()

    out_df = main(args.filename)
    out_df.to_csv(f"{args.filename}_appended.tsv", index=False, sep="\t")
