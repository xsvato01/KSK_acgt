
process CRAWL_ACGT {
	tag "Getting KSK on $Chr using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/perChr_allelicFreq", mode:'copy'
	input:
	val Chr

 output:
 path '*.tsv'

	script:
	"""
	touch ${Chr}.tsv
 python $params.get_csvs --chr ${Chr} --bed_genome ${params.bedpath} --gene_list ${params.geny_ksk} --hdf5_dir /mnt/shared/MedGen/ACGTdatabase/data/hdf5_variants_673samp/
	"""
}

process CRAWL_ACGT {
	tag "Getting KSK on $Chr using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/perChr_allelicFreq", mode:'copy'
	input:
	val Chr

 output:
 path '*.tsv'

	script:
	"""
	touch ${Chr}_appended.tsv
 python $params.get_additional --input_filename ${Chr}
	"""
}

process MERGE_FILES {
	tag "Merging KSK files using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}", mode:'copy'
	input:
	val CSVs

 output:
 path '*.csv'

	script:
"""
 python $params.cat_csvs --save_filename KSK_resul --csvs $CSVs
	ls
"""
}


workflow {
	chromosome_ch = Channel.of(1..22, 'X', 'Y')
 files = CRAWL_ACGT(chromosome_ch)

 APPEND_INFO(files)

	//files.collect().view()
	//MERGE_FILES(files.collect().view())
}
