
process CRAWL_ACGT {
	tag "Getting KSK on $Chr using $task.cpus CPUs and $task.memory memory"
	//publishDir  "${launchDir}/fastQC", mode:'copy'
	input:
	val Chr

 output:
 path '*.csv'

	script:
	"""
	touch ${Chr}.csv
 python $params.get_csvs --hdf5_dir /mnt/shared/MedGen/ACGTdatabase/data/hdf5/all_chr_10k/ --chr ${Chr}
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
"""
}


workflow {
	chromosome_ch = Channel.of(1..22, 'X', 'Y')
 files = CRAWL_ACGT(chromosome_ch)
	files.collect().view()
	MERGE_FILES(files.collect().view())
}
