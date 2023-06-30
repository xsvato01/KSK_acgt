
process CRAWL_ACGT {
	tag "Getting KSK on $Chr using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/perChr_Allele_appended_info", mode:'copy'
	input:
	val Chr

 output:
 tuple val (Chr), path ("*.tsv")

	script:
	"""
	touch ${Chr}.tsv
 python $params.get_csvs --chr ${Chr} --bed_genome ${params.bedpath} --gene_list ${params.geny_ksk} --hdf5_dir /mnt/shared/MedGen/ACGTdatabase/data/hdf5_variants_673samp/
	"""
	// /mnt/shared/MedGen/ACGTdatabase/data/hdf5_variants_673samp/
}

process APPEND_INFO {
	tag "Getting KSK on $Name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/perChr_Allele_appended_info", mode:'copy'
	input:
 tuple val(Name), path(FilePath)

 output:
 path "*"

	script:
	"""
 python $params.get_additional --filename ${Name} --filepath ${FilePath}
	"""
}


//not used anymore:
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
