run = "${params.datain}".split("/")
run = run[run.size()-1]
launchDir = "${launchDir}/${run}"


process CRAWL_ACGT {
	tag "Getting KSK on $name using $task.cpus CPUs and $task.memory memory"
	//publishDir  "${launchDir}/fastQC", mode:'copy'

	script:
	"""
python $params.py
	"""
}


workflow {
 CRAWL_ACGT()
}
