process SALMON_INDEX {
    tag "${transcripts.name}"
    label 'process_medium'

    conda 'bioconda::salmon=1.10.3'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'

    input:
    path transcripts
    path genome
    val  use_decoys
    val  gencode

    output:
    path 'salmon_index', emit: index
    tuple val("${task.process}"), val('salmon'), eval('salmon --version | sed "s/salmon //"'), topic: versions

    script:
    def args = task.ext.args ?: ''
    def gc   = gencode ? '--gencode' : ''
    if (use_decoys)
        """
        # decoy-aware index: the whole genome is appended as decoy sequence so that reads
        # originating from unannotated loci are not forced onto similar transcripts
        grep '^>' ${genome} | cut -d ' ' -f 1 | sed 's/^>//' > decoys.txt
        cat ${transcripts} ${genome} > gentrome.fa
        salmon index --threads ${task.cpus} -t gentrome.fa -d decoys.txt ${gc} ${args} -i salmon_index
        """
    else
        """
        salmon index --threads ${task.cpus} -t ${transcripts} ${gc} ${args} -i salmon_index
        """
}
