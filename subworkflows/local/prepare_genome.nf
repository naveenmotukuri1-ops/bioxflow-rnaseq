/*
 * Prepare reference files: decompress, extract transcripts, build indices.
 */

include { GUNZIP as GUNZIP_FASTA       } from '../../modules/local/gunzip'
include { GUNZIP as GUNZIP_GTF         } from '../../modules/local/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPTS } from '../../modules/local/gunzip'
include { GFFREAD                      } from '../../modules/local/gffread'
include { TX2GENE                      } from '../../modules/local/tx2gene'
include { SALMON_INDEX                 } from '../../modules/local/salmon_index'
include { STAR_GENOMEGENERATE          } from '../../modules/local/star_genomegenerate'

workflow PREPARE_GENOME {
    main:
    def fasta_in = file(params.fasta, checkIfExists: true)
    def gtf_in   = file(params.gtf, checkIfExists: true)

    ch_fasta = fasta_in.name.endsWith('.gz') ? GUNZIP_FASTA(fasta_in).gunzip : channel.value(fasta_in)
    ch_gtf   = gtf_in.name.endsWith('.gz')   ? GUNZIP_GTF(gtf_in).gunzip     : channel.value(gtf_in)

    TX2GENE(ch_gtf)

    if (params.salmon_index) {
        ch_salmon_index = channel.value(file(params.salmon_index, checkIfExists: true))
    }
    else {
        if (params.transcript_fasta) {
            def tx_in = file(params.transcript_fasta, checkIfExists: true)
            ch_transcripts = tx_in.name.endsWith('.gz') ? GUNZIP_TRANSCRIPTS(tx_in).gunzip : channel.value(tx_in)
        }
        else {
            ch_transcripts = GFFREAD(ch_fasta, ch_gtf).transcripts
        }
        ch_salmon_index = SALMON_INDEX(ch_transcripts, ch_fasta, params.salmon_decoys, params.gencode).index
    }

    ch_star_index = channel.empty()
    if (!params.skip_alignment) {
        ch_star_index = params.star_index
            ? channel.value(file(params.star_index, checkIfExists: true))
            : STAR_GENOMEGENERATE(ch_fasta, ch_gtf).index
    }

    emit:
    tx2gene      = TX2GENE.out.tsv
    salmon_index = ch_salmon_index
    star_index   = ch_star_index
}
