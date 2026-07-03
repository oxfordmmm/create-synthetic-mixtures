// mix_pipeline/modules/mix.nf — standalone synthetic read mixing
// Lifted from the main pipeline's modules/create_mixes/mix.nf + modules/strain_detection/mapping.nf,
// trimmed to just what's needed to go from (challenge reads, background reads) -> mixed BAMs.

def pStr(p) { String.format('%.6f', p as double).replaceAll(/0+$/, '0').replaceAll(/\.$/, '.0') }


process filter_reads {
    // Keep only reads that map to the given reference genome.

    label 'ont_map'
    tag { name }

    storeDir "${params.outdir}/filtered_reads/${name}"

    input:
        tuple val(name), path(reads), path(fasta)

    output:
        tuple val(name), path("${name}_filtered.fastq.gz"), emit: reads

    script:
        """
        minimap2 -ax map-ont --secondary=no -t ${task.cpus} ${fasta} ${reads} \\
            | samtools view -F 4 -u \\
            | samtools fastq - \\
            | gzip > ${name}_filtered.fastq.gz
        echo "Kept \$(zcat ${name}_filtered.fastq.gz | awk 'NR%4==1' | wc -l) reads mapping to ${name}"
        """
}


process index_reference {

    label 'mapping'
    tag { reference_name }

    storeDir "${params.outdir}/index/${reference_name}"

    input:
        tuple val(reference_name), path(fasta)

    output:
        path("${reference_name}.mmi"), emit: indexed

    script:
        """
        minimap2 -d ${reference_name}.mmi -x map-ont ${fasta}
        """
}


// Filename stem carrying both strains and both proportions, e.g.
// barcode03_p0.1_barcode01_p0.9_rep1
def mixStem(chal_name, proportion, bg_name, replicate) {
    "${chal_name}_p${pStr(proportion)}_${bg_name}_p${pStr(1 - (proportion as double))}_rep${replicate}"
}


process create_mix_reads {
    // Subsample challenge and background reads at target proportion and concatenate.

    label 'mix_reads'
    tag { "${chal_name}_vs_${bg_name} p=${proportion} rep=${replicate}" }

    storeDir "${params.outdir}/mixes/${chal_name}_vs_${bg_name}/reads"

    input:
        tuple val(exp_label), val(proportion), val(replicate), val(seed),
              path(chal_reads), path(bg_reads), val(chal_name), val(bg_name)

    output:
        tuple val(exp_label), val(proportion), val(replicate), val(seed),
              path("${mixStem(chal_name, proportion, bg_name, replicate)}.fastq.gz"),
              val(chal_name), val(bg_name), emit: reads

    script:
        def total_reads = params.mix_total_reads as int
        def chal_n      = (proportion * total_reads).round() as int
        def bg_n        = total_reads - chal_n
        def stem        = mixStem(chal_name, proportion, bg_name, replicate)
        """
        if [ ${chal_n} -gt 0 ]; then
            seqtk sample -s${seed} ${chal_reads} ${chal_n} > mix_chal.fastq
        else
            touch mix_chal.fastq
        fi
        if [ ${bg_n} -gt 0 ]; then
            seqtk sample -s${seed} ${bg_reads} ${bg_n} > mix_bg.fastq
        else
            touch mix_bg.fastq
        fi
        cat mix_chal.fastq mix_bg.fastq | pigz -p ${task.cpus} > ${stem}.fastq.gz
        echo "Mix ${stem}: chal=${chal_n} bg=${bg_n}"
        """
}


process map_mix_reads {
    // Map mixed reads to challenge reference. All reads retained (including unmapped).

    label 'ont_map'
    tag { "${chal_name}_vs_${bg_name} p=${proportion} rep=${replicate}" }

    storeDir "${params.outdir}/mixes/${chal_name}_vs_${bg_name}/bam"

    input:
        tuple val(exp_label), val(proportion), val(replicate), val(seed),
              path(mix_reads), val(chal_name), val(bg_name), path(index)

    output:
        tuple val(exp_label), val(proportion), val(replicate), val(seed),
              path("${mixStem(chal_name, proportion, bg_name, replicate)}.bam"),
              path("${mixStem(chal_name, proportion, bg_name, replicate)}.bam.bai"), emit: bam

    script:
        def stem = mixStem(chal_name, proportion, bg_name, replicate)
        """
        minimap2 -ax map-ont --cs -t ${task.cpus} ${index} ${mix_reads} \\
            | samtools sort -@ 4 -o ${stem}.bam
        samtools index ${stem}.bam
        TOTAL=\$(samtools view -c ${stem}.bam)
        MAPPED=\$(samtools view -c -F 4 ${stem}.bam)
        echo "Mix ${stem}: total=\${TOTAL} mapped=\${MAPPED} unmapped=\$((\${TOTAL}-\${MAPPED}))"
        """
}
