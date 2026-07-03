#!/usr/bin/env nextflow
// mix_pipeline/main.nf — standalone synthetic read-mix generator
//
// Takes one challenge strain (params.challenge: {name, fasta, reads}) and one or more
// background strains (params.backgrounds: [{name, fasta?, reads}, ...]) and produces
// mixed FASTQs/BAMs at a range of proportions. Extracted from the main pipeline's
// Stage C so mixes can be generated without running reference-panel/marker/detection
// stages.
//
// Each background is mixed against the challenge as its own independent pairwise
// experiment (challenge_vs_<background>) — this is N separate 2-strain mix sets,
// not one N-way combined mix.
//
// Usage: bash run.sh [-resume] [--param value ...]
nextflow.enable.dsl=2

include { filter_reads; index_reference; create_mix_reads; map_mix_reads } from './modules/mix'

workflow {

    if (!params.backgrounds || params.backgrounds.size() < 1) {
        error "params.backgrounds must contain at least one background strain."
    }

    def chal_fasta_file = file(params.challenge.fasta)
    def chal_name        = params.challenge.name ?: chal_fasta_file.baseName.tokenize('_')[0]
    def no_file           = file('NO_FILE')

    def ch_chal_fasta = channel.value(chal_fasta_file)
    def ch_chal_reads = channel.value(file(params.challenge.reads))

    // tuple(name, reads, fasta_or_NO_FILE) per background
    def ch_backgrounds = channel.fromList(params.backgrounds)
        .map { b -> tuple(b.name, file(b.reads), b.fasta ? file(b.fasta) : no_file) }

    def props = ([params.mix_proportions].flatten()
        .collect { v -> v.toString().trim().toDouble() })

    def combos = [props, (1..(params.mix_replicates as int))].combinations()
        .collect { p, r -> tuple(p as double, r as int, (r * 1000 + (p * 100).toInteger()) as int) }

    // ── Filter reads to their own reference (removes off-target contamination) ──
    def ch_chal_in = ch_chal_reads.map { r -> tuple(chal_name, r, chal_fasta_file) }
    def ch_bg_in   = ch_backgrounds.filter { _n, _r, f -> f.name != 'NO_FILE' }

    filter_reads(ch_chal_in.mix(ch_bg_in))

    def ch_chal_prepared = filter_reads.out.reads
        .filter { name, _r -> name == chal_name }
        .map    { _name, r -> r }

    // tuple(bg_name, prepared_reads) — filtered ones for backgrounds that had a fasta,
    // raw reads for the rest
    def ch_bg_filtered = filter_reads.out.reads
        .filter { name, _r -> name != chal_name }

    def ch_bg_unfiltered = ch_backgrounds
        .filter { _n, _r, f -> f.name == 'NO_FILE' }
        .map    { n, r, _f -> tuple(n, r) }

    def ch_bg_prepared = ch_bg_filtered.mix(ch_bg_unfiltered)

    // ── Index reference & build mixes ────────────────────────────────────────
    index_reference(ch_chal_fasta.map { fa -> tuple(chal_name, fa) })

    channel.fromList(combos)
        .combine(ch_bg_prepared)
        .combine(ch_chal_prepared)
        .map { p, rep, seed, bg_name, bg_r, chal_r ->
            tuple("${chal_name}_vs_${bg_name}", p, rep, seed, chal_r, bg_r, chal_name, bg_name) }
        .set { ch_mix_inputs }

    create_mix_reads(ch_mix_inputs)

    create_mix_reads.out.reads
        .combine(index_reference.out.indexed)
        .set { ch_map_inputs }

    map_mix_reads(ch_map_inputs)
}
