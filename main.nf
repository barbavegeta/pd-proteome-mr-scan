#!/usr/bin/env nextflow
/*
 * Proteome-wide cis-MR + colocalisation scan: plasma proteins -> Parkinson's disease.
 * One RUN_PROTEIN task per protein, run in parallel; results are FDR-corrected,
 * tiered and rendered into an HTML report.
 */
nextflow.enable.dsl = 2

process PREPARE_OUTCOME {
    tag "${gwas.simpleName}"
    input:
    path gwas
    output:
    path "outcome_by_chr", emit: dir
    script:
    """
    prepare_outcome.py --gwas ${gwas} --format ${params.outcome_format} --outdir outcome_by_chr
    """
}

process RUN_PROTEIN {
    tag "${meta.gene}"
    input:
    tuple val(meta), path(exposure, stageAs: 'exposure/*')
    path outcome_dir
    path ld_files
    output:
    path "${meta.key}.summary.json", emit: summary
    path "${meta.key}.mr.tsv", emit: mr
    tuple path("${meta.key}.region.tsv.gz"), path("${meta.key}.instruments.tsv"), emit: region
    script:
    def ld = params.ld_ref ? "--ld-ref ${file(params.ld_ref).name}" : ""
    """
    run_protein.py \\
        --protein-id '${meta.protein_id}' --gene ${meta.gene} --chr ${meta.chr} \\
        --start ${meta.start} --end ${meta.end} \\
        --exposure ${exposure} --exposure-format ${params.exposure_format} \\
        --outcome-dir ${outcome_dir} ${ld} --prefix ${meta.key} \\
        --window-kb ${params.window_kb} --info-min ${params.info_min} --maf-min ${params.maf_min} \\
        --p-threshold ${params.p_threshold} --f-min ${params.f_min} \\
        --clump-r2 ${params.clump_r2} --clump-kb ${params.clump_kb} \\
        --p1 ${params.p1} --p2 ${params.p2} --p12 ${params.p12} \\
        --outcome-ncase ${params.outcome_ncase} --outcome-ncontrol ${params.outcome_ncontrol}
    """
}

process AGGREGATE {
    publishDir params.outdir, mode: 'copy'
    input:
    path summaries
    path mrs
    output:
    path "scan_results.tsv", emit: results
    path "mr_all_methods.tsv", optional: true
    script:
    """
    aggregate.py --summaries *.summary.json --mr *.mr.tsv --fdr ${params.fdr} --out scan_results.tsv
    """
}

process REPORT {
    publishDir params.outdir, mode: 'copy'
    input:
    path results
    path regions
    path params_json
    output:
    path "pd_proteome_mr_report.html"
    script:
    """
    make_report.py --results ${results} --region-dir . --highlight ${params.highlight} \\
        --params ${params_json} --out pd_proteome_mr_report.html
    """
}

workflow {
    if (!params.samplesheet || !params.outcome) {
        error "Provide --samplesheet and --outcome (see README)"
    }
    def sheet_dir = file(params.samplesheet).parent
    proteins = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .filter { row -> !params.genes || params.genes.tokenize(',').contains(row.gene) }
        .map { row ->
            // relative exposure paths are resolved against the samplesheet's folder
            def exp = row.exposure.startsWith('/') || row.exposure.contains('://') ? file(row.exposure) : sheet_dir.resolve(row.exposure)
            if (!exp.exists()) error "Exposure file not found: ${exp}"
            def key = row.protein_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
            tuple([protein_id: row.protein_id, key: key, gene: row.gene, chr: row.chr,
                   start: row.start, end: row.end], exp)
        }

    // LD reference: one genome-wide PLINK set, or per-chromosome sets via '{chr}' in the prefix
    ld_files = params.ld_ref
        ? Channel.fromPath("${params.ld_ref.replace('{chr}', '*')}.{bed,bim,fam}", checkIfExists: true).collect()
        : Channel.value([])

    outcome = PREPARE_OUTCOME(file(params.outcome, checkIfExists: true))
    RUN_PROTEIN(proteins, outcome.dir, ld_files)

    AGGREGATE(RUN_PROTEIN.out.summary.collect(), RUN_PROTEIN.out.mr.collect())

    def keys = ['window_kb', 'info_min', 'maf_min', 'p_threshold', 'f_min', 'clump_r2', 'clump_kb',
                'p1', 'p2', 'p12', 'fdr', 'outcome_ncase', 'outcome_ncontrol']
    def pj = file("${workflow.workDir}/run_params.json")
    pj.text = groovy.json.JsonOutput.toJson(keys.collectEntries { [(it): params[it]] } +
        [outcome: file(params.outcome).name, ld_ref: params.ld_ref ? file(params.ld_ref).name : 'none (lead variant)'])
    REPORT(AGGREGATE.out.results, RUN_PROTEIN.out.region.flatten().collect(), pj)
}
