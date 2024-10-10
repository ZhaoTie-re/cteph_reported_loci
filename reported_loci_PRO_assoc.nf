params.vcfPath = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/CTEPH/08.addTommo_HGVD_AF_vcf/"
params.nagavcfPath = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/NAGAHAMA/08.addTommo_HGVD_AF_vcf/"
params.protPath = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Proteome/case_control_prepare"
params.scriptPath = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_reported_loci/scripts"

params.annoFile = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Proteome/case_control_prepare/annot_protein.csv"
params.lociKnown = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_reported_loci/datasrc/cteph_reported_loci.xlsx"
params.outDir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_reported_loci/result_loci_PRO_assoc"

process reported_loci_withdrawn {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/01.reported_loci", mode: 'symlink'
    
    input:
    val vcfPath from params.vcfPath
    path lociKnown from params.lociKnown
    
    output:
    tuple file(vcf), file(vcf_tbi) into reported_loci_ch, case_stats_ch
    
    script:
    vcf = 'reported_loci.vcf.gz'
    vcf_tbi = 'reported_loci.vcf.gz.tbi'
    """
    source activate cteph_geno_pro
    python ${params.scriptPath}/loci_withdrawn.py --vcfPath ${params.vcfPath} --lociKnown ${params.lociKnown} --outputFile ${vcf}
    """
}

process GT_prepare {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/02.GT_prepare", mode: 'symlink'
    
    input:
    tuple file(vcf), file(vcf_tbi) from reported_loci_ch
    
    output:
    file(GTMatrix) into GTMatrix_ch
    file(nocalFigure)
    file('reported_loci_ls.txt') into lociList_ch
    
    script:
    GTMatrix = 'GTmatrix_df.csv'
    nocalFigure = 'nocall.pdf'
    """
    source activate cteph_geno_pro
    python ${params.scriptPath}/GT_prepare.py --vcfPath $vcf --GTMatrix $GTMatrix --nocallFigure $nocalFigure
    """
}

lociList_ch
    .splitText()
    .combine(GTMatrix_ch)
    .set { combined_ch }

process PRO_GT_assoc {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/03.PRO_GT_assoc", mode: 'symlink'
    
    input:
    tuple loci, file(GTMatrix) from combined_ch

    output:
    tuple loci, file('*_PRO_result.csv') into assocTable
    file('*_QQ_plot.pdf')
    
    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptPath}/PRO_GT_assoc.py --GTMatrix $GTMatrix --protPath ${params.protPath} --loci $loci
    """
}

process PRO_GT_assoc_summary {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/04.PRO_GT_assoc_summary", mode: 'symlink'
    
    input:
    tuple loci, file(Pro_csv) from assocTable
    path lociKnown from params.lociKnown
    path annoFile from params.annoFile
    
    output:
    tuple loci, file('*_PRO_result_anno.csv') into assocSummary
    file('*_PRO_result_qq_plot.html')
    
    script:
    """
    source activate cteph_geno_pro

    python ${params.scriptPath}/PRO_GT_assoc_summary.py \\
    --loci ${loci.trim()} \\
    --Pro_ass_file ${Pro_csv} \\
    --reported_loci ${lociKnown} \\
    --anno_file ${annoFile}
    """
}

process reported_loci_naga {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/05.reported_loci_naga", mode: 'symlink'
    
    input:
    val vcfPath from params.nagavcfPath
    path lociKnown from params.lociKnown
    
    output:
    tuple file(vcf), file(vcf_tbi) into reported_loci_naga_ch, control_stats_ch
    
    script:
    vcf = 'reported_loci_naga.vcf.gz'
    vcf_tbi = 'reported_loci_naga.vcf.gz.tbi'
    """
    source activate cteph_geno_pro
    python ${params.scriptPath}/loci_withdrawn.py --vcfPath ${vcfPath} --lociKnown ${lociKnown} --outputFile ${vcf}
    """
}

process GT_prepare_naga {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/06.GT_prepare_naga", mode: 'symlink'
    
    input:
    tuple file(vcf), file(vcf_tbi) from reported_loci_naga_ch
    
    output:
    file(GTMatrix_naga) into GTMatrix_naga_ch
    file(nocalFigure_naga)
    file('reported_loci_ls.txt') into lociList_naga_ch
    
    script:
    GTMatrix_naga = 'GTmatrix_df_naga.csv'
    nocalFigure_naga = 'nocall_naga.pdf'
    """
    source activate cteph_geno_pro
    python ${params.scriptPath}/GT_prepare.py --vcfPath $vcf --GTMatrix $GTMatrix_naga --nocallFigure $nocalFigure_naga
    """
}

process case_control_stats {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/07.case_control_stats", mode: 'symlink'
    
    input:
    tuple file(case_vcf), file(case_vcf_tbi) from case_stats_ch
    tuple file(control_vcf), file(control_vcf_tbi) from control_stats_ch
    path lociKnown from params.lociKnown

    output:
    file('case_control_stats.csv')
    
    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptPath}/CaseControl_stats.py --case_vcfPath $case_vcf --control_vcfPath $control_vcf --lociKnown ${lociKnown}
    """
}

lociList_naga_ch
    .splitText()
    .combine(GTMatrix_naga_ch)
    .set{ combined_naga_ch }


process PRO_GT_naga_assoc {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/08.PRO_GT_naga_assoc", mode: 'symlink'
    
    input:
    tuple loci, file(GTMatrix_naga) from combined_naga_ch

    output:
    tuple loci, file('*_PRO_result.csv') into assocTable_naga
    file('*_QQ_plot.pdf')
    
    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptPath}/PRO_GT_assoc.py --GTMatrix $GTMatrix_naga --protPath ${params.protPath} --loci $loci
    """
}

process PRO_GT_naga_assoc_summary {
    
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/09.PRO_GT_naga_assoc_summary", mode: 'symlink'
    
    input:
    tuple loci, file(Pro_csv) from assocTable_naga
    path lociKnown from params.lociKnown
    path annoFile from params.annoFile
    
    output:
    tuple loci, file('*_PRO_result_anno.csv') into assocSummary_naga
    file('*_PRO_result_qq_plot.html')
    
    script:
    """
    source activate cteph_geno_pro

    python ${params.scriptPath}/PRO_GT_assoc_summary.py \\
    --loci ${loci.trim()} \\
    --Pro_ass_file ${Pro_csv} \\
    --reported_loci ${lociKnown} \\
    --anno_file ${annoFile}
    """
}

assocSummary
    .map { it -> ["\"${it[0].trim()}\"", "\"${it[1]}\""] }
    .collect(flat:false)
    .map { it -> ["phom", it] }
    .set { formatted_ass_ph }

assocSummary_naga
    .map { it -> ["\"${it[0].trim()}\"", "\"${it[1]}\""] }
    .collect(flat:false)
    .map { it -> ["naga", it] }
    .set { formatted_ass_naga }

formatted_ass_ph
    .mix(formatted_ass_naga)
    .set{ formatted_ass }


process pca_beta {

    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/10.beta_pca", mode: 'symlink'

    input:
    tuple val(group), val(assocSummaryLsit) from formatted_ass
    path(info_path) from params.lociKnown

    output:
    tuple val(group), file(beta_csv) into loci_Pro_beta
    file(pca_html)

    script:
    beta_csv = "${group}_beta_values.csv"
    pca_html = "${group}_pca_beta.html"
    """
    source activate cteph_geno_pro

    python ${params.scriptPath}/pca_beta_group.py \\
    --group ${group} \\
    --assocSummary '${assocSummaryLsit}' \\
    --info_path ${info_path}
    """
}

loci_Pro_beta.view()
