# PERFORME VARIANT QUALITY SCORE RECALIBRATION
# Load the predefined resource requirments ("config/resources.yaml")
import yaml

with open("config/resources.yaml", "r") as f:
    resources = yaml.safe_load(f)

rule VQSR_indel:
    input:
         sites_only_vcf=cohort_dir + "/cohort.sites.only.vcf.gz",
         GATK=GATK,
         mills_indels=mills_indels,
         dbSNP=dbSNP,
         axiom_poly=axiom_poly

    output:
          cohort_indel_recal=cohort_dir + "/filtration/cohort_indel.recal",
          tranche_files=cohort_dir + "/filtration/cohort_indel.tranches"

    benchmark:
             benchmark=cohort_dir + "/benchmark/VQSR_indel.txt"

    threads: resources["threads"]["VQSR_indel"]

    resources:
             runtime_min=resources["runtime"]["VQSR_indel"],
             mem_mb=resources["memory"]["VQSR_indel"]

    params:
          truth_sensitivity_tranches="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 "
                                     "-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 "
                                     "-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0",
          cohort_dir=cohort_dir

    log:
       VQSR_indel=cohort_dir + "/logs/VQSR_indel.out"

    shell:
         "{input.GATK} VariantRecalibrator -V {input.sites_only_vcf} --trust-all-polymorphic "
         "{params.truth_sensitivity_tranches} -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP "
         "-mode INDEL "
         "--max-gaussians 4 "
         "-AS "
         "--resource:mills,known=false,training=true,truth=true,prior=12 {input.mills_indels} "
         "--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 {input.axiom_poly} "
         "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbSNP} "
         "-O {output.cohort_indel_recal} "
         "--tranches-file {output.tranche_files} &> {log.VQSR_indel} \n"

rule VQSR_snp:
    input:
         sites_only_vcf=cohort_dir + "/cohort.sites.only.vcf.gz",
         GATK=GATK,
         dbSNP=dbSNP,
         hapmap=hapmap,
         omni=omni,
         thousandG_phase1_snps_high_confidence=thousandG_phase1_snps_high_confidence

    output:
          cohort_snp_recal=cohort_dir + "/filtration/cohort_snp.recal",
          tranche_files=cohort_dir + "/filtration/cohort_snp.tranches"

    benchmark:
             benchmark=cohort_dir + "/benchmark/VQSR_snp.txt"

    threads: resources["threads"]["VQSR_snp"]

    resources:
             runtime_min=resources["runtime"]["VQSR_snp"],
             mem_mb=resources["memory"]["VQSR_snp"]

    params:
          truth_sensitivity_tranches="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 "
                                     "-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 "
                                     "-tranche 97.0 -tranche 90.0",
          cohort_dir=cohort_dir

    log:
       VQSR_snp=cohort_dir + "/logs/VQSR_snp.out"

    shell:
         "{input.GATK} VariantRecalibrator -V {input.sites_only_vcf} --trust-all-polymorphic "
         "{params.truth_sensitivity_tranches} -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP "
         "-mode SNP "
         "--max-gaussians 6 "
         "-AS "
         "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} "
         "--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} "
         "--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG_phase1_snps_high_confidence} "
         "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbSNP} "
         #"--rscript-file {params.cohort_dir}/filtration/diagnostic.plots.snp.R "
         "-O {output.cohort_snp_recal} "
         "--tranches-file {output.tranche_files} &> {log.VQSR_snp}"

rule ApplyVQSR_indel:
    input:
         merged_cohort_vcf=cohort_dir + "/cohort.vcf.gz",
         GATK=GATK,
         cohort_indel_recal=cohort_dir + "/filtration/cohort_indel.recal",
         tranche_files=cohort_dir + "/filtration/cohort_indel.tranches",

    output:
          indel_recalibrated_cohort=cohort_dir + "/cohort.indel.recalibrated.vcf.gz"

    benchmark:
             benchmark=cohort_dir + "/benchmark/ApplyVQSR_indel.txt"

    threads: resources["threads"]["ApplyVQSR_indel"]

    resources:
             runtime_min=resources["runtime"]["ApplyVQSR_indel"],
             mem_mb=resources["memory"]["ApplyVQSR_indel"]

    log:
       ApplyVQSR_indel=cohort_dir + "/logs/ApplyVQSR_indel.out"

    params:
          cohort_dir=cohort_dir

    shell:
         "{input.GATK} ApplyVQSR -V {input.merged_cohort_vcf} --recal-file {input.cohort_indel_recal} "
         "--tranches-file {input.tranche_files} --truth-sensitivity-filter-level 97 "
         "--create-output-variant-index true -mode INDEL "
         "-AS "
         "-O {output.indel_recalibrated_cohort} &> {log.ApplyVQSR_indel} \n"
         "rm {cohort_dir}/cohort.sites.only.vcf.gz.tbi"

rule ApplyVQSR_snp:
    input:
         indel_recalibrated_cohort=cohort_dir + "/cohort.indel.recalibrated.vcf.gz",
         GATK=GATK,
         cohort_snp_recal=cohort_dir + "/filtration/cohort_snp.recal",
         tranche_files=cohort_dir + "/filtration/cohort_snp.tranches",

    output:
          recalibrated_cohort=cohort_dir + "/cohort.recalibrated.vcf.gz"

    benchmark:
             benchmark=cohort_dir + "/benchmark/ApplyVQSR_snp.txt"

    threads: resources["threads"]["ApplyVQSR_snp"]

    resources:
             runtime_min=resources["runtime"]["ApplyVQSR_snp"],
             mem_mb=resources["memory"]["ApplyVQSR_snp"]

    log:
       ApplyVQSR_snp=cohort_dir + "/logs/ApplyVQSR_snp.out"

    shell:
         "{input.GATK} ApplyVQSR -V {input.indel_recalibrated_cohort} --recal-file {input.cohort_snp_recal} "
         "--tranches-file {input.tranche_files} --truth-sensitivity-filter-level 97 "
         "--create-output-variant-index true -mode SNP "
         "-AS "
         "-O {output.recalibrated_cohort} &> {log.ApplyVQSR_snp}"

rule SelectVariants:
    input:
         bcftools=bcftools,
         filtered_vcf=cohort_dir + "/cohort.recalibrated.vcf.gz"

    output:
          pass_only_vcf=temp(cohort_dir + "/cohort.recalibrated.pass.vcf.gz")

    log:
       err_SelectVariants=cohort_dir + "/logs/SelectVariants.err"

    benchmark:
             benchmark=cohort_dir + "/benchmark/SelectVariants.txt"

    threads: resources["threads"]["SelectVariants"]

    resources:
             runtime_min=resources["runtime"]["SelectVariants"],
             mem_mb=resources["memory"]["SelectVariants"]

    shell:
         "{input.bcftools} view -f PASS {input.filtered_vcf} -o {output.pass_only_vcf} -O z &> {log.err_SelectVariants} \n"