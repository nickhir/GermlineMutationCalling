# RUN THE VARIANT EFFECT PREDICTOR
# Load the predefined resource requirments ("config/resources.yaml")
import yaml
with open("config/resources.yaml","r") as f:
    resources = yaml.safe_load(f)

rule VEP:
    input:
         vep=vep,
         reference=reference,
         pass_only_vcf=cohort_dir + "/cohort.recalibrated.pass.vcf.gz"

    output:
          vep_vcf=cohort_dir + "/cohort.recalibrated.pass.vep.vcf.gz"

    log:
       out_VEP=cohort_dir + "/logs/VEP.out"

    benchmark:
             benchmark=cohort_dir + "/benchmark/VEP.txt"

    threads: resources["threads"]["VEP"]

    resources:
             runtime_min=resources["runtime"]["VEP"],
             mem_mb=resources["memory"]["VEP"]

    params:
          vep_cache=vep_cache

    shell:
         "{input.vep} --cache --dir_cache {params.vep_cache} --offline --fasta {input.reference} --vcf --buffer_size 55000 "
         "--pick --fork {threads} --sift b --variant_class --regulatory --check_existing --af_1kg --compress_output bgzip "
         "-i {input.pass_only_vcf} -o {output.vep_vcf} &> {log.out_VEP} \n"
