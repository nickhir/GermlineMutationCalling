# USED TO PERFORM JOINT GENOTYPING
# -------------------------------------------------
# Load the predefined resource requirments ("config/resources.yaml")
import yaml
with open("config/resources.yaml","r") as f:
    resources = yaml.safe_load(f)

# GATKs multiprocessing version of HaplotypeCaller is still in testing.
# GenotypeGVCFs can be performed chromosome by chromosome.
# We can parallize it by ourself, by spliting the genome into multiple regions (25) and analyzing them simultaneously.
# For this we have to define some helper functions which gather the results of the distributed analysis.


def get_gvcfs(wildcards):
    '''
    We have to combine the scattered gvcfs. Seperate them with a '-I' because
    MergeVcfs expects the input like this
    '''
    import re
    with open(intervals, "r") as f:
        # get the interval name (everything after the last "/")
        content = f.read().splitlines()
        content = [re.search(r"([^\/]+$)", i).group(1) for i in content]
    recal_tables = [f"-I {tmp_dir}/{wildcards.patient}/{wildcards.patient}.{interval}.g.vcf.gz " for interval in
                    content]
    return "".join(recal_tables)  # creates one large input string

def get_CombineGVCFs_files(wildcards):
    '''
    We have to combine the g.vcf files per sample in one command. Seperate them with a '-V' because
    CombineGVCFs expects the input like this
    '''
    files = expand(path + "/{patient}/{patient}.germline.merged.g.vcf.gz", patient=patient_id)
    files = [f"-V {file}" for file in files]
    return files

def get_CombineGVCFs_decoys(wildcards):
    '''
    Now we have to combine the g.vcf files per sample in one command. For that we use all the decoys which were created
    in the preceding step. Only if they are all present, the step gets executed.
    '''
    files = expand(path + "/{patient}/{patient}.MergeHaplotypeCaller_finished", patient=patient_id)
    return files

def get_cohort_gvcfs(wildcards):
    '''
    We have to combine the scattered cohort vcfs (we splitted them by chromosome). Seperate them with a '-I' because
    MergeVcfs expects the input like this
    '''
    with open(chromosomes, "r") as f:
        content = f.read().splitlines()
    cohort_gvcfs = [f"-I {cohort_dir}/{chrom}.cohort.vcf.gz " for chrom in content]
    return "".join(cohort_gvcfs)  # creates one large input string

rule HaplotypeCaller:
    input:
         recal_bam=tmp_dir + "/{patient}/{patient}.normal.recal.sorted.bam",
         recal_bam_index=tmp_dir + "/{patient}/{patient}.normal.recal.sorted.bai",
         reference=reference,
         intervals=intervals,
         GATK=GATK

    output:
          decoy=temp(tmp_dir + "/{patient}/HaplotypeCaller_finished")

    benchmark:
             benchmark=path + "/{patient}/benchmark/HaplotypeCaller.txt"

    threads: resources["threads"]["HaplotypeCaller"]

    resources:
             runtime_min=resources["runtime"]["HaplotypeCaller"],
             mem_mb=resources["memory"]["HaplotypeCaller"]

    priority: 50

    params:
          tmp_dir=tmp_dir

    log:
       HaplotypeCaller=path + "/{patient}/logs/HaplotypeCaller"

    shell:
         # determine how many intervals we can process in parallel
         "parallel_jobs=$(( {threads} / 3)) \n"
         "export parallel_jobs \n"
         
         "mem_per_job=$(({resources.mem_mb} / ${{parallel_jobs}})) \n"
         "export mem_per_job \n"

         "cat {input.intervals} | parallel --plus -j $parallel_jobs '{input.GATK} --java-options -Xmx${{mem_per_job}}m HaplotypeCaller  "
         "-R {input.reference} "
         "-I {input.recal_bam} -L {{}} "
         "-O {params.tmp_dir}/{wildcards.patient}/{wildcards.patient}.{{/}}.g.vcf.gz "
         "-ERC GVCF "
         "--tmp-dir {params.tmp_dir}/{wildcards.patient}/ "
         "--native-pair-hmm-threads 5 "
         "-G StandardAnnotation -G AS_StandardAnnotation "
         "-G StandardHCAnnotation &> {log.HaplotypeCaller}.{{/}}.out' \n"

         "touch {output.decoy}"

rule MergeHaplotypeCaller:
    input:
         decoy=tmp_dir + "/{patient}/HaplotypeCaller_finished",
         GATK=GATK,
         intervals=intervals

    output:
          merged_gvcf=path + "/{patient}/{patient}.germline.merged.g.vcf.gz",
          
          # also create a decoy which ensures that CombineGVCFs is only run if all decoys are there
          decoy=temp(path + "/{patient}/{patient}.MergeHaplotypeCaller_finished")

    benchmark:
             benchmark=path + "/{patient}/benchmark/MergeHaplotypeCaller.txt"

    threads: resources["threads"]["MergeHaplotypeCaller"]

    resources:
             runtime_min=resources["runtime"]["MergeHaplotypeCaller"],
             mem_mb=resources["memory"]["MergeHaplotypeCaller"]

    params:
          tmp_dir=tmp_dir,
          gvcfs=get_gvcfs,

    log:
       MergeHaplotypeCaller=path + "/{patient}/logs/MergeHaplotypeCaller.out"

    shell:
         "{input.GATK} MergeVcfs {params.gvcfs} -O {output.merged_gvcf} &> {log.MergeHaplotypeCaller} \n"
         "rm -rf {params.tmp_dir}/{wildcards.patient} \n"
         "touch {output.decoy}"

rule CombineGVCFs:
    input:
         input_decoy=get_CombineGVCFs_decoys,
         GATK=GATK,
         reference=reference,
         chromosomes=chromosomes

    output:
          # we use a decoy output again to foul snakemake
          decoy=temp(cohort_dir + "/CombineGVCFs_DONE")

    benchmark:
             benchmark=cohort_dir + "/benchmark/CombineGVCFs.txt"

    log:
       CombineGVCFs_log=cohort_dir + "/logs/CombineGVCFs"

    params:
          annotation="-G StandardAnnotation -G AS_StandardAnnotation",
          cohort_dir=cohort_dir,
          input_files=get_CombineGVCFs_files,
          tmp_dir=tmp_dir

    threads: resources["threads"]["CombineGVCFs"]

    resources:
             runtime_min=resources["runtime"]["CombineGVCFs"],
             mem_mb=resources["memory"]["CombineGVCFs"]

    shell:
         "mem_per_job=$(expr {resources.mem_mb} / {threads}) \n"
         "export mem_per_job \n"
         
         "cat {input.chromosomes} | parallel -j {threads} '{input.GATK} CombineGVCFs --java-options -Xmx${{mem_per_job}}m "
         "{params.annotation} -R {input.reference} {params.input_files} "
         "-O {params.cohort_dir}/{{}}.cohort.g.vcf.gz -L {{}} &> {log.CombineGVCFs_log}.{{}}.out' \n"

         "touch {output.decoy}"

rule GenotypeGVCFs:
    input:
         GATK=GATK,
         decoy=cohort_dir + "/CombineGVCFs_DONE",
         reference=reference,
         chromosomes=chromosomes

    output:
          decoy=temp(cohort_dir + "/GenotypeGVCFs_DONE")

    benchmark:
             benchmark=cohort_dir + "/benchmark/GenotypeGVCFs.txt"

    log:
       GenotypeGVCFs_log=cohort_dir + "/logs/GenotypeGVCFs"
    params:
          annotation="-G StandardAnnotation -G AS_StandardAnnotation",
          cohort_dir=cohort_dir

    threads: resources["threads"]["GenotypeGVCFs"]

    resources:
             runtime_min=resources["runtime"]["GenotypeGVCFs"],
             mem_mb=resources["memory"]["GenotypeGVCFs"]

    shell:
         "mem_per_job=$(expr {resources.mem_mb} / {threads}) \n"
         "export mem_per_job \n"
         
         "cat {input.chromosomes} | parallel -j {threads} '{input.GATK} --java-options -Xmx${{mem_per_job}}m GenotypeGVCFs  "
         "{params.annotation} -R {input.reference} -V {params.cohort_dir}/{{}}.cohort.g.vcf.gz "
         "-O {params.cohort_dir}/{{}}.cohort.vcf.gz &> {log.GenotypeGVCFs_log}.{{}}.out' \n"

         "rm {params.cohort_dir}/*g.vcf* \n"
         "touch {output.decoy}"

rule MergeCohortVCFs:
    input:
         decoy=cohort_dir + "/GenotypeGVCFs_DONE",
         GATK=GATK,
         chromosomes=chromosomes

    output:
          merged_cohort_vcf=cohort_dir + "/cohort.vcf.gz",
          sites_only_vcf=(cohort_dir + "/cohort.sites.only.vcf.gz")

    benchmark:
             benchmark=cohort_dir + "/benchmark/MergeCohortVCFs.txt"

    threads: resources["threads"]["MergeCohortVCFs"]

    resources:
             runtime_min=resources["runtime"]["MergeCohortVCFs"],
             mem_mb=resources["memory"]["MergeCohortVCFs"]

    params:
          cohort_dir=cohort_dir,
          cohort_gvcfs=get_cohort_gvcfs


    log:
       MergeCohortVCFs_log=cohort_dir + "/logs/MergeCohortVCFs.out",
       MakeSitesOnlyVcf_log=cohort_dir + "/logs/MakeSitesOnly.out",

    shell:
         "{input.GATK} MergeVcfs {params.cohort_gvcfs} -O {output.merged_cohort_vcf} &> {log.MergeCohortVCFs_log} \n"
         "{input.GATK} MakeSitesOnlyVcf -I {output.merged_cohort_vcf} -O {output.sites_only_vcf} &> {log.MakeSitesOnlyVcf_log} \n"
         "rm {params.cohort_dir}/{{1..22}}.coh* {params.cohort_dir}/X.coh* {params.cohort_dir}/Y.coh*"
