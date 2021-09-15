# USED TO PERFORME BASE QUALITY SCORE RECALIBRATION
# -------------------------------------------------
# Load the predefined resource requirments ("config/resources.yaml")
import yaml
with open("config/resources.yaml","r") as f:
    resources = yaml.safe_load(f)

# GATKs multiprocessing BaseRecalibrator is still in testing.
# We can parallize it by ourself, by spliting the genome into multiple regions (25) and analyzing them simultaneously.
# For this we have to define some helper functions which gather the results of the distributed analysis.

def get_GatherBQSRReports_files(wildcards):
    '''
    We have to combine the scattered BQSR reports. Seperate them with a '-I' because
    GatherBQSRReports expects the input like this
    '''
    import re
    with open(intervals, "r") as f:
        # get the interval name (everything after the last "/")
        content = f.read().splitlines()
        content = [re.search(r"([^\/]+$)", i).group(1) for i in content]
    recal_tables = [f"-I {tmp_dir}/{wildcards.patient}/normal_recal_bqsr.{interval}.table " for interval in content]
    return "".join(recal_tables)  # creates one large input string


def get_GatherRecalBamFiles_files(wildcards):
    '''
    We have to combine the scattered BAM files. Seperate them with a '-I' because
    GatherBamFiles expects the input like this
    '''
    import re
    with open(intervals, "r") as f:
        # get the interval name (everything after the last "/")
        content = f.read().splitlines()
        content = [re.search(r"([^\/]+$)", i).group(1) for i in content]
    bams = [f"-I {tmp_dir}/{wildcards.patient}/recal.{interval}.bam " for interval in content]
    return "".join(bams)  # creates one large input string


rule BaseRecalibrator:
    input:
         GATK=GATK,
         reference=reference,
         dbSNP=dbSNP,
         intervals=intervals,
         marked_dup_bam=lambda wildcards: samples[wildcards.patient][0]

    output:
          decoy=temp(tmp_dir + "/{patient}/BaseRecalibrator_finished")

    benchmark:
             benchmark=path + "/{patient}/benchmark/BaseRecalibrator.txt"

    priority: 30

    threads: resources["threads"]["BaseRecalibrator"]

    resources:
             runtime_min=resources["runtime"]["BaseRecalibrator"],
             mem_mb=resources["memory"]["BaseRecalibrator"]

    params:
          recal_data=tmp_dir + "/{patient}/normal_recal_bqsr",

    log:
       out_recal=path + "/{patient}/logs/BaseRecalibrator"

    shell:
         "mem_per_job=$(expr {resources.mem_mb} / {threads}) \n"
         "export mem_per_job \n"

         "cat {input.intervals} | parallel --plus -j {threads} '{input.GATK} --java-options -Xmx${{mem_per_job}}m BaseRecalibrator "
         "-I {input.marked_dup_bam} "
         "-R {input.reference} "
         "--known-sites {input.dbSNP} "
         "-L {{}} "
         "-O {params.recal_data}.{{/}}.table &> {log.out_recal}.{{/}}.out' \n"

         "touch {output.decoy}"

rule GatherBQSRReports:
    input:
         decoy=tmp_dir + "/{patient}/BaseRecalibrator_finished",
         GATK=GATK

    output:
          recal_data=temp(tmp_dir + "/{patient}/normal_recal_bqsr.combined.table"),

    benchmark:
             benchmark=path + "/{patient}/benchmark/GatherBQSRReports.txt"


    priority: 35

    threads: resources["threads"]["GatherBQSRReports"]

    resources:
             runtime_min=resources["runtime"]["GatherBQSRReports"],
             mem_mb=resources["memory"]["GatherBQSRReports"]

    params:
          recal_tables=get_GatherBQSRReports_files,
          tmp_dir=tmp_dir

    log:
       out_GatherBQSRReports=path + "/{patient}/logs/GatherBQSRReports.out",

    shell:
         "{input.GATK} GatherBQSRReports "
         "{params.recal_tables} "
         "-O {output.recal_data} &> {log.out_GatherBQSRReports} \n"

         "rm {params.tmp_dir}/{wildcards.patient}/normal_recal_bqsr.0*table"

rule ApplyBQSR:
    input:
         GATK=GATK,
         intervals=intervals,
         reference=reference,
         recal_data=tmp_dir + "/{patient}/normal_recal_bqsr.combined.table",
         marked_dup_bam=lambda wildcards: samples[wildcards.patient][0],

    output:
          decoy=temp(tmp_dir + "/{patient}/ApplyBQSR_finished")

    benchmark:
             benchmark=path + "/{patient}/benchmark/ApplyBQSR.txt"

    threads: resources["threads"]["ApplyBQSR"]

    resources:
             runtime_min=resources["runtime"]["ApplyBQSR"],
             mem_mb=resources["memory"]["ApplyBQSR"]
    params:
          recal_bam=tmp_dir + "/{patient}/recal",

    priority: 40

    log:
       out_ApplyBQSR=path + "/{patient}/logs/ApplyBQSR"

    shell:
         "mem_per_job=$(expr {resources.mem_mb} / {threads}) \n"
         "export mem_per_job \n"
         "cat {input.intervals} | parallel --plus -j {threads} '{input.GATK} --java-options -Xmx${{mem_per_job}}m ApplyBQSR "
         "-R {input.reference} "
         "-L {{}} "
         "-I {input.marked_dup_bam} "
         "--bqsr-recal-file {input.recal_data} "
         "-O {params.recal_bam}.{{/}}.bam &> {log.out_ApplyBQSR}.{{/}}.out' \n"

         "touch {output.decoy}"

rule GatherRecalBamFiles:
    input:
         decoy=tmp_dir + "/{patient}/ApplyBQSR_finished",
         GATK=GATK

    output:
          final_recal_bam=temp(tmp_dir + "/{patient}/{patient}.normal.recal.bam"),

    benchmark:
             benchmark=path + "/{patient}/benchmark/GatherRecalBamFiles.txt"

    threads: resources["threads"]["GatherRecalBamFiles"]

    resources:
             runtime_min=resources["runtime"]["GatherRecalBamFiles"],
             mem_mb=resources["memory"]["GatherRecalBamFiles"]

    priority: 46

    params:
          bams=get_GatherRecalBamFiles_files,
          path=path,
          tmp_dir=tmp_dir

    log:
       out_GatherRecalBamFiles=path + "/{patient}/logs/GatherRecalBamFiles.out"

    shell:
         "{input.GATK} GatherBamFiles "
         "{params.bams} "
         "-O {output.final_recal_bam} &> {log.out_GatherRecalBamFiles} \n"

         # remove the intermediate files
         "rm {params.tmp_dir}/{wildcards.patient}/recal.0*-scattered.interval_list.ba*"

rule SortBam:
    input:
         final_recal_bam=tmp_dir + "/{patient}/{patient}.normal.recal.bam",
         samtools=samtools

    output:
          recal_sorted=(tmp_dir + "/{patient}/{patient}.normal.recal.sorted.bam"),

    benchmark:
             benchmark=path + "/{patient}/benchmark/SortBam.txt"

    threads: resources["threads"]["SortBam"]*3

    resources:
             runtime_min=resources["runtime"]["SortBam"],
             mem_mb=resources["memory"]["SortBam"]

    priority: 47

    log:
       out_SortBam=path + "/{patient}/logs/SortBam.out"

    shell:
         "{input.samtools} sort -@ {threads} -o {output.recal_sorted} {input.final_recal_bam} &> {log.out_SortBam}"

rule IndexBam:
    input:
         recal_sorted=tmp_dir + "/{patient}/{patient}.normal.recal.sorted.bam",
         sambamba=sambamba

    output:
          recal_sorted_index=(tmp_dir + "/{patient}/{patient}.normal.recal.sorted.bai"),

    benchmark:
             benchmark=path + "/{patient}/benchmark/IndexBam.txt"

    threads: resources["threads"]["IndexBam"]

    resources:
             runtime_min=resources["runtime"]["IndexBam"],
             mem_mb=resources["memory"]["IndexBam"]

    priority: 48

    log:
       out_IndexBam=path + "/{patient}/logs/IndexBam.out"

    shell:
         "{input.sambamba} index -p -t {threads} {input.recal_sorted} {output.recal_sorted_index} &> {log.out_IndexBam}"
