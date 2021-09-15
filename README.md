# Germline Mutation Calling
This Snakemake workflow follows the
[GATK best-practice recommandations](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) 
to call small germline variants.

The pipeline requires as inputs aligned BAM files (e.g. with [BWA](http://bio-bwa.sourceforge.net/bwa.shtml)) 
where the duplicates are already marked (e.g. with [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
or [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-markdup.html)).
It then performed [Base Quality Score Recalibration](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) 
and [joint genotyping of multiple samples](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-),
which is automatically parallized over user defined intervals (for examples see [intervals.txt](inputs/WGS-interval-files-excluding-supercontigs/intervals.txt)) and chromosomes. 

Filtering is performed using GATKs state-of-the-art [Variant Quality Score Recalibration](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)
 
At the end of the worklow, the [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) is used to annotate the identified germline mutations.
 

A high level overview of the performed steps can be seen below: 

![DAG](imgs/dag.svg)

As seen by the execution graph, an arbitrary number of samples/BAM files 
can be processed in parallel up to the joint variant calling.

## Installation 
Required tools:
- [GATK4](https://github.com/broadinstitute/gatk/) 
- [samtools](http://www.htslib.org/download/) 
- [bcftools](http://www.htslib.org/download/) 
- [sambamba](https://lomereiter.github.io/sambamba/index.html)
- [vep](https://m.ensembl.org/info/docs/tools/vep/script/vep_download.html)
- [snakemake](https://snakemake.readthedocs.io/en/stable/) 

The majority of the listed tools can be quite easily installed with [conda](https://docs.conda.io/en/latest/) which is recommanded. 

## Usage
First, modify the [config_wgs.yaml](config/config_wgs.yaml) and [resources.yaml](config/resources.yaml) files.
Both files contain detailed description what is expected. The [config_wgs.yaml](config/config_wgs.yaml) also contains 
links to some reference resources. Be careful, they are all specific for the GRCh37/hg19/b37 genome assembly. 

After setting up all the config files and installing all tools, you can simply run: 
```bash
snakemake --latency-wait 300 -j 5 --cluster "sbatch --mem={resources.mem_mb} --time {resources.runtime_min} --cpus-per-task {threads} --job-name={rule}.%j --output snakemake_cluster_submit/{rule}.%j.out --mail-type=FAIL"
```
This assumes that the cluster you are using is running [SLURM](https://slurm.schedmd.com/documentation.html).
If this is not the case, you have to adjust the command after `--cluster`.

`-j` specifies the number of jobs/rules should be submitted in parallel.
## Output