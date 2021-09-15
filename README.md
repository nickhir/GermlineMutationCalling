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

At the end of the worklow, the [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) is used to annotate the identified germline mutations.
 

A high level overview of the performed steps can be seen below: 
![DAG](imgs/dag.svg)
As seen by the execution graph, an arbitrary number of samples/BAM files 
can be processed in parallel up to the joint variant calling.
## Installation 
-> required packages

## Usage

## Output