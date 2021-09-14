import yaml
import subprocess

# load in the config file which contains all the information needed for the workflow
with open("config/config_wgs.yaml") as f:
    config = yaml.safe_load(f)

### PROGRAMS USED
GATK = config["GATK"]
bcftools = config["bcftools"]
sambamba = config["sambamba"]
samtools = config["samtools"]
vep = config["vep"]
vep_cache = config["vep_cache"]

### OUTPUT DIRECTORIES
tmp_dir = config["tmp_dir"]
path = config["path"]
cohort_dir = config["cohort_dir"]

### GENOMIC INTERVALS OVER WHICH OPERATIONS ARE PARALLIZED
intervals = config["intervals"]
chromosomes = config["chromosomes"]

### REFERENCE RESOURCES
reference = config["reference"]
reference_dict = config["reference_dict"]
dbSNP = config["dbSNP"]
hapmap = config["hapmap"]
omni = config["omni"]
thousandG_phase1_snps_high_confidence = config["thousandG_phase1_snps_high_confidence"]
mills_indels = config["mills_indels"]
axiom_poly = config["axiom_poly"]

### ACTUAL INPUT DATA
patient_id = config["patient_id"]
bam_files = config["bam_files"]

# Each bam file belongs to one patient. We create a dictionary which does that mapping
samples = {}
for i in range(len(patient_id)):
    samples[patient_id[i]] = [bam_files[i]]


rule all:
    input:
         cohort_dir + "/cohort.recalibrated.pass.vep.vcf.gz"

include: "rules/BaseQualityScoreRecalibration.smk"
include: "rules/JointGenotyping.smk"
include: "rules/VQSR.smk"
include: "rules/VEP.smk"
