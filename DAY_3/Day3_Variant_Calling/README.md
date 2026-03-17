# 🧬 SNP Calling Pipeline using BCFtools

This pipeline describes how to perform read alignment and variant (SNP) calling from raw FASTQ files using bwa, samtools, and bcftools.

## 📦 Requirements
```
sra-tools
bwa  
samtools  
bcftools  
tabix
```

## 📁 Input Data
We'll be downloading both our reference genome and our test data.

#### Reference genome
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz 
gunzip GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz    
```
#### Paired-end reads 
```
fasterq-dump --split-files SRR37561623/SRR37561623.sra
```
That will give you paired-end Illumina files, containing SARS-CoV-2 reads, and your directory should consist of 3 files:
```
ls
GCF_009858895.2_ASM985889v3_genomic.fna  SRR37561623_1.fastq  SRR37561623_2.
```
We'll rename our reference genome to make it easier for use.
```
mv GCF_009858895.2_ASM985889v3_genomic.fna ref.fa
```

## 🧭 Workflow Overview

1. Index reference genome  
2. Align reads to reference  
3. Convert SAM → BAM, sort, index
4. Mark PCR duplicates
5. Generate pileup  
6. Call variants  
7. Filter SNPs  

## 🔧 Step-by-Step Pipeline

### 1️⃣ Index the Reference Genome

```
bwa index ref.fa  
samtools faidx ref.fa  
```
### 2️⃣ Align Reads with BWA-MEM

```
bwa mem -t 8 ref.fa SRR37561623_1.fastq SRR37561623_2.fastq > covid.sam  
```

### 3️⃣ Convert SAM to BAM, Sort and Index

```
samtools view -bS sample.sam > sample.bam  
samtools sort -@ 8 -o sample.sorted.bam sample.bam  
samtools index sample.sorted.bam
```

### 4️⃣ Generate Pileup and Call Variants

bcftools mpileup -Ou -f ref.fa sample.sorted.bam | bcftools call -mv -Oz -o sample.vcf.gz  

### 5️⃣ Index the VCF File

bcftools index sample.vcf.gz  

### 6️⃣ Filter SNPs

bcftools view -v snps sample.vcf.gz -Oz -o sample.snps.vcf.gz  
bcftools index sample.snps.vcf.gz  

bcftools filter -e 'QUAL<30 || DP<10' sample.snps.vcf.gz -Oz -o sample.filtered.vcf.gz  
bcftools index sample.filtered.vcf.gz  

## 📊 Optional: Multi-sample Variant Calling

bcftools mpileup -Ou -f ref.fa *.sorted.bam | bcftools call -mv -Oz -o cohort.vcf.gz  

## 🧪 Output Files

- sample.sorted.bam  
- sample.vcf.gz  
- sample.snps.vcf.gz  
- sample.filtered.vcf.gz  

## 🧹 Cleanup (Optional)

rm *.bam.bai  

## 📌 Notes

- Suitable for haploid and diploid organisms  
- Adjust filtering thresholds as needed  
