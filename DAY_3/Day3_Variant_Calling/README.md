
# 🧬 SNP Calling Pipeline using BCFTools

This workflow performs **read alignment and SNP calling** from FASTQ files using:

## 📦 Dependencies required

```
- sra-tools
- bwa
- samtools
- bcftools
- tabix
```

---

## Conda Environment Activation
```
conda activate bioset1
```

## 📁 Input Data

#### Reference genome:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
gunzip GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
```
We will also rename our reference genome to make it easier to process.
```
mv GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna ref.fa
```
#### Paired-end reads:
```
prefetch SRR37561623
fasterq-dump --split-files SRR37561623/SRR37561623.sra
```
---

## 🧭 Workflow Overview

1. Index reference genome  
2. Align reads  
3. Convert + sort BAM  
4. Fix mate information  
5. Mark duplicates  
6. Variant calling  
7. Filtering  

---
---

### 1️⃣ Index the Reference Genome

```
bwa index ref.fa
samtools faidx ref.fa
```
---

### 2️⃣ Align Reads with BWA-MEM

```bash
bwa mem -t 8 -R '@RG\tID:covid\tSM:covid\tPL:illumina' ref.fa SRR37561623_1.fastq SRR37561623_2.fastq > covid.sam
```
---

### 3️⃣ Convert SAM → BAM and Sort

```
samtools view -@ 8 -bS covid.sam > covid.bam
samtools sort -@ 8 -o covid.querysort.bam covid.bam -n
```
---

### 4️⃣ Fix Mate Information

```
samtools fixmate -@ 8 -m covid.querysort.bam covid.fixmate.bam
```
---

### 5️⃣ Re-sort by Position

```
samtools sort -@ 8 -o covid.positionsort.bam covid.fixmate.bam
```
---

### 6️⃣ Mark PCR Duplicates

```bash
samtools markdup -@ 8 covid.positionsort.bam covid.markdup.bam
samtools index covid.markdup.bam
```
---

### 7️⃣ Variant Calling

```bash
bcftools mpileup --threads 8 -Ou -f ref.fa covid.markdup.bam | bcftools call --threads 8 -mv -Oz -o covid.vcf.gz
bcftools index covid.vcf.gz

```

### 8️⃣ Extract SNPs

```bash
bcftools view -v snps covid.vcf.gz -Oz -o covid.snps.vcf.gz
bcftools index covid.snps.vcf.gz
```
---

### 9️⃣ Filter SNPs

```bash
bcftools filter -e 'QUAL<30 || DP<10' covid.snps.vcf.gz -Oz -o covid.filtered.vcf.gz
bcftools index covid.filtered.vcf.gz
```

## 🧪 Output Files

| File | Description |
|------|-------------|
| covid.markdup.bam | Final processed alignment |
| covid.vcf.gz | Raw variants |
| covid.snps.vcf.gz | SNP-only variants |
| covid.filtered.vcf.gz | High-quality SNPs |

---
