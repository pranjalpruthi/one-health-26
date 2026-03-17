
# 🧬 SNP Calling Pipeline using BCFtools (with Duplicate Marking)

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
bwa mem -t 8 ref.fa SRR37561623_1.fastq SRR37561623_2.fastq > covid.sam
```
---

### 3️⃣ Convert SAM → BAM and Sort

```
samtools view -@ 8 -bS covid.sam > covid.bam
samtools sort -@ 8 -o covid.sorted.bam covid.bam
```
---

### 4️⃣ Fix Mate Information

```
samtools fixmate -m sample.sorted.bam sample.fixmate.bam
```

**Why?**
- Adds mate-pair information needed for accurate duplicate marking  
- Ensures proper pairing statistics  

---

### 5️⃣ Re-sort by Position

```bash
samtools sort -@ 8 -o sample.positionsort.bam sample.fixmate.bam
```

**Why?**
- Required before duplicate marking  
- Ensures reads are ordered by genomic position  

---

### 6️⃣ Mark PCR Duplicates

```bash
samtools markdup -@ 8 sample.positionsort.bam sample.markdup.bam
samtools index sample.markdup.bam
```

**Why?**
- Removes PCR amplification bias  
- Prevents false SNP calls due to artificially duplicated reads  

---

### 7️⃣ Variant Calling

```bash
bcftools mpileup -Ou -f ref.fa sample.markdup.bam | \
bcftools call -mv -Oz -o sample.vcf.gz
```

**Why?**
- `mpileup`: summarizes base information at each position  
- `call`: identifies SNPs and variants  

---

### 8️⃣ Index the VCF

```bash
bcftools index sample.vcf.gz
```

**Why?**
- Enables fast querying of the VCF  

---

### 9️⃣ Extract SNPs

```bash
bcftools view -v snps sample.vcf.gz -Oz -o sample.snps.vcf.gz
bcftools index sample.snps.vcf.gz
```

**Why?**
- Keeps only SNPs (removes indels and other variants)  

---

### 🔟 Filter SNPs

```bash
bcftools filter -e 'QUAL<30 || DP<10' sample.snps.vcf.gz -Oz -o sample.filtered.vcf.gz
bcftools index sample.filtered.vcf.gz
```

**Why?**
- Removes low-quality variants  
- Improves reliability of downstream analysis  

---

## 📊 Optional: Multi-sample Calling

```bash
bcftools mpileup -Ou -f ref.fa *.markdup.bam | \
bcftools call -mv -Oz -o cohort.vcf.gz
```

**Why?**
- Calls variants jointly across samples  
- Improves consistency across datasets  

---

## 🧪 Output Files

| File | Description |
|------|-------------|
| sample.markdup.bam | Final processed alignment |
| sample.vcf.gz | Raw variants |
| sample.snps.vcf.gz | SNP-only variants |
| sample.filtered.vcf.gz | High-quality SNPs |

---

## ⚠️ Notes

- Duplicate marking is **critical for diploid genomes**
- Filtering thresholds (`QUAL`, `DP`) should be adjusted based on dataset
- For very large datasets, use pipes (`-Ou`) to reduce I/O overhead

---

## 🧹 Cleanup (Optional)

```bash
rm sample.sam sample.bam sample.sorted.bam sample.fixmate.bam sample.positionsort.bam
```

---

## 🚀 Tips

- Use multiple threads (`-t`, `-@`)
- Avoid writing intermediate files when possible
- For HPC: parallelize across samples, not steps
