
# 🧬 SNP Calling Pipeline using BCFtools (with Duplicate Marking)

This pipeline performs **read alignment and SNP calling** from raw FASTQ files using:
- `bwa`
- `samtools`
- `bcftools`

Each step includes **commands + explanation** so you know *why* you're doing it.

---

## 📦 Requirements

```
bwa
samtools
bcftools
tabix (optional)
```

---

## 📁 Input Data

- Reference genome: `ref.fa`
- Reads:
  - `sample_R1.fastq.gz`
  - `sample_R2.fastq.gz`

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

## 🔧 Step-by-Step Pipeline

---

### 1️⃣ Index the Reference Genome

```bash
bwa index ref.fa
samtools faidx ref.fa
```

**Why?**
- `bwa index`: creates an index so reads can be aligned efficiently  
- `samtools faidx`: allows fast random access to the reference  

---

### 2️⃣ Align Reads with BWA-MEM

```bash
bwa mem -t 8 ref.fa sample_R1.fastq.gz sample_R2.fastq.gz > sample.sam
```

**Why?**
- Maps sequencing reads to the reference genome  
- Produces a SAM file with alignment information  

---

### 3️⃣ Convert SAM → BAM and Sort

```bash
samtools view -bS sample.sam > sample.bam
samtools sort -@ 8 -o sample.sorted.bam sample.bam
```

**Why?**
- BAM is compressed → faster and smaller  
- Sorting is required for downstream tools  

---

### 4️⃣ Fix Mate Information

```bash
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
