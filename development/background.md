下面按你给的三个流程来拆：**WARP 的 `RNAWithUMIsPipeline` 更偏 RNA-seq/UMI 预处理和 QC；真正的 GATK RNA SNP/Indel calling 是 `gatk4-rnaseq-germline-snps-indels` 或 nf-core/rnavar；产出的 VCF 可以接 PrecisionProDB**，因为 PrecisionProDB v2 的输入本来就包括 VCF/TSV 变异、GTF/GFF3 注释、基因组 FASTA 和蛋白 FASTA，用于生成个体化蛋白数据库。

## 1. 三个流程的定位

**WARP `RNAWithUMIsPipeline.wdl`：适合有 UMI 的 RNA-seq 预处理。**
它支持输入 uBAM 或 paired FASTQ，并要求 `read1Structure/read2Structure` 来描述 UMI/read 结构；流程包括 FASTQ 转 uBAM、UMI 提取、`fastp` adapter trimming、STAR 比对、UMI-aware duplicate marking、RNASeQC/Picard QC、污染估计，输出 aligned BAM、transcriptome BAM 和 QC 文件。这个 WDL 本身**不做 HaplotypeCaller 变异检测**。  

**GATK `gatk4-rnaseq-germline-snps-indels`：这是你要的 RNA SNP/Indel calling 核心。**
官方说明它是 “processing RNA data for germline short variant discovery”，输入期望是 BAM/uBAM，输出包括 BAM、VCF 和 filtered VCF。 其核心 WDL 步骤是：GTF 生成 exon interval → RevertSam/SamToFastq → STAR → MergeBamAlignment → MarkDuplicates → SplitNCigarReads → BaseRecalibrator/ApplyBQSR → HaplotypeCaller → MergeVCFs → VariantFiltration。  

**nf-core/rnavar：更适合批量、可复现、容器化运行。**
nf-core/rnavar 1.3.0 明确是 “GATK4 RNA variant calling pipeline”，流程总结包括 FastQC、可选 UMI extraction、STAR、SAMtools sort/index、Picard MarkDuplicates、GATK SplitNCigarReads、BQSR、HaplotypeCaller、MergeVCFs、Tabix、VariantFiltration，以及 VEP/snpEff/bcftools annotation 和 MultiQC。([nf-core][1]) 它的使用方式是准备 `samplesheet.csv` 后运行 `nextflow run nf-core/rnavar ... --input samplesheet.csv --outdir ... --genome GRCh38`。([nf-core][1])

## 2. GATK RNA SNP/Indel calling 的核心步骤与关键参数

GATK 官方 best-practice 对 RNA-seq 的关键点是：STAR 推荐 two-pass mode；`SplitNCigarReads` 是 RNA 特异性的关键步骤，因为 RNA aligner 会用 CIGAR `N` 表示跨内含子比对，必须转换成 HaplotypeCaller 更适合处理的形式；RNA per-sample calling 不支持 joint calling；HaplotypeCaller 在 `SplitNCigarReads` 后基本不需要特殊 RNA 改造，但推荐把 calling confidence 调到 20；RNA 缺少合适训练 truth set，因此推荐 hard filter，而不是 VQSR/CNNScoreVariants。([GATK][2])

| 阶段          | 工具                               | 关键参数                                                                                                   | 说明                                                                     |
| ----------- | -------------------------------- | ------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------- |
| 下载/准备 FASTQ | SRA Toolkit                      | `fasterq-dump --split-files`                                                                           | 从 SRR run 得到 R1/R2 FASTQ                                               |
| 比对          | STAR                             | `--twopassMode Basic`, `--sjdbGTFfile`, `--sjdbOverhang readLength-1`, `--outSAMattrRGline`            | two-pass 对新剪接位点更好                                                      |
| 去重复         | Picard/GATK `MarkDuplicates`     | `--CREATE_INDEX true`, `--VALIDATION_STRINGENCY SILENT`                                                | 普通 RNA-seq 用 MarkDuplicates；UMI RNA-seq 应用 UMI-aware duplicate marking |
| RNA 特异处理    | GATK `SplitNCigarReads`          | `-R`, `-I`, `-O`                                                                                       | 必跑；处理 CIGAR 中的 `N`                                                     |
| BQSR        | `BaseRecalibrator` + `ApplyBQSR` | `--known-sites dbSNP`, `--known-sites known_indels`, `--use-original-qualities`                        | 参考、GTF、dbSNP/known indels 的 contig 命名必须完全一致                            |
| calling 区域  | GTF exon → interval list         | `BedToIntervalList`, `-L exons.interval_list`                                                          | RNA-seq 建议限制到 exon/transcribed regions                                 |
| 变异检测        | `HaplotypeCaller`                | `-L`, `--dont-use-soft-clipped-bases`, `--standard-min-confidence-threshold-for-calling 20`, `--dbsnp` | WDL 中使用 exon interval、soft-clip 不参与 calling、call confidence 默认 20      |
| 硬过滤         | `VariantFiltration`              | `--window 35 --cluster 3`, `FS > 30.0`, `QD < 2.0`                                                     | GATK RNA WDL/nf-core 的核心 hard filter                                   |
| 注释          | VEP/snpEff/bcftools              | 按参考版本匹配                                                                                                | 用于下游蛋白影响判断；PrecisionProDB 也可直接吃 VCF                                    |

nf-core/rnavar 1.3.0 的 schema 也体现了这些默认值：`gatk_hc_call_conf` 默认 20，interval scatter 默认 25；过滤参数默认 `QD=2`、`FS=30`、`window=35`、`cluster=3`。 

## 3. 假设你有一个人的 SRR paired-end FASTQ：核心运行代码

下面是**简化版本地命令**。它不完全复刻 WDL 的 uBAM→MergeBamAlignment 结构，但保留了 RNA GATK calling 的核心参数。实际项目中，参考 FASTA、GTF、dbSNP、known indels 必须来自同一 genome build，例如全都 GRCh38/GATK resource bundle 或全都同一套 GENCODE/Ensembl 体系。

```bash
set -euo pipefail

THREADS=16
SAMPLE=SRRXXXXXXX
READLEN=150

OUT=gatk_rna_${SAMPLE}
mkdir -p ${OUT}

# 参考文件：示例变量，请替换为你的实际路径
REF=/path/to/Homo_sapiens_assembly38.fasta
DICT=/path/to/Homo_sapiens_assembly38.dict
GTF=/path/to/gencode.v47.annotation.gtf
STAR_INDEX=/path/to/STAR_GRCh38_GENCODE_v47

DBSNP=/path/to/Homo_sapiens_assembly38.dbsnp138.vcf.gz
KNOWN_INDEL1=/path/to/Homo_sapiens_assembly38.known_indels.vcf.gz
KNOWN_INDEL2=/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

### 3.1 下载 SRR FASTQ

```bash
fasterq-dump --split-files --threads ${THREADS} -O ${OUT} ${SAMPLE}
pigz -p ${THREADS} ${OUT}/${SAMPLE}_*.fastq

R1=${OUT}/${SAMPLE}_1.fastq.gz
R2=${OUT}/${SAMPLE}_2.fastq.gz
```

如果是 single-end SRR，把后面的 STAR `--readFilesIn` 改成只输入一个 FASTQ。

### 3.2 准备参考索引

```bash
# FASTA index
samtools faidx ${REF}

# sequence dictionary
gatk CreateSequenceDictionary \
  -R ${REF} \
  -O ${DICT}

# VCF index
gatk IndexFeatureFile -I ${DBSNP}
gatk IndexFeatureFile -I ${KNOWN_INDEL1}
gatk IndexFeatureFile -I ${KNOWN_INDEL2}
```

### 3.3 建 STAR index

GATK WDL 中 STAR index 的关键参数是 `--sjdbGTFfile` 和 `--sjdbOverhang read_length-1`。

```bash
mkdir -p ${STAR_INDEX}

STAR \
  --runThreadN ${THREADS} \
  --runMode genomeGenerate \
  --genomeDir ${STAR_INDEX} \
  --genomeFastaFiles ${REF} \
  --sjdbGTFfile ${GTF} \
  --sjdbOverhang $((READLEN - 1))
```

### 3.4 STAR two-pass 比对，加 read group

GATK 必须有 read group，尤其是 `SM`。没有 `SM` 会导致后续 GATK/Picard 报错。

```bash
STAR \
  --runThreadN ${THREADS} \
  --genomeDir ${STAR_INDEX} \
  --readFilesIn ${R1} ${R2} \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 45000000000 \
  --outSAMattrRGline ID:${SAMPLE}.rg1 SM:${SAMPLE} LB:${SAMPLE}.lib1 PL:ILLUMINA PU:${SAMPLE}.unit1 \
  --outFileNamePrefix ${OUT}/${SAMPLE}.

samtools index ${OUT}/${SAMPLE}.Aligned.sortedByCoord.out.bam
```

WDL 中 STAR align 也用了 two-pass、sorted BAM、`limitBAMsortRAM` 和 `limitOutSJcollapsed`。

### 3.5 MarkDuplicates

```bash
gatk --java-options "-Xmx8g" MarkDuplicates \
  -I ${OUT}/${SAMPLE}.Aligned.sortedByCoord.out.bam \
  -O ${OUT}/${SAMPLE}.dedup.bam \
  -M ${OUT}/${SAMPLE}.markdup.metrics.txt \
  --CREATE_INDEX true \
  --VALIDATION_STRINGENCY SILENT
```

WDL 对普通 RNA-seq 用 `MarkDuplicates --CREATE_INDEX true`；如果是 UMI RNA-seq，不建议用普通 MarkDuplicates 作为最终去重复逻辑，应参考 WARP 里的 UMI-aware duplicate marking。 

### 3.6 SplitNCigarReads：RNA-seq 必跑

```bash
gatk --java-options "-Xmx8g" SplitNCigarReads \
  -R ${REF} \
  -I ${OUT}/${SAMPLE}.dedup.bam \
  -O ${OUT}/${SAMPLE}.split.bam
```

WDL 的 `SplitNCigarReads` 核心就是 `-R ref -I input.bam -O split.bam`。

### 3.7 BQSR

```bash
gatk --java-options "-Xmx8g" BaseRecalibrator \
  -R ${REF} \
  -I ${OUT}/${SAMPLE}.split.bam \
  --use-original-qualities \
  --known-sites ${DBSNP} \
  --known-sites ${KNOWN_INDEL1} \
  --known-sites ${KNOWN_INDEL2} \
  -O ${OUT}/${SAMPLE}.recal.table

gatk --java-options "-Xmx8g" ApplyBQSR \
  -R ${REF} \
  -I ${OUT}/${SAMPLE}.split.bam \
  --use-original-qualities \
  --bqsr-recal-file ${OUT}/${SAMPLE}.recal.table \
  -O ${OUT}/${SAMPLE}.recal.bam

samtools index ${OUT}/${SAMPLE}.recal.bam
```

GATK WDL 在 `BaseRecalibrator` 和 `ApplyBQSR` 中都用了 `--use-original-qualities`，并输入 dbSNP 和 known indels。

### 3.8 从 GTF 生成 exon interval

GATK WDL 先取 GTF 中 `exon`，转 BED，再用 `BedToIntervalList` 生成 calling interval。

```bash
awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {print $1, $4-1, $5}' ${GTF} \
  | sort -k1,1 -k2,2n \
  > ${OUT}/exons.bed

gatk BedToIntervalList \
  -I ${OUT}/exons.bed \
  -O ${OUT}/exons.interval_list \
  -SD ${DICT}
```

### 3.9 HaplotypeCaller

```bash
gatk --java-options "-Xmx8g" HaplotypeCaller \
  -R ${REF} \
  -I ${OUT}/${SAMPLE}.recal.bam \
  -L ${OUT}/exons.interval_list \
  -O ${OUT}/${SAMPLE}.raw.vcf.gz \
  --dont-use-soft-clipped-bases \
  --standard-min-confidence-threshold-for-calling 20 \
  --dbsnp ${DBSNP} \
  --native-pair-hmm-threads 4
```

WDL 中 HaplotypeCaller 的核心参数就是 `-R`、`-I`、`-L interval_list`、`-O`、`-dont-use-soft-clipped-bases`、`--standard-min-confidence-threshold-for-calling 20` 和 `--dbsnp`。 GATK 文档也说明 `--dont-use-soft-clipped-bases` 的作用是“不分析 reads 中 soft clipped bases”，`--standard-min-confidence-threshold-for-calling` 是 variant calling 的 phred confidence 阈值。([GATK][3]) ([GATK][3])

### 3.10 VariantFiltration

```bash
gatk --java-options "-Xmx4g" VariantFiltration \
  -R ${REF} \
  -V ${OUT}/${SAMPLE}.raw.vcf.gz \
  --window 35 \
  --cluster 3 \
  --filter-name "FS" \
  --filter "FS > 30.0" \
  --filter-name "QD" \
  --filter "QD < 2.0" \
  -O ${OUT}/${SAMPLE}.filtered.vcf.gz
```

这和 GATK WDL 的硬过滤一致：`--window 35`、`--cluster 3`、`FS > 30.0`、`QD < 2.0`。

提取 PASS 变异：

```bash
bcftools view -f PASS -Oz \
  -o ${OUT}/${SAMPLE}.PASS.vcf.gz \
  ${OUT}/${SAMPLE}.filtered.vcf.gz

tabix -p vcf ${OUT}/${SAMPLE}.PASS.vcf.gz
```

## 4. nf-core/rnavar 一键运行版本

样本表：

```bash
cat > samplesheet.csv <<EOF
sample,fastq_1,fastq_2
${SAMPLE},${R1},${R2}
EOF
```

用内置 GRCh38：

```bash
nextflow run nf-core/rnavar -r 1.3.0 \
  -profile docker \
  --input samplesheet.csv \
  --outdir results_rnavar \
  --genome GRCh38
```

如果用自定义参考：

```bash
nextflow run nf-core/rnavar -r 1.3.0 \
  -profile docker \
  --input samplesheet.csv \
  --outdir results_rnavar \
  --igenomes_ignore \
  --fasta ${REF} \
  --gtf ${GTF} \
  --star_index ${STAR_INDEX} \
  --dbsnp ${DBSNP} \
  --known_indels "${KNOWN_INDEL1},${KNOWN_INDEL2}" \
  --read_length ${READLEN} \
  --gatk_hc_call_conf 20 \
  --gatk_interval_scatter_count 25 \
  --gatk_vf_qd_filter 2 \
  --gatk_vf_fs_filter 30 \
  --gatk_vf_window_size 35 \
  --gatk_vf_cluster_size 3
```

nf-core/rnavar 的参数 schema 显示：`--fasta`、`--gtf`、`--star_index`、`--exon_bed`、`--read_length`、`--known_indels`、`--dbsnp` 都是核心参考/资源参数；如果没有 `--genome`，`--fasta` 和 `--gtf` 就是必需参考输入。 

## 5. 给 PrecisionProDB 的建议用法

对 proteogenomics，不建议盲目把 RNA-only VCF 全部放入数据库。更稳妥是：

1. **有 WGS/WES VCF 时优先用 DNA VCF**，RNA-seq 只用来支持表达、剪接和 allele-specific expression。
2. **只有 RNA-seq 时**，用 `PASS.vcf.gz`，再加 DP/AD/VAF 门槛，例如 `DP >= 10`、ALT reads ≥ 3、VAF ≥ 0.1 或更严格。
3. **保留 raw/filtered 两版**：

   * `raw.vcf.gz` 用于追踪召回；
   * `PASS.vcf.gz` 用于 PrecisionProDB 主数据库；
   * 低置信度变异单独做敏感性分析。
4. **参考必须一致**：`REF/GTF/protein FASTA/VCF contig names` 必须同一体系；PrecisionProDB v2 支持 RefSeq、GENCODE 和自定义 GTF，所以 RNA-seq calling 的 GTF 最好和后续蛋白数据库构建的 GTF 保持一致。PrecisionProDB v2 草稿中也强调了它可支持 RefSeq、GENCODE 和 custom GTF annotation。

[1]: https://nf-co.re/rnavar/1.3.0 "rnavar: Introduction"
[2]: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels "RNAseq short variant discovery (SNPs + Indels) – GATK"
[3]: https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller "HaplotypeCaller – GATK"


