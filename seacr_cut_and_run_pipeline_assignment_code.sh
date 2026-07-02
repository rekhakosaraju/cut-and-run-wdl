#!/bin/bash

cat <<'BANNER'

# ------------------------------------------------------------------------------
# SEACR CUT&RUN PIPELINE WITH E. COLI SPIKE‑IN SCALING (PE SEQUENCING ONLY)
# ------------------------------------------------------------------------------
# This script processes CUT&RUN paired‑end FASTQs through trimming, alignment,
# duplicate removal, fragment filtering, peak calling, BigWig generation,
# gene‑feature annotation, and reporting. This version includes alignment of
# trimmed reads to the standard E. coli genome for normalization and scaling.
# ------------------------------------------------------------------------------

################################################################################
#              STEPS                                                           #
################################################################################
#  ---- Preparatory steps  
#  0. Load parameters from config.json                                         #
#  1. Define paths to tools                                                    #
#  2. Create output directories                                                #
#  3. Define utility functions                                                 #
#  4. Derive trimmed filenames                                                 #
#  5. Set genome size                                                          #                                                    # 
#                                                                           #
#  ---- Computational steps                                                    #
#  6.  FastQC — raw read quality control                                       #
#  7.  Trim Galore — adapter & quality trimming                                #
#  8.  Control sample logic                                                    #
#  9.  STAR spike-in alignment (E. coli)                                       #
#  10. STAR host genome alignment                                              #
#  11. Picard AddRG + MarkDuplicates                                           #
#  12. Fragment-length filtering                                               #
#  13. SEACR peak calling                                                      #
#  14. Spike-in scale factor calculation                                       #
#  15. BigWig generation (with spike-in scaling)                               #
#  16. Peak annotation with bedtools intersect                                 #
#                                                                              #
################################################################################

BANNER

set -o pipefail

# ------------------------------------------------------------------------------
# 0  Load parameters from config.json
# ------------------------------------------------------------------------------
CONFIG_FILE="/path/to/config.json"

RAW_FASTQ_DIR=$(jq -r '.raw_fastq_dir'              "$CONFIG_FILE")
ALIGNMENT_DIR=$(jq -r '.alignment_dir'              "$CONFIG_FILE")
OUTPUT_DIR=$(jq  -r '.output_dir'                   "$CONFIG_FILE")
LOG_DIR=$(jq     -r '.log_dir'                      "$CONFIG_FILE")

TREATMENT_R1=$(jq -r '.samples.treatment.r1'        "$CONFIG_FILE")
TREATMENT_R2=$(jq -r '.samples.treatment.r2'        "$CONFIG_FILE")
CONTROL_R1=$(  jq -r '.samples.control.r1 // empty' "$CONFIG_FILE")
CONTROL_R2=$(  jq -r '.samples.control.r2 // empty' "$CONFIG_FILE")

REFERENCE_GENOME=$(jq -r '.reference_genome'        "$CONFIG_FILE")
ANNOTATION_GENES=$(jq -r '.annotation_genes'        "$CONFIG_FILE")
CHROM_SIZE=$(      jq -r '.chrom_sizes'             "$CONFIG_FILE")
ECOLI_INDEX=$(     jq -r '.ecoli_index'             "$CONFIG_FILE")

# Dbug paths
echo "DEBUG paths: T1=$TREATMENT_R1  T2=$TREATMENT_R2  C1=$CONTROL_R1  C2=$CONTROL_R2"

# Optional explicit spike‑in BAM directory; default → <alignment>/spikein
SPIKE_BAM_DIR=$(jq -r '.spike_bam_dir // empty'     "$CONFIG_FILE")
[ -z "$SPIKE_BAM_DIR" ] || [ "$SPIKE_BAM_DIR" = "null" ] && \
  SPIKE_BAM_DIR="$ALIGNMENT_DIR/spikein"

GENOME_SIZE_STRING=$(jq -r '.genome_size'            "$CONFIG_FILE")
FRAGMENT_SIZE_FILTER=$(jq -r '.fragment_size_filter' "$CONFIG_FILE")
CUSTOM_GENOME_SIZE=$( jq -r '.custom_genome_size'    "$CONFIG_FILE")
NUM_THREADS=$(       jq -r '.num_threads'            "$CONFIG_FILE")

# When read count is low or fragment length is narrow (as expected in histone or TF targeting experiments), MACS2 simply falls back to non-model-based peak calling and recommends these extsizes. Used in SE only.
BROAD_EXTSIZE=$( jq -r '.broad_peak_extsize'  "$CONFIG_FILE")
NARROW_EXTSIZE=$(jq -r '.narrow_peak_extsize' "$CONFIG_FILE")

# ------------------------------------------------------------------------------
# Build FASTQ_FILES array 
# ------------------------------------------------------------------------------
FASTQ_FILES=("$TREATMENT_R1" "$TREATMENT_R2")

# Only push control files if they exist
if [[ -n "$CONTROL_R1" && -n "$CONTROL_R2" ]]; then
  FASTQ_FILES+=("$CONTROL_R1" "$CONTROL_R2")
fi

# ------------------------------------------------------------------------------
# 1 Paths to tools and software
# ------------------------------------------------------------------------------
# Picard tools path
PICARD_PATH="/path/to/picard.jar"  # Path to the Picard jar file (e.g., picard.jar)
# FastQC tools path
FASTQC_PATH="/path/to/FastQC/fastqc"
# STAR path
STAR_PATH="/path/to/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR"
# SEACR path
SEACR_SH="$SEACR_DIR/SEACR_1.3.sh"              # shell wrapper

# -----------------------------------------------------------------------------
# 2  Create required directories
# -----------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$ALIGNMENT_DIR" "$SPIKE_BAM_DIR"
FASTQC_DIR="$OUTPUT_DIR/fastqc_reports";    mkdir -p "$FASTQC_DIR"
PEAK_DIR="$OUTPUT_DIR/macs2_peaks";         mkdir -p "$PEAK_DIR"
BW_DIR="$OUTPUT_DIR/bigwig_bedgraphs";      mkdir -p "$BW_DIR"
ANN_DIR="$OUTPUT_DIR/annotated_peaks";      mkdir -p "$ANN_DIR"
MULTIQC_DIR="$OUTPUT_DIR/multiqc_reports";  mkdir -p "$MULTIQC_DIR"

# ------------------------------------------------------------------------------
# 3  Utility functions
# ------------------------------------------------------------------------------
get_sample_basename() {
  local f=$1
  local b=${f##*/}            # strip directory
  b=${b%.fastq.gz}; b=${b%.fq.gz}
  b=${b%_R1_001};  b=${b%_R2_001}   # ← strip the full _R#_001 suffix
  b=${b%_R1};      b=${b%_R2}       # fallback for plain _R1/_R2
  b=${b%_1};       b=${b%_2}        # fallback for _1/_2
  echo "$b" | sed 's/[^A-Za-z0-9._-]//g'
}

# Usage: run_star <STAR_index> <R1> <R2> <outPrefix> <logBase>
# Runs STAR → sorted BAM + logs.
run_star() {
  local index=$1 r1=$2 r2=$3 prefix=$4 logbase=$5
  "$STAR_PATH" --runThreadN "$NUM_THREADS" \
               --genomeDir  "$REFERENCE_GENOME" \
               --readFilesIn "$r1" "$r2" \
               --readFilesCommand zcat \
               --outSAMtype BAM SortedByCoordinate \
               --outFileNamePrefix "$prefix" \
               > "$LOG_DIR/${logbase}.log"  \
               2> "$LOG_DIR/${logbase}_err.log"
}

run_spikein_align() {
  # Align trimmed reads to E. coli spike‑in genome with STAR.
  # Usage: run_spikein_align R1 R2 SAMPLE_BASENAME
  local r1=$1 r2=$2 sample=$3
  echo "[SPIKE‑IN] Aligning $sample to E. coli genome…" | tee -a "$LOG_DIR/pipeline.log"
  "$STAR_PATH" --runThreadN "$NUM_THREADS" \
               --genomeDir  "$ECOLI_INDEX" \
               --readFilesIn "$r1" "$r2" \
               --readFilesCommand zcat \
               --outSAMtype BAM SortedByCoordinate \
               --outFileNamePrefix "$SPIKE_BAM_DIR/${sample}_ecoli_" \
               > "$LOG_DIR/STAR_${sample}_ecoli.log"  \
               2> "$LOG_DIR/STAR_${sample}_ecoli_error.log"

  local tmp="$SPIKE_BAM_DIR/${sample}_ecoli_Aligned.sortedByCoord.out.bam"
  if [[ -f "$tmp" ]]; then
    mv "$tmp" "$SPIKE_BAM_DIR/${sample}.ecoli.sorted.bam"
    samtools index "$SPIKE_BAM_DIR/${sample}.ecoli.sorted.bam"
  else
    echo "❌ Spike‑in alignment failed for $sample — BAM not produced." | tee -a "$LOG_DIR/pipeline.log"
    exit 1
  fi
}

bam_to_bedgraph () {
  local inbam=$1 outbg=$2
  # fragment-level coverage (-pc) → BedGraph sorted by chrom,start
  bedtools genomecov -ibam "$inbam" -bg -pc | \
    sort -k1,1 -k2,2n > "$outbg"
}

# ------------------------------------------------------------------------------
# 4  Derive filenames for downstream steps
# ------------------------------------------------------------------------------
TREATMENT_BASE=$(get_sample_basename "$TREATMENT_R1")
TREATMENT_TRIMMED_R1="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R1.fq.gz"
TREATMENT_TRIMMED_R2="$ALIGNMENT_DIR/${TREATMENT_BASE}_trimmed_R2.fq.gz"

if [[ -n "$CONTROL_R1" ]]; then
  CONTROL_BASE=$(get_sample_basename "$CONTROL_R1")
  CONTROL_TRIMMED_R1="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz"
  CONTROL_TRIMMED_R2="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz"
fi

# ------------------------------------------------------------------------------
# 5  Compute numeric genome size
# ------------------------------------------------------------------------------
GENOME_SIZE_HUMAN=2913022398; GENOME_SIZE_MOUSE=2652783500
GENOME_SIZE_DROSOPHILA=165000000; GENOME_SIZE_CELEGANS=1000000000
GENOME_SIZE_YEAST=12000000

case "$GENOME_SIZE_STRING" in
  hs) GENOME_SIZE=$GENOME_SIZE_HUMAN ;;
  mm) GENOME_SIZE=$GENOME_SIZE_MOUSE ;;
  dm) GENOME_SIZE=$GENOME_SIZE_DROSOPHILA ;;
  ce) GENOME_SIZE=$GENOME_SIZE_CELEGANS ;;
  sc) GENOME_SIZE=$GENOME_SIZE_YEAST ;;
  *)  if [[ "$CUSTOM_GENOME_SIZE" != "null" ]]; then
        GENOME_SIZE=$CUSTOM_GENOME_SIZE
      else
        echo "Error: Invalid genome size string $GENOME_SIZE_STRING" | tee -a "$LOG_DIR/pipeline.log"
        exit 1
      fi ;;
esac

echo "Using genome size $GENOME_SIZE for host genome: $GENOME_SIZE_STRING" | tee -a "$LOG_DIR/pipeline.log"

# ------------------------------------------------------------------------------
# 6  FASTQC (raw reads)
# ------------------------------------------------------------------------------
FASTQ_FILES=("$TREATMENT_R1" "$TREATMENT_R2")
if [[ -n "$CONTROL_R1" ]]; then FASTQ_FILES+=("$CONTROL_R1" "$CONTROL_R2"); fi

echo "Running FastQC…" | tee -a "$LOG_DIR/pipeline.log"
for fq in "${FASTQ_FILES[@]}"; do
  $FASTQC_PATH --extract -o "$FASTQC_DIR" "$fq" >> "$LOG_DIR/pipeline.log" 2>&1
done

# ------------------------------------------------------------------------------
# 7  Adapter trimming (Trim Galore!)  –  trim ALL declared FASTQ pairs
# ------------------------------------------------------------------------------
echo "[Trim Galore] starting…" | tee -a "$LOG_DIR/pipeline.log"
echo "DEBUG: FASTQ_FILES = ${FASTQ_FILES[*]}" | tee -a "$LOG_DIR/pipeline.log"

i=0
while [[ $i -lt ${#FASTQ_FILES[@]} ]]; do
  R1=${FASTQ_FILES[$i]}
  R2=${FASTQ_FILES[$((i+1))]}
  BASE=$(get_sample_basename "$R1")

  echo "  ↳ trimming $BASE" | tee -a "$LOG_DIR/pipeline.log"
  trim_galore --paired --quality 20 --phred33 \
            --output_dir "$ALIGNMENT_DIR" "$R1" "$R2" \
            > "$LOG_DIR/trim_${BASE}.log" 2>&1

  VAL1=$(find "$ALIGNMENT_DIR" -name "*_val_1.fq.gz" | grep "$BASE" | head -n1)
  VAL2=$(find "$ALIGNMENT_DIR" -name "*_val_2.fq.gz" | grep "$BASE" | head -n1)

  if [[ -f "$VAL1" && -f "$VAL2" ]]; then
    mv "$VAL1" "$ALIGNMENT_DIR/${BASE}_trimmed_R1.fq.gz"
    mv "$VAL2" "$ALIGNMENT_DIR/${BASE}_trimmed_R2.fq.gz"
  else
    echo "❌ Trim Galore did not produce trimmed files for $BASE — skipping." | tee -a "$LOG_DIR/pipeline.log"
  fi
  i=$((i+2))
done

# DEBUGGING
echo "DEBUG: FASTQ_FILES = ${FASTQ_FILES[*]}"
echo "DEBUG: R1=$R1  R2=$R2  BASE=$(get_sample_basename "$R1")" \
  | tee -a "$LOG_DIR/pipeline.log"

# ------------------------------------------------------------------------------
# 8 Include/exclude control sample logic (for downstream steps)
# ------------------------------------------------------------------------------
USE_CONTROL=0
if [[ -n "${CONTROL_BASE:-}" ]] && \
   [[ -f "$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz" ]] && \
   [[ -f "$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz" ]]; then
  USE_CONTROL=1
  CONTROL_TRIMMED_R1="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R1.fq.gz"
  CONTROL_TRIMMED_R2="$ALIGNMENT_DIR/${CONTROL_BASE}_trimmed_R2.fq.gz"
else
  echo "⚠️  No trimmed control FASTQs found — proceeding without control." | tee -a "$LOG_DIR/pipeline.log"
fi

# Build sample list for the rest of the script
SAMPLES=("$TREATMENT_BASE")
[[ $USE_CONTROL -eq 1 ]] && SAMPLES+=("$CONTROL_BASE")

# ------------------------------------------------------------------------------
# 9  Spike-in alignment (E. coli)
# ------------------------------------------------------------------------------
echo "Aligning to the E. coli genome with STAR for subsequent spike-in scaling…" | tee -a "$LOG_DIR/pipeline.log"
run_spikein_align "$TREATMENT_TRIMMED_R1" "$TREATMENT_TRIMMED_R2" "$TREATMENT_BASE"

if [[ $USE_CONTROL -eq 1 ]]; then
  run_spikein_align "$CONTROL_TRIMMED_R1" "$CONTROL_TRIMMED_R2" "$CONTROL_BASE"
fi

# ------------------------------------------------------------------------------
# 10  Host‑genome alignment (STAR)
# ------------------------------------------------------------------------------
echo "Aligning to the host genome with STAR…" | tee -a "$LOG_DIR/pipeline.log"
run_star "$REFERENCE_GENOME" \
         "$TREATMENT_TRIMMED_R1" "$TREATMENT_TRIMMED_R2" \
         "$ALIGNMENT_DIR/${TREATMENT_BASE}." "STAR_${TREATMENT_BASE}"

if [[ $USE_CONTROL -eq 1 ]]; then
  run_star "$REFERENCE_GENOME" \
           "$CONTROL_TRIMMED_R1" "$CONTROL_TRIMMED_R2" \
           "$ALIGNMENT_DIR/${CONTROL_BASE}." "STAR_${CONTROL_BASE}"
fi

# ------------------------------------------------------------------------------
# 11  Picard AddRG + MarkDuplicates          (only current-run BAMs)            
# ------------------------------------------------------------------------------
echo "[Picard] processing ${SAMPLES[*]}" | tee -a "$LOG_DIR/pipeline.log"

for samp in "${SAMPLES[@]}"; do
  in_bam="$ALIGNMENT_DIR/${samp}.Aligned.sortedByCoord.out.bam"
  [[ -s "$in_bam" ]] || { echo "❌ BAM not found: $in_bam" | tee -a "$LOG_DIR/pipeline.log"; exit 1; }

  java -jar "$PICARD_PATH" AddOrReplaceReadGroups \
       I="$in_bam" \
       O="$ALIGNMENT_DIR/${samp}.rg.bam" \
       RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM="$samp" \
       VALIDATION_STRINGENCY=LENIENT

  java -jar "$PICARD_PATH" MarkDuplicates \
       I="$ALIGNMENT_DIR/${samp}.rg.bam" \
       O="$ALIGNMENT_DIR/${samp}.dedup.bam" \
       M="$LOG_DIR/${samp}.metrics.txt" \
       REMOVE_DUPLICATES=true \
       VALIDATION_STRINGENCY=LENIENT

  samtools index "$ALIGNMENT_DIR/${samp}.dedup.bam"
done

# ------------------------------------------------------------------------------
# 12  Fragment‑length filtering
# ------------------------------------------------------------------------------
case "$FRAGMENT_SIZE_FILTER" in
  histones)              FRAG_CMD='{if ($9 >= 130 && $9 <= 300 || $1 ~ /^@/) print $0}' ;;
  transcription_factors) FRAG_CMD='{if ($9 < 130 || $1 ~ /^@/) print $0}' ;;
  *)                     FRAG_CMD='{if ($9 < 1000 || $1 ~ /^@/) print $0}' ;;
esac

echo "Filtering fragments by range $FRAGMENT_SIZE_FILTER…" | tee -a "$LOG_DIR/pipeline.log"
for bam in "$ALIGNMENT_DIR"/*.dedup.bam; do
  base=$(basename "$bam" .dedup.bam)
  samtools view -h "$bam" | awk "$FRAG_CMD" | samtools view -bS - > "$ALIGNMENT_DIR/${base}.dedup.filtered.bam"
done

# ------------------------------------------------------------------------------
# 13  Peak calling (SEACR)
# ------------------------------------------------------------------------------
echo "[SEACR] calling peaks" | tee -a "$LOG_DIR/pipeline.log"

# ---- make BedGraph for treatment --------------------------------------------
TREAT_BG="$PEAK_DIR/${TREATMENT_BASE}.bedgraph"
bam_to_bedgraph "$ALIGNMENT_DIR/${TREATMENT_BASE}.dedup.filtered.bam" "$TREAT_BG"

# ---- if control exists, make BedGraph for control ---------------------------
if [[ $USE_CONTROL -eq 1 ]]; then
  CONTROL_BG="$PEAK_DIR/${CONTROL_BASE}.bedgraph"
  bam_to_bedgraph "$ALIGNMENT_DIR/${CONTROL_BASE}.dedup.filtered.bam" "$CONTROL_BG"

  # ---------- SEACR control-subtracted peaks (two thresholds) ----------------
  # stringent (FDR-style 0.01)
  bash "$SEACR_SH" "$TREAT_BG" "$CONTROL_BG" non stringent \
       "$PEAK_DIR/${TREATMENT_BASE}_seacr_ctrl_0.01"

  # relaxed (0.05)
  bash "$SEACR_SH" "$TREAT_BG" "$CONTROL_BG" non 0.05 \
       "$PEAK_DIR/${TREATMENT_BASE}_seacr_ctrl_0.05"

else
  # ---------- SEACR with no control (background = treatment) -----------------
  bash "$SEACR_SH" "$TREAT_BG" 0.01 non stringent \
       "$PEAK_DIR/${TREATMENT_BASE}_seacr_bg_0.01"

  bash "$SEACR_SH" "$TREAT_BG" 0.05 non stringent \
       "$PEAK_DIR/${TREATMENT_BASE}_seacr_bg_0.05"
fi

# ---- convert SEACR BEDs → narrowPeak (6-column dummy) -----------------------
for bed in "$PEAK_DIR"/${TREATMENT_BASE}_seacr_*".bed"; do
  [[ -f "$bed" ]] || continue
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",1000,"."}' "$bed" \
    > "${bed%.bed}.narrowPeak"
done


# ------------------------------------------------------------------------------
# 14  Spike-in scaling factors (per-run samples only)                         
# ------------------------------------------------------------------------------
echo "[Spike-in] calculating scale factors" | tee -a "$LOG_DIR/pipeline.log"

declare -A SCALE   # sample → factor

count_reads() { samtools view -c -F 2304 "$1"; }

for samp in "${SAMPLES[@]}"; do
  host_bam="$ALIGNMENT_DIR/${samp}.dedup.filtered.bam"
  ecoli_bam="$SPIKE_BAM_DIR/${samp}.ecoli.sorted.bam"

  if [[ -f "$host_bam" && -f "$ecoli_bam" ]]; then
    host_reads=$(count_reads "$host_bam")
    spike_reads=$(count_reads "$ecoli_bam")

    if (( spike_reads > 0 )); then
      factor=$(awk -v h=$host_reads -v s=$spike_reads 'BEGIN{printf "%.6f", h/s}')
      echo "  ↳ $samp : host=$host_reads  spike=$spike_reads  scale=$factor" \
           | tee -a "$LOG_DIR/pipeline.log"
      SCALE["$samp"]=$factor
    else
      echo "  ↳ $samp : spike reads = 0 — no scaling applied" \
           | tee -a "$LOG_DIR/pipeline.log"
    fi
  else
    echo "  ↳ $samp : missing host or spike BAM — skipped" \
         | tee -a "$LOG_DIR/pipeline.log"
  fi
done

# ------------------------------------------------------------------------------
# 15  BigWig generation (with optional scaling)                                
# ------------------------------------------------------------------------------
echo "[BigWig] generating coverage tracks" | tee -a "$LOG_DIR/pipeline.log"

for samp in "${SAMPLES[@]}"; do
  host_bam="$ALIGNMENT_DIR/${samp}.dedup.filtered.bam"
  [[ -f "$host_bam" ]] || { echo "❌ BAM missing for $samp — skipping" | tee -a "$LOG_DIR/pipeline.log"; continue; }

  scale_opt=""
  if [[ -n "${SCALE[$samp]:-}" ]]; then
    echo "  ↳ $samp : applying scaleFactor ${SCALE[$samp]}" | tee -a "$LOG_DIR/pipeline.log"
    # bedtools genomecov lacks a scale flag, so we multiply depth via awk
    scaled_bg="$BW_DIR/${samp}.scaled.bedgraph"
    bedtools genomecov -ibam "$host_bam" -bg -pc | \
      awk -v f="${SCALE[$samp]}" '{ $4=$4*f; print }' > "$scaled_bg"
    bedGraphToBigWig "$scaled_bg" "$CHROM_SIZE" "$BW_DIR/${samp}.bw"
  else
    echo "  ↳ $samp : no scaleFactor (spike-in absent)" | tee -a "$LOG_DIR/pipeline.log"
    bedtools genomecov -ibam "$host_bam" -bg -pc > "$BW_DIR/${samp}.bedgraph"
    bedGraphToBigWig "$BW_DIR/${samp}.bedgraph" "$CHROM_SIZE" "$BW_DIR/${samp}.bw"
  fi
done

# ------------------------------------------------------------------------------
# 16  Peak annotation
# ------------------------------------------------------------------------------
echo "[Peak annotation] intersecting peaks with gene features" | tee -a "$LOG_DIR/pipeline.log"

for samp in "${SAMPLES[@]}"; do
  peak_file="$PEAK_DIR/${TREATMENT_BASE}_seacr_0.05.narrowPeak"
  [[ -f "$peak_file" ]] || { echo "  ↳ $samp : no narrowPeak file — skipping" | tee -a "$LOG_DIR/pipeline.log"; continue; }

  full_out="$ANN_DIR/${samp}_peaks.annotated.bed"
  tsv_out="$ANN_DIR/${samp}_peaks.annotated.tsv"

  echo "  ↳ annotating $samp" | tee -a "$LOG_DIR/pipeline.log"

  # full BED12 style intersect (peak + gene feature columns)
  bedtools intersect -a "$peak_file" -b "$ANNOTATION_GENES" -wa -wb > "$full_out"

  # concise TSV: peak coords + gene name + strand
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$10,$11,$12}' "$full_out" > "$tsv_out"
done

echo "🎉  CUT&RUN pipeline complete!  Results are in: $OUTPUT_DIR" | tee -a "$LOG_DIR/pipeline.log"

