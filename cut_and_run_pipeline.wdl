version 1.0

## =============================================================================
## CUT&RUN Pipeline WDL Workflow
## Wraps the SEACR CUT&RUN bash pipeline for execution on AWS via WDL/Cromwell
## or AWS HealthOmics.
## 
## Inputs:  Paired-end FASTQ files (treatment + optional control),
##          reference genome index, E. coli spike-in index, annotation files
## Outputs: Peaks (BED/narrowPeak), BigWig tracks, annotated peak TSVs
## =============================================================================

workflow CutAndRunPipeline {

  input {
    # ---- Sample FASTQs ----
    File treatment_r1          # Treatment sample, read 1
    File treatment_r2          # Treatment sample, read 2
    File? control_r1           # Optional IgG/control sample, read 1
    File? control_r2           # Optional IgG/control sample, read 2

    # ---- Reference files ----
    File reference_genome_tar  # Tarball of STAR genome index (host)
    File ecoli_genome_tar      # Tarball of STAR E. coli spike-in index
    File chrom_sizes           # Chromosome sizes file for BigWig generation
    File annotation_genes      # Gene annotation BED file for peak annotation

    # ---- Parameters ----
    String genome_size_string          # hs, mm, dm, ce, sc, or custom
    String fragment_size_filter        # histones | transcription_factors | all
    Int    num_threads         = 8
    String docker_image        = "your-ecr-repo/cutandrun:latest"

    # ---- Runtime ----
    Int    memory_gb           = 32
    Int    disk_gb             = 200
  }

  # ---------------------------------------------------------------------------
  # Step 1: FastQC — quality check raw reads
  # ---------------------------------------------------------------------------
  call FastQC {
    input:
      r1            = treatment_r1,
      r2            = treatment_r2,
      docker_image  = docker_image,
      memory_gb     = 8,
      disk_gb       = 50
  }

  # ---------------------------------------------------------------------------
  # Step 2: Trim Galore — remove adapters and low-quality bases
  # ---------------------------------------------------------------------------
  call TrimGalore {
    input:
      r1           = treatment_r1,
      r2           = treatment_r2,
      docker_image = docker_image,
      memory_gb    = 8,
      disk_gb      = 100
  }

  # Optional: trim control reads if provided
  if (defined(control_r1) && defined(control_r2)) {
    call TrimGalore as TrimControl {
      input:
        r1           = select_first([control_r1]),
        r2           = select_first([control_r2]),
        docker_image = docker_image,
        memory_gb    = 8,
        disk_gb      = 100
    }
  }

  # ---------------------------------------------------------------------------
  # Step 3: Spike-in alignment — align to E. coli for normalization
  # ---------------------------------------------------------------------------
  call SpikeInAlign {
    input:
      r1              = TrimGalore.trimmed_r1,
      r2              = TrimGalore.trimmed_r2,
      ecoli_index_tar = ecoli_genome_tar,
      num_threads     = num_threads,
      docker_image    = docker_image,
      memory_gb       = memory_gb,
      disk_gb         = disk_gb
  }

  # ---------------------------------------------------------------------------
  # Step 4: Host genome alignment — align to reference genome with STAR
  # ---------------------------------------------------------------------------
  call HostAlign {
    input:
      r1                   = TrimGalore.trimmed_r1,
      r2                   = TrimGalore.trimmed_r2,
      reference_genome_tar = reference_genome_tar,
      num_threads          = num_threads,
      docker_image         = docker_image,
      memory_gb            = memory_gb,
      disk_gb              = disk_gb
  }

  # ---------------------------------------------------------------------------
  # Step 5: Picard — add read groups and remove duplicates
  # ---------------------------------------------------------------------------
  call MarkDuplicates {
    input:
      input_bam    = HostAlign.aligned_bam,
      sample_name  = "treatment",
      docker_image = docker_image,
      memory_gb    = 16,
      disk_gb      = 100
  }

  # ---------------------------------------------------------------------------
  # Step 6: Fragment size filtering
  # ---------------------------------------------------------------------------
  call FilterFragments {
    input:
      input_bam           = MarkDuplicates.dedup_bam,
      fragment_size_filter = fragment_size_filter,
      docker_image        = docker_image,
      memory_gb           = 8,
      disk_gb             = 100
  }

  # ---------------------------------------------------------------------------
  # Step 7: SEACR peak calling
  # ---------------------------------------------------------------------------
  call CallPeaks {
    input:
      treatment_bam = FilterFragments.filtered_bam,
      docker_image  = docker_image,
      memory_gb     = 16,
      disk_gb       = 50
  }

  # ---------------------------------------------------------------------------
  # Step 8: Compute spike-in scale factor and generate BigWig
  # ---------------------------------------------------------------------------
  call GenerateBigWig {
    input:
      filtered_bam  = FilterFragments.filtered_bam,
      ecoli_bam     = SpikeInAlign.ecoli_bam,
      chrom_sizes   = chrom_sizes,
      sample_name   = "treatment",
      docker_image  = docker_image,
      memory_gb     = 8,
      disk_gb       = 50
  }

  # ---------------------------------------------------------------------------
  # Step 9: Peak annotation
  # ---------------------------------------------------------------------------
  call AnnotatePeaks {
    input:
      peaks_bed        = CallPeaks.peaks_stringent,
      annotation_genes = annotation_genes,
      sample_name      = "treatment",
      docker_image     = docker_image,
      memory_gb        = 8,
      disk_gb          = 20
  }

  # ---------------------------------------------------------------------------
  # Workflow outputs
  # ---------------------------------------------------------------------------
  output {
    File fastqc_html           = FastQC.report_html
    File trimmed_r1            = TrimGalore.trimmed_r1
    File trimmed_r2            = TrimGalore.trimmed_r2
    File dedup_bam             = MarkDuplicates.dedup_bam
    File filtered_bam          = FilterFragments.filtered_bam
    File peaks_stringent       = CallPeaks.peaks_stringent
    File peaks_relaxed         = CallPeaks.peaks_relaxed
    File bigwig                = GenerateBigWig.bigwig
    File annotated_peaks_tsv   = AnnotatePeaks.annotated_tsv
  }
}

# =============================================================================
# TASK DEFINITIONS
# =============================================================================

## FastQC: assess raw read quality
task FastQC {
  input {
    File   r1
    File   r2
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    fastqc --extract -o . ~{r1} ~{r2}
  >>>
  output {
    File report_html = glob("*_fastqc.html")[0]
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    2
  }
}

## TrimGalore: remove Illumina adapters and bases with quality < 20
task TrimGalore {
  input {
    File   r1
    File   r2
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    trim_galore --paired --quality 20 --phred33 \
                --output_dir . ~{r1} ~{r2}
    # Rename to consistent output names
    mv *_val_1.fq.gz trimmed_R1.fq.gz
    mv *_val_2.fq.gz trimmed_R2.fq.gz
  >>>
  output {
    File trimmed_r1 = "trimmed_R1.fq.gz"
    File trimmed_r2 = "trimmed_R2.fq.gz"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    4
  }
}

## SpikeInAlign: align trimmed reads to E. coli genome for normalization
task SpikeInAlign {
  input {
    File   r1
    File   r2
    File   ecoli_index_tar
    Int    num_threads
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    # Unpack the E. coli STAR index
    mkdir -p ecoli_index
    tar -xzf ~{ecoli_index_tar} -C ecoli_index

    STAR --runThreadN ~{num_threads} \
         --genomeDir ecoli_index \
         --readFilesIn ~{r1} ~{r2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix spikein_

    mv spikein_Aligned.sortedByCoord.out.bam sample.ecoli.sorted.bam
    samtools index sample.ecoli.sorted.bam
  >>>
  output {
    File ecoli_bam     = "sample.ecoli.sorted.bam"
    File ecoli_bam_bai = "sample.ecoli.sorted.bam.bai"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    num_threads
  }
}

## HostAlign: align trimmed reads to reference genome with STAR
task HostAlign {
  input {
    File   r1
    File   r2
    File   reference_genome_tar
    Int    num_threads
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    mkdir -p ref_index
    tar -xzf ~{reference_genome_tar} -C ref_index

    STAR --runThreadN ~{num_threads} \
         --genomeDir ref_index \
         --readFilesIn ~{r1} ~{r2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix host_
  >>>
  output {
    File aligned_bam = "host_Aligned.sortedByCoord.out.bam"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    num_threads
  }
}

## MarkDuplicates: add read groups and remove PCR duplicates with Picard
task MarkDuplicates {
  input {
    File   input_bam
    String sample_name
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    java -jar /opt/picard.jar AddOrReplaceReadGroups \
         I=~{input_bam} O=rg.bam \
         RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=~{sample_name} \
         VALIDATION_STRINGENCY=LENIENT

    java -jar /opt/picard.jar MarkDuplicates \
         I=rg.bam O=dedup.bam M=metrics.txt \
         REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

    samtools index dedup.bam
  >>>
  output {
    File dedup_bam     = "dedup.bam"
    File dedup_bam_bai = "dedup.bam.bai"
    File metrics       = "metrics.txt"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    2
  }
}

## FilterFragments: keep only fragments in the expected size range
task FilterFragments {
  input {
    File   input_bam
    String fragment_size_filter   # histones | transcription_factors | default
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    case "~{fragment_size_filter}" in
      histones)
        AWK_CMD='{if ($9 >= 130 && $9 <= 300 || $1 ~ /^@/) print $0}';;
      transcription_factors)
        AWK_CMD='{if ($9 < 130 || $1 ~ /^@/) print $0}';;
      *)
        AWK_CMD='{if ($9 < 1000 || $1 ~ /^@/) print $0}';;
    esac

    samtools view -h ~{input_bam} \
      | awk "$AWK_CMD" \
      | samtools view -bS - > filtered.bam
  >>>
  output {
    File filtered_bam = "filtered.bam"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    2
  }
}

## CallPeaks: convert BAM to BedGraph and call peaks with SEACR
task CallPeaks {
  input {
    File   treatment_bam
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    # BAM → sorted BedGraph (fragment-level coverage)
    bedtools genomecov -ibam ~{treatment_bam} -bg -pc \
      | sort -k1,1 -k2,2n > treatment.bedgraph

    # SEACR: stringent peak calling (top 1%)
    bash /opt/SEACR/SEACR_1.3.sh treatment.bedgraph 0.01 non stringent peaks_stringent

    # SEACR: relaxed peak calling (top 5%)
    bash /opt/SEACR/SEACR_1.3.sh treatment.bedgraph 0.05 non stringent peaks_relaxed
  >>>
  output {
    File peaks_stringent = "peaks_stringent.bed"
    File peaks_relaxed   = "peaks_relaxed.bed"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    2
  }
}

## GenerateBigWig: spike-in-scaled coverage track for genome browser visualization
task GenerateBigWig {
  input {
    File   filtered_bam
    File   ecoli_bam
    File   chrom_sizes
    String sample_name
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    # Count host and spike-in reads to compute scale factor
    host_reads=$(samtools view -c -F 2304 ~{filtered_bam})
    spike_reads=$(samtools view -c -F 2304 ~{ecoli_bam})

    if (( spike_reads > 0 )); then
      factor=$(awk -v h=$host_reads -v s=$spike_reads \
               'BEGIN{printf "%.6f", h/s}')
      bedtools genomecov -ibam ~{filtered_bam} -bg -pc \
        | awk -v f="$factor" '{$4=$4*f; print}' > scaled.bedgraph
    else
      bedtools genomecov -ibam ~{filtered_bam} -bg -pc > scaled.bedgraph
    fi

    bedGraphToBigWig scaled.bedgraph ~{chrom_sizes} ~{sample_name}.bw
  >>>
  output {
    File bigwig = "~{sample_name}.bw"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    2
  }
}

## AnnotatePeaks: intersect peaks with gene annotation to assign biological context
task AnnotatePeaks {
  input {
    File   peaks_bed
    File   annotation_genes
    String sample_name
    String docker_image
    Int    memory_gb
    Int    disk_gb
  }
  command <
    bedtools intersect -a ~{peaks_bed} -b ~{annotation_genes} -wa -wb \
      > full_annotation.bed

    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$10,$11,$12}' full_annotation.bed \
      > ~{sample_name}_annotated.tsv
  >>>
  output {
    File annotated_bed = "full_annotation.bed"
    File annotated_tsv = "~{sample_name}_annotated.tsv"
  }
  runtime {
    docker: docker_image
    memory: "~{memory_gb}G"
    disks:  "local-disk ~{disk_gb} SSD"
    cpu:    2
  }
}