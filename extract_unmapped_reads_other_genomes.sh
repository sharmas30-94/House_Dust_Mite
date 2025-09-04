#!/bin/bash

# Allocate a chunk of local disc for the job mentioned in SLURM_JOBID

#########Copy bwa-mem2 references into Job folders#######
export TMPDIR=/lscratch/$SLURM_JOB_ID


module load samtools
module load bwa-mem2/2.2.1

sample=${1}
df_refile="/data/Sherlock_Lung2/Sunandini_HDM/script/references/DF"
dp_refile="/data/Sherlock_Lung2/Sunandini_HDM/script/references/DP"
fastq_files="/data/Sherlock_Lung/JohnMce/Microbiome/Kraken_WGS3"
Sherlock="/data/Sherlock_Lung2/Sunandini_HDM/Sherlock_WGS3"


mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome
mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp
mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_PE
mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_PE/input
mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_PE/output
mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_Singletons
mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_Singletons/input
mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_Singletons/output


#mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_40-49bp
#mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_40-49bp/DF_PE
#mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_40-49bp/DF_PE/input
#mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_40-49bp/DF_PE/output
#mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_40-49bp/DF_Singletons
#mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_40-49bp/DF_Singletons/input
#mkdir -p /data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_40-49bp/DF_Singletons/output


###### Use Dermatophagoides farinae Chromosome  ##########

# Align paired reads using bwa-mem2
bwa-mem2 mem -t $SLURM_CPUS_ON_NODE $df_refile/GCA_024713945.1_ASM2471394v1_genomic.fa \
    $fastq_files/${sample}/${sample}_1.fq.gz $fastq_files/${sample}/${sample}_2.fq.gz > $TMPDIR/${sample}_GCA_024713945.1_DF.sam

# Convert paired SAM to BAM
samtools view -bS $TMPDIR/${sample}_GCA_024713945.1_DF.sam > $TMPDIR/${sample}_GCA_024713945.1_DF.bam

# Sort paired BAM
samtools sort $TMPDIR/${sample}_GCA_024713945.1_DF.bam -o $Sherlock/${sample}/DF_Alignment/${sample}_GCA_024713945.1_DF_sorted.bam

# Extract mapped reads from paired BAM
samtools view -b -F 4 $Sherlock/${sample}/DF_Alignment/${sample}_GCA_024713945.1_DF_sorted.bam > $Sherlock/${sample}/DF_Alignment/${sample}_GCA_024713945.1_DF_mapped_only.bam


# Align paired reads using bwa-mem2
bwa-mem2 mem -t $SLURM_CPUS_ON_NODE $df_refile/GCA_024713945.1_ASM2471394v1_genomic.fa \
    $fastq_files/${sample}/${sample}_single.fq.gz  > $TMPDIR/${sample}_GCA_024713945.1_DF_Singletons.sam

# Convert paired SAM to BAM
samtools view -bS $TMPDIR/${sample}_GCA_024713945.1_DF_Singletons.sam > $TMPDIR/${sample}_GCA_024713945.1_DF_Singletons.bam

# Sort paired BAM
samtools sort $TMPDIR/${sample}_GCA_024713945.1_DF_Singletons.bam -o $Sherlock/${sample}/singletons/${sample}_GCA_024713945.1_DF_singletons_sorted.bam

# Extract mapped reads from paired BAM
samtools view -b -F 4 $Sherlock/${sample}/singletons/${sample}_GCA_024713945.1_DF_singletons_sorted.bam > $Sherlock/${sample}/singletons/${sample}_GCA_024713945.1_DF_singletons_mapped_only.bam



##### Get the cigar READS >=50 BP and save fasta file ####
# Output FASTA path
fasta_out="/data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_PE/input/${sample}_matched_Monly_50bp.fasta"

# Extract matched M sequences ≥50 bp from BAM and save as FASTA
samtools view -h -q 30 "$Sherlock/${sample}/DF_Alignment/${sample}_GCA_024713945.1_DF_mapped_only.bam" | \
awk '
function extract_match(seq, cigar) {
  match_seq = ""
  pos = 1
  while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
    len = substr(cigar, RSTART, RLENGTH - 1) + 0
    op = substr(cigar, RLENGTH, 1)
    if (op == "M") {
      match_seq = match_seq substr(seq, pos, len)
      pos += len
    } else if (op == "I" || op == "S") {
      pos += len
    }
    # D/N/H/P: reference gaps or clipping, no pos change
    cigar = substr(cigar, RLENGTH + 1)
  }
  return match_seq
}
$1 !~ /^@/ {
  mseq = extract_match($10, $6)
  if (length(mseq) >= 50)
    print ">" $1 "\n" mseq
}' > "$fasta_out"


##### Get the cigar READS >=50 BP and save fasta file ####
# Output FASTA path
fasta_out="/data/Sherlock_Lung2/Sunandini_HDM/DF_Chromosome/blast_matched_50-150bp/DF_Singletons/input/${sample}_matched_Monly_50bp.fasta"

# Extract matched M sequences ≥50 bp from BAM and save as FASTA
samtools view -h -q 30 "$Sherlock/${sample}/singletons/${sample}_GCA_024713945.1_DF_singletons_mapped_only.bam" | \
awk '
function extract_match(seq, cigar) {
  match_seq = ""
  pos = 1
  while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
    len = substr(cigar, RSTART, RLENGTH - 1) + 0
    op = substr(cigar, RLENGTH, 1)
    if (op == "M") {
      match_seq = match_seq substr(seq, pos, len)
      pos += len
    } else if (op == "I" || op == "S") {
      pos += len
    }
    # D/N/H/P: reference gaps or clipping, no pos change
    cigar = substr(cigar, RLENGTH + 1)
  }
  return match_seq
}
$1 !~ /^@/ {
  mseq = extract_match($10, $6)
  if (length(mseq) >= 50)
    print ">" $1 "\n" mseq
}' > "$fasta_out"






































































