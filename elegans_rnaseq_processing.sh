#!/bin/bash

#SBATCH --partition=128x24
#SBATCH --mail-type=NONE
#SBATCH --mail-user=ccockrum@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G
#SBATCH --output="./logs/%A-%a.out"

#### load modules ####

module unload python
module load samtools 
module load miniconda3

#### $1 is a project directory to process from sbatch submission. Raw reads are stored in subdirectory 'reads' ####

dir=$1

if [[ ! -d "${dir}/map" ]]; then
   mkdir "${dir}/map"
fi

if [[ ! -d "${dir}/counts" ]]; then
   mkdir "${dir}/counts"
fi

if [[ ! -d "${dir}/bigWigs" ]]; then
   mkdir "${dir}/bigWigs"
fi

data=${dir}/reads
genomeIndex=/hb/home/ccockrum/omes/ws220/genome/ws220_ercc
transcriptome=/hb/home/ccockrum/omes/ws220/transcriptome/ws220_ercc.gtf
map=${dir}/map
counts=${dir}/counts
bigWigs=${dir}/bigWigs


#### activate conda environment ####

. /hb/software/apps/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# Define array

        cd $data
        fileArray=(`ls *R1*.fastq.gz`)

# Define read files

        read1=${fileArray[$SLURM_ARRAY_TASK_ID]}
        echo 'Read 1 file is' $read1

        read2=${read1/R1/R2}
        echo 'Read 2 file is' $read2

# Define name of sample to be used for entire script

        name=${read1%%_S*_*}
        echo 'Name of sample is' $name

  if [[ ! -d "${data}/fastqc" ]]; then
     mkdir "${data}/fastqc"
  fi

#### qc reads with fastqc

  	fastqc $read1 -o ./fastqc
 	fastqc $read2 -o ./fastqc

# Map reads using hisat2

        hisat2 --new-summary -p ${SLURM_CPUS_PER_TASK} --no-unal --no-mixed --no-discordant  -x $genomeIndex -1 $read1 -2 $read2 -S ${map}/${name}.sam &> ${map}/${name}_hisat.log

# Sort and fix mates

        cd $map
        mappedFile=${name}.sam
        echo 'Mapped file is' $mappedFile
        samtools view -q 10 -b -h $mappedFile | samtools sort -n -@ ${SLURM_CPUS_PER_TASK} -O BAM -o ${name}_sorted.bam


# Coordinate sort for duplicate removal

        samtools fixmate -r -m ${name}_sorted.bam ${name}_fixMate.bam
        samtools sort -@ 12 -T ${name}_tmp -O BAM -o ${name}_fixMate_positionSorted.bam ${name}_fixMate.bam


# Remove duplicates

        samtools markdup -@ ${SLURM_CPUS_PER_TASK} -T ${name}_tmp -r -s ${name}_fixMate_positionSorted.bam ${name}_dedup.bam

# Remove intermediate files

        rm ${name}_fixMate*
        rm ${name}*.sam

# Map QC

          if [[ ! -d "${map}/qualimap" ]]; then
            mkdir "${map}/qualimap"
          fi

         qualimap bamqc -bam ${name}_dedup.bam -gff $transcriptome -outdir ./qualimap/${name}_bamqc
         qualimap rnaseq -bam ${name}_dedup.bam -pe -s -gtf $transcriptome -outdir ./qualimap/${name}_rnaseq


# Count reads over transcriptome annotation

        featureCounts -T ${SLURM_CPUS_PER_TASK} -B -C -t exon -g gene_id --extraAttributes gene_biotype -p -a $transcriptome -o ${counts}/${name}_count.txt ${name}_dedup.bam

        if [[ ! -d "${counts}/summaries" ]]; then
           mkdir "${counts}/summaries"
        fi

        mv ${counts}/*summary ${counts}/summaries

# Make coverage tracks

	cd $map	

	samtools index ${name}_dedup.bam

	SCALEFACTOR=$(awk -v a="$name" 'BEGIN {OFS='\t'} ; $1 == a {print $3}' ${dir}/sizeFacs.txt)
      	INVERSE_SCALEFACTOR=1/$SCALEFACTOR
    	bamCoverage -b ${name}_dedup.bam -p 12 --exactScaling --centerReads --scaleFactor $INVERSE_SCALEFACTOR -o ${bigWigs}/${name}.bw --binSize 5 --smoothLength 15

      



