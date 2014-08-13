#!/bin/sh
# pipeline for snp calling from WES illumina data
# Matt Bawn 13-08-14 mattbawn@gmail.com
# 1 run takes around 600 mins

name=NIH-Exome

fastq=$name-filtered.fq
align=$name-aligned.sai
sam=$name.sam
bam=$name.bam
sortedBam=$name-sorted.bam
reorder=$name-reorder.bam
dedup=$name-dedup.bam
finalbam=$name-final.bam
readgroup=$name-readgroup.bam
rawbcf=$name-raw.bcf
snpvcf=$name-snp.vcf

# Check reference genome present and if not download

if  [ ! -f wg.fa ]
    then wget http://hgdownload.cse.ucsc.edu//goldenPath/hg19/bigZips/chromFa.tar.gz

    tar zvfx chromFa.tar.gz
    cat *.fa > wg.fa
    rm chr*.fa
fi

# FastQC Analysis of raw input data

if [ ! -d fastqc_test ]
    then mkdir fastqc_test
    cd fastqc_test
    cp ../$fastq .
    fastqc $fastq
    cd ..
fi

# Create BWA index file

if [ ! -f hg19bwaidx.sa ]
    then bwa index -p hg19bwaidx -a bwtsw wg.fa 
fi

# Align using BWA

if [ ! -f $align ]
    then bwa aln hg19bwaidx $fastq > $align
fi

# Convert to SAM file

if [ ! -f $sam ]
    then bwa samse hg19bwaidx $align $fastq > $sam
fi

# SAM To BAM conversion and SORTING in SAMTOOLS

if [ ! -f $bam ]
    then samtools view -bS  -o $bam $sam
fi

# Sort BAM file SAMTOOLS

if [ ! -f $sortedBam ]
    then samtools sort -@ 5 -m 4G  $bam -f $sortedBam
fi

# Reorder SAM using picard-tools

if [ ! -f $reorder ]
    then java -Xmx1g -jar /home/mattbawn/picard-tools-1.118/SortSam.jar \
                            INPUT=$sortedBam \
                            OUTPUT=$reorder \
                            SORT_ORDER=coordinate 
fi

# Find and Remove Duplicates picard-tools

if [ ! -f $dedup ]
    then java -Xmx1g -jar /home/mattbawn/picard-tools-1.118/MarkDuplicates.jar \
                            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                            METRICS_FILE=out.metrics \
                            REMOVE_DUPLICATES=true \
                            ASSUME_SORTED=true  \
                            VALIDATION_STRINGENCY=LENIENT \
                            INPUT=$reorder \
                            OUTPUT=$dedup
fi

# Add readgroups using picard-tools

if [ ! -f $readgroup ]
    then java -Xmx1g -jar /home/mattbawn/picard-tools-1.118/AddOrReplaceReadGroups.jar \
			    INPUT=$dedup \
			    OUTPUT=$readgroup \
			    RGID=group1 RGLB= lib1 RGPL=illumina RGPU=unit1 RGSM=sample1
fi

# Reorder Sam to UCSC reference using picard-tools

if [ ! -f $finalbam ]
    then java -Xmx1g -jar /home/mattbawn/picard-tools-1.118/ReorderSam.jar \
                            INPUT=$readgroup \
                            OUTPUT=$finalbam \
                            REFERENCE=/home/mattbawn/Documents/ucsc.hg19.fasta
fi


# Index bam file

if [ ! -f *.bai ]
    then java -Xmx1g -jar /home/mattbawn/picard-tools-1.118/BuildBamIndex.jar \
                           INPUT=$finalbam 
fi

# Create GATK dictionary

if [ ! wg.fa.dict ]
    then java -Xmx1g -jar /home/mattbawn/picard-tools-1.118/CreateSequenceDictionary.jar \
                            R=wg.fa \
                            O=wg.fa.dict
fi

# Index Reference again

if [ ! -f *.fai ]
    then samtools faidx wg.fa
fi


# Call SNPs using bcftools

if [ ! -f $snpvcf ]
    then java -Xmx1g -jar /home/mattbawn/GATK/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \
                   -R /home/mattbawn/Documents/ucsc.hg19.fasta \
                   -T HaplotypeCaller \
                   -I $finalbam -minPruning 3 \
                   -o $snpvcf 
fi


