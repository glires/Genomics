#!/bin/sh -x
#
# ngsplg - NGS pipeline for genomic sequences
#            (e.g. exome and target re-sequencing)
#
# usage:
#   $ cd Documents/
#   $ mkdir ProjectName
#   $ cd ProjectName
#   $ cp -p /mnt/data/Unaligned/Project/*_R[12]_*.fastq.gz ./
#   $ nohup ngsplg 1> log01.txt 2> log02.txt &
#
# Note that ./seq/, ./qc/, ./aln/, ./log/, and ./var/ will be overwritten.

echo $SHELL

kit=SureSelect
	# kit=Nextera
	# kit should be SureSelect, TruSeq, or Nextera
shortest=36
hit=unique
	# hit=multiple
	# hit should be either unique or multiple
biol=/usr/local/biol
reference=$biol/References/hg19/hs37d5.fa
baits=$biol/References/baits/SureSelect_v4_80Mb.bed
	# baits=$biol/References/baits/SureSelect_v4_51Mb.bed
	# baits=$biol/References/baits/baits20120903.bed
padding=0
	# Note that Agilent SureSelect BED files have been already
	# expanded 100 bp from the borders, so padding=100 means
	# expanding 200 bp.
project=`basename $PWD`

export	LANG=C
date 1>&2	# from stdout to stderr
rm -fr log/
mkdir log

ls -l *_R1_*.fastq.gz 1>&2
ls -l *_R2_*.fastq.gz 1>&2

gzip -d *_R1_*.fastq.gz
gzip -d *_R2_*.fastq.gz

cat *_R1_*.fastq > ${project}_R1_0.fastq
cat *_R2_*.fastq > ${project}_R2_0.fastq

date 1>&2
rm -fr seq/
mkdir seq/

mv ${project}_R1_0.fastq ${project}_R2_0.fastq seq/

# for fastq in *_R[12]_*.fastq
# do
#   gzip $fastq
# done
# mv -i *_R[12]_*.fastq.gz seq/

date 1>&2
cd seq/

# rm -fr ../qc/
# mkdir ../qc

# $biol/FastQC/fastqc --outdir ../qc ¥
#   ${project}_R1_0.fastq ${project}_R2_0.fastq 1>&2

## Check FASTQ
date 1>&2
$biol/NGS/jf0010.pl ¥
  -a ${project}_R1_0.fastq -b ${project}_R2_0.fastq ¥
  -c ${project}_R1_1.fastq -d ${project}_R2_1.fastq -e ${project}_R3_1.fastq
date 1>&2
gzip ${project}_R1_0.fastq
gzip ${project}_R2_0.fastq

#### Trimmomatic
date 1>&2
$biol/NGS/make_adaptor $kit
  if [ $? -ne 0 ]; then	# check the exit status of the previous command
    rm -f adaptor.fa
    cd ../
    exit
  fi
java -classpath $biol/Trimmomatic-0.22/trimmomatic-0.22.jar ¥
  org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 ¥
  -trimlog ../log/log03.txt ${project}_R1_1.fastq ${project}_R2_1.fastq ¥
  ${project}_R1_2.fastq ${project}_R3_2.fastq ¥
  ${project}_R2_2.fastq ${project}_R4_2.fastq ¥
  ILLUMINACLIP:adaptor.fa:1:40:10 LEADING:8 TRAILING:8 SLIDINGWINDOW:8:16 ¥
  MINLEN:$shortest 1>&2
rm -f adaptor.fa

#### compress sequence files
date 1>&2
# gzip ${project}_R1_1.fastq ${project}_R2_1.fastq ${project}_R3_1.fastq
# gzip ${project}_R3_2.fastq ${project}_R4_2.fastq
mv ${project}_R1_2.fastq ${project}_R1.fastq
mv ${project}_R2_2.fastq ${project}_R2.fastq
# gzip ${project}_R?_?.fastq

#### BWA aln
date 1>&2
rm -fr ../aln/
mkdir ../aln
$biol/bwa-0.6.2/bwa aln -t 2 $reference ¥
  ${project}_R1.fastq -f ../aln/${project}_R1.sai ¥
  1> ../log/log04.txt 2> ../log/log05.txt
$biol/bwa-0.6.2/bwa aln -t 2 $reference ¥
  ${project}_R2.fastq -f ../aln/${project}_R2.sai ¥
  1> ../log/log06.txt 2> ../log/log07.txt

#### BWA sampe
date 1>&2
$biol/bwa-0.6.2/bwa sampe -n 4 -N 8 -a 400 -o 80000 ¥
  -r "@RG¥tID:MFB¥tSM:${project}¥tPL:Illumina" ¥
  -P -f ../aln/${project}.sam ¥
  $reference ¥
  ../aln/${project}_R1.sai ../aln/${project}_R2.sai ¥
  ${project}_R1.fastq ${project}_R2.fastq ¥
  1> ../log/log08.txt 2> ../log/log09.txt

#### select unique hits and SAMtools
date 1>&2
cd ../aln/
if [ $hit = multiple ]; then
  cp -p ${project}.sam ${project}_1.sam
else
  $biol/NGS/jf0005.pl ${project}.sam > ${project}_1.sam
fi
$biol/samtools-0.1.18/samtools view -bS ${project}_1.sam ¥
  1> ${project}_0.bam 2> ../log/log10.txt
$biol/samtools-0.1.18/samtools sort ${project}_0.bam ¥
  ${project}_1 1> ../log/log11.txt 2> ../log/log12.txt

#### Picard MarkDuplicates
date 1>&2
java -Xmx6G -jar $biol/picard-tools-1.83/MarkDuplicates.jar ¥
  ASSUME_SORTED=true REMOVE_DUPLICATES=true ¥
  INPUT=${project}_1.bam OUTPUT=${project}_2.bam ¥
  METRICS_FILE=../log/log13.txt VALIDATION_STRINGENCY=LENIENT ¥
  1> ../log/log14.txt 2> ../log/log15.txt
$biol/samtools-0.1.18/samtools index ${project}_2.bam

#### GATK RealignerTargetCreator
date 1>&2
java -Xmx6G -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type RealignerTargetCreator -R $reference ¥
  --input_file ${project}_2.bam --known $biol/References/hg19/indels.vcf ¥
  --intervals $baits --interval_padding $padding -o ${project}.intervals
	# bait file should not be in hg19 but in b37
ls -l ${project}.intervals

#### GATK IndelRealigner
date 1>&2
rm -fr tmp/
mkdir tmp
java -Xmx6G -Djava.io.tmpdir=./tmp ¥
  -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type IndelRealigner -R $reference ¥
  --input_file ${project}_2.bam -log ../log/log16.txt ¥
  -targetIntervals ${project}.intervals -o ${project}_3.bam ¥
  1> ../log/log17.txt 2> ../log/log18.txt

#### Picard FixMateInformation; SAMtools
date 1>&2
java -Xmx6G -jar $biol/picard-tools-1.83/FixMateInformation.jar ¥
  INPUT=${project}_3.bam SO=coordinate VALIDATION_STRINGENCY=SILENT ¥
  TMP_DIR=./tmp 1> ../log/log19.txt 2> ../log/log20.txt
	# ${project}_3.bam is changed
$biol/samtools-0.1.18/samtools index ${project}_3.bam
rm ${project}_3.bai
rm -fr tmp

#### GATK BaseRecalibrator and PrintReads; SAMtools
date 1>&2
java -Xmx6G -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type BaseRecalibrator --input_file ${project}_3.bam ¥
  -R $reference ¥
  --knownSites $biol/References/hg19/dbSNP_201206.vcf -log ../log/log21.txt ¥
  --covariate ContextCovariate --out recal_data.grp
java -Xmx6G -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type PrintReads -R $reference ¥
  --input_file ${project}_3.bam ¥
  -BQSR recal_data.grp --out ${project}.bam
	# ${project}.bai is simultaniously created.

#### GATK DepthOfCoverage
date 1>&2
java -Xmx6G -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type DepthOfCoverage -R $reference --input_file ${project}.bam ¥
  -o ${project} --intervals $baits --interval_padding $padding ¥
  -omitBaseOutput -ct 10 -ct 20 -ct 30 -ct 40
cd ../

#### GATK UnifiedGenotyper and VariantFiltration for SNV
date 1>&2
rm -fr var/
mkdir var
cd var/
rm -fr tmp/
mkdir tmp
java -Djava.io.tmpdir=./tmp -Xmx6G ¥
  -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type UnifiedGenotyper -R $reference ¥
  --input_file ../aln/${project}.bam --num_threads 4 ¥
  --intervals $baits --interval_padding $padding ¥
  --genotype_likelihoods_model SNP --out ${project}_0.vcf ¥
  --annotation AlleleBalance --annotation DepthOfCoverage ¥
  --standard_min_confidence_threshold_for_calling 40 ¥
  --standard_min_confidence_threshold_for_emitting 40 ¥
  --downsample_to_coverage 400
java -Djava.io.tmpdir=./tmp -Xmx6G ¥
  -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type VariantFiltration -R $reference ¥
  -V ${project}_0.vcf --out ${project}_1.vcf --clusterWindowSize 10 ¥
  --filterExpression "QUAL < 30.0 || QD < 5.0" ¥
                                             --filterName "QLflt" ¥
  --filterExpression "QD < 2.0"              --filterName "QDflt" ¥
  --filterExpression "MQ < 40.0"             --filterName "MQflt" ¥
  --filterExpression "FS > 60.0"             --filterName "FSflt" ¥
  --filterExpression "HaplotypeScore > 13.0" --filterName "HTflt"

#### GATK UnifiedGenotyper and VariantFiltration for indel
date 1>&2
java -Djava.io.tmpdir=./tmp -Xmx6G ¥
  -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type UnifiedGenotyper -R $reference ¥
  --input_file ../aln/${project}.bam --num_threads 4 ¥
  --intervals $baits --interval_padding $padding ¥
  --genotype_likelihoods_model INDEL --out ${project}_2.vcf ¥
  --annotation AlleleBalance --annotation DepthOfCoverage ¥
  --standard_min_confidence_threshold_for_calling 40 ¥
  --standard_min_confidence_threshold_for_emitting 40 ¥
  --downsample_to_coverage 400
java -Djava.io.tmpdir=./tmp -Xmx6G ¥
  -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type VariantFiltration -R $reference ¥
  -V ${project}_2.vcf --out ${project}_3.vcf ¥
  --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" ¥
                                              --filterName "MDflt" ¥
  --filterExpression "QUAL < 10"              --filterName "QLflt" ¥
  --filterExpression "QD < 2.0"               --filterName "QDflt" ¥
  --filterExpression "FS > 200.0"             --filterName "FSflt"

#### GATK CombineVariants
date 1>&2
java -Djava.io.tmpdir=./tmp -Xmx6G ¥
  -jar $biol/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar ¥
  --analysis_type CombineVariants -R $reference ¥
  --variant:SNV ${project}_1.vcf --variant:INDEL ${project}_3.vcf ¥
  --filteredAreUncalled --genotypemergeoption PRIORITIZE ¥
  --rod_priority_list SNV,INDEL --out ${project}_4.vcf
rm -fr tmp

#### snpEff
date 1>&2
cp -p $biol/snpEff_3_1/snpEff.config ./
java -Xmx6G -jar $biol/snpEff_3_1/snpEff.jar eff -v -i vcf -o vcf ¥
  GRCh37.69 ${project}_4.vcf > ${project}.vcf
rm snpEff.config

#### IGVTools
$biol/IGVTools/igvtools index ${project}.vcf 1> ../log/log23.txt ¥
  2> ../log/log24.txt
mv -f igv.log ../log/log25.txt
cd ../

$biol/NGS/mv_logs &
date 1>&2

# log/log01.txt	ngsplg stdout
# log/log02.txt	ngsplg stderr
# log/log03.txt	Trimmomatic log
# log/log04.txt	bwa aln R1 stdout
# log/log05.txt	bwa aln R1 stderr
# log/log06.txt	bwa aln R2 stdout
# log/log07.txt	bwa aln R2 stderr
# log/log08.txt	bwa sampe stdout
# log/log09.txt	bwa sampe stderr
# log/log10.txt	samtools view stderr
# log/log11.txt	samtools sort stdout
# log/log12.txt	samtools sort stderr
# log/log13.txt	picard metrics file
# log/log14.txt	picard stdout
# log/log15.txt	picard stderr
# log/log16.txt	GATK IndelRealigner log
# log/log17.txt	GATK IndelRealigner stdout
# log/log18.txt	GATK IndelRealigner stderr
# log/log19.txt	picard FixMateInformation stdout
# log/log20.txt	picard FixMateInformation stderr
# log/log21.txt	GATK BaseRecalibrator log
# log/log23.txt	IGVTools stdout
# log/log24.txt	IGVTools stderr
# log/log25.txt	IGVTools log
