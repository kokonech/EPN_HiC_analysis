# Prepare for HiC BreakFinder calling

#PBS -l nodes=1:ppn=8
#PBS -l mem=50g
#PBS -l walltime=48:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /icgc/dkfzlsdf/analysis/dktk/Ependymoma/external_data/logFolder/epnHicBreakFinder

#Script for aligning Hi-C data using BWA
BWA=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/svAnalysis/jesse_pipeline_hic/bwa_mem_hic_aligner.pl

#Script for manually pairing Hi-C reads
COMBINE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/svAnalysis/jesse_pipeline_hic/two_read_bam_combiner.pl

#Some things that should be hardcoded in, the location of the reference, location of samtools, and the location for picard
#REF=/pbld/netapp/database/human/b38_decoy/hs38d5.fa

module load SAMtools/1.3.1
module load BWA/0.7.15
module load Java/1.8.0_172

# using hg19
REF=/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/indexes/bwa/bwa06_1KGRef/hs37d5.fa
SAMTOOLS=samtools
PICARD=/home/okonechn/B080/okonechn/tools/picard-2.7.1/picard.jar

#Parameters for Java for running picard. If somehow you need more memory, you cna just changed the -Xmx10g to a larger number
TMPDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/external_data/tmpData
JAVA="java -Xmx36g -jar -Djava.io.tmpdir=$TMPDIR"

# Test the procedure

#INPUT=SRR6213067
#Directory where the fastq files are located
#DATADIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/external_data/HiC_public/Caki2/rep2/
#RESDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/Caki2/hicBreakFinderRes2


echo "Processing $INPUT"
echo "Result dir: $RESDIR"
echo "Data dir: $DATADIR"

i=$INPUT

echo "R1: $DATADIR/$i\_1.fastq.gz"
echo "R2: $DATADIR/$i\_2.fastq.gz"

mkdir -p $RESDIR
cd $RESDIR

echo "START ANALYSIS"


#Align read one. You will want to change the "_R1_001.fastq.gz" to match whatever the suffix is for your file
echo
echo "Align read one"
$BWA $DATADIR/$i\_1.fastq.gz $REF 8 |\
$SAMTOOLS view -bS -o $i\_hg19_read1.bam -

echo
echo "Align read two..."
$BWA $DATADIR/$i\_2.fastq.gz $REF 8 |\
$SAMTOOLS view -bS -o $i\_hg19_read2.bam -

echo
echo "Combine samples.."

#Merge files together, sort out reads that have less than MQ=30 (--qual 30), also keeps single end mapping reads (--single, can affect normalization a bit)
$COMBINE --qual 30 --single --file1 $i\_hg19_read1.bam --file2 $i\_hg19_read2.bam |\
$SAMTOOLS view -u -o - - |\
$SAMTOOLS sort -T $i -o $i\_hg19.sorted.bam -

echo
echo "Mark duplicates..."


$JAVA $PICARD MarkDuplicates INPUT=$i\_hg19.sorted.bam OUTPUT=$i\_hg19.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=metrics.$i.txt


echo "DONE!"

