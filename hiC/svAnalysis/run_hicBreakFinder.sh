# Main HiC BreakFinder calling

#PBS -l nodes=1:ppn=8
#PBS -l mem=50g
#PBS -l walltime=24:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /icgc/dkfzlsdf/analysis/dktk/Ependymoma/external_data/logFolder/epnHicBreakFinder

TOOL=/home/okonechn/B080/okonechn/tools/hic_breakfinder-master/src/hic_breakfinder

INTEREXP=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/hicBreakFinderTest/inter_expect_1Mb.hg19.no_chr.txt

INTRAEXP=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/hicBreakFinderTest/intra_expect_100kb.hg19.no_chr.txt

module load BamTools/2.4.1-foss-2017a

echo "Processing $INPUT"
echo "Result dir: $RESDIR"

INFILE=$RESDIR/${INPUT}_hg19.nodup.bam

cd $RESDIR

echo "START ANALYSIS"

cmd="$TOOL --bam-file $INFILE --exp-file-inter $INTEREXP --exp-file-intra $INTRAEXP --name $INPUT"

echo $cmd
$cmd
echo

echo "DONE!"


