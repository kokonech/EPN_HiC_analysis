# Test submitting jobs

#PBS -l nodes=1:ppn=8
#PBS -l mem=40g
#PBS -l walltime=24:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /omics/odcf/analysis/OE0290_projects/Ependymoma/external_data/logFolder/epnFitHiC

module load Python/3.6.1-foss-2017a
module load R-bundle/20180209-foss-2017a-R-3.4.3

SID=$INPUT
echo "Processing $SID"

###########
# ANALYSIS

# main
#cd /icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/$SID/hicpro_result

# update with Capture
cd /omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/CaptureHiC/per_sample/$SID/hicpro_result

RESDIR=fithic2_res/loops5KB
mkdir -p $RESDIR

PREP_SCRIPT="/home/okonechn/tools/fithic2/fithic/utils/HiCPro2FitHiC.py"

echo 
echo "Step1"
# 5K
cmd1="python $PREP_SCRIPT -i hic_results/matrix/${SID}/raw/5000/${SID}_5000.matrix  -b hic_results/matrix/${SID}/raw/5000/${SID}_5000_abs.bed -s hic_results/matrix/${SID}/iced/5000/${SID}_5000_iced.matrix.biases -o $RESDIR -r 5000"


# adapt for FFPE
#cmd1="python $PREP_SCRIPT -i hic_results/matrix/${SID}/raw/5000/rawdata_5000.matrix  -b hic_results/matrix/${SID}/raw/5000/rawdata_5000_abs.bed -s hic_results/matrix/${SID}/iced/5000/rawdata_5000_iced.matrix.biases -o $RESDIR -r 5000"


echo $cmd1
$cmd1


echo 
echo "Step2"
cd $RESDIR

# 5K 
cmd2="fithic -p 2 -r 5000 -L 10000 -U 5000000 -f fithic.fragmentMappability.gz -i fithic.interactionCounts.gz -t fithic.biases.gz  -o fithicoutput_low_bound_10K_res"


echo $cmd2
$cmd2


echo "Done!"




