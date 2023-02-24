# this script calls TopDom on JuiceBox HiC data

#module load R-bundle/20180209-foss-2017a-R-3.4.3

# stable parameters

JUICEBOX=/home/okonechn/tools/juicebox/juicer_tools_0.7.5.jar
SCRIPT1=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/prepareJuiceBoxOutputForTopDom.py
SCRIPT2=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/callTopDom.R
SCRIPT3=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/convertTopDomToJuiceBox.py

# specific options

# RELA_BT
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_BT/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_BT/hicpro_result/RELPOS_BT_allValidPairs.hic

# RELA_EP
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_EP/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_EP/hicpro_result/RELPOS_EP_allValidPairs.hic

# IMR90
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/IMR90_RenLab_rep1/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/IMR90_RenLab_rep1/hicpro_result/IMR90_rep1_allValidPairs.hic

# EPD210FH
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result/PFA_EPD210FH_allValidPairs.hic

#Cerebellum Astro rep1

#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/astroscyte_cerebellum/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/astroscyte_cerebellum/hicpro_result/AstroCb_rep1_allValidPairs.hic

# Cerebellum Astro rep2
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/astroscyte_cerebellum_rep2/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/astroscyte_cerebellum_rep2/hicpro_result/AstroCb_rep2_allValidPairs.hic

# SK-N-SH
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/SK-N-SH/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/SK-N-SH/hicpro_result/SK_N_SH_allValidPairs.hic

# EPD210FH rebuild ref
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result/PFA_EPD210FH_chr1_chr8_ConvPairsAlleleProp.hic
# for manual
#chr=8_EPD210FH_REGION_F

# chr11 4EP53 ref rebuild
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/4EP53/hicpro_result/TopDom_results
#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/4EP53/hicpro_result/4EP53_chr11_ConvPairsAlleleProp.hic

# chr11 11EP22 ref rebuild
TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/11EP22/hicpro_result/TopDom_results
HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/11EP22/hicpro_result/11EP22_chr11_ConvPairsAlleleProp.hic


############
# use input

#HICFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/${SID}/hicpro_result/${SID}_allValidPairs.hic
#TOPDOMDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/${SID}/hicpro_result/TopDom_results


echo "Input file: $HICFILE"
echo "Result dir: $TOPDOMDIR"

mkdir -p $TOPDOMDIR

# size of the bin
#BINSIZE=100000
#BINSIZE=50000
BINSIZE=25000
#BINSIZE=10000
#BINSIZE=5000

# normalization
NORM_TYPE="KR"

# 2 chromsomes: development test
#declare -a CHR_IDS=("chr11" "chr20" )

# all chromosomes (wtihout chrX and chrY)
#declare -a CHR_IDS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" )

# typical analysis
#HICBASE=`basename $HICFILE`
#echo $HICBASE
#SID="${HICBASE/_allValidPairs.hic/}"

# manual for paricular case
#SID=PFA_EPD210FH_chr1_chr8
#declare -a CHR_IDS=("chr8_EPD210FH_region_f")

#SID=4EP53_c11orf95
#declare -a CHR_IDS=("chr11_4EP53_c11orf95_f")

SID=11EP22_c11orf95
declare -a CHR_IDS=("chr11_11EP22_c11orf95_f")


echo "SID: $SID"

RESDIR=$TOPDOMDIR/norm${NORM_TYPE}_${BINSIZE}

if [ ! -d "$RESDIR" ]; then
    echo "Creating directory..."
    mkdir $RESDIR
fi

##  loop through the chromosome array

for chr in "${CHR_IDS[@]}"
do
    
    echo
    echo "Processing $chr"
    echo

    echo "Extract connections from JuiceBox result..."    
    RES1=$RESDIR/${SID}.${chr}.norm_${NORM_TYPE}.$BINSIZE.txt
    if [ ! -f "$RES1" ];
    then
        cmd1="java -jar $JUICEBOX dump observed $NORM_TYPE $HICFILE $chr $chr BP $BINSIZE $RES1"
        echo $cmd1
        $cmd1
    else
        echo "Result exists, skipping"
    fi
 
    echo "Convert to TopDom format..."
    RES2=${RES1/.txt/.TopDom_input.txt}
    if [ ! -f "$RES2" ];
    then
        cmd2="python $SCRIPT1 $RES1"
        echo $RES2
        echo $cmd2
        $cmd2
    else
        echo "Result exists, skipping"
    fi

    echo "Run TopDom..."
    RES3=${RES2/TopDom_input.txt/TopDom_result.domain}
    if [ ! -f "$RES3" ];
    then
        # old version does not work
        #cmd3="/home/okonechn/tools/R-3.3.1/bin/Rscript $SCRIPT2 $RES2"
        
        # R-bundle/20180209-foss-2017a-R-3.4.3
        cmd3="Rscript $SCRIPT2 $RES2"
        echo $RES3
        echo $cmd3
        $cmd3
    else
        echo "Result exists, skipping"
        echo $RES3
    fi

    echo "Convert to JuiceBox 2D annotation..."
    RES4=${RES3/.domain/.juicebox.txt}
    if [ ! -f "$RES4" ];
    then
        cmd4="python $SCRIPT3 $RES3"
        echo $RES4
        echo $cmd4
        $cmd4
    else
        echo "Result exists, skipping"
    fi
    
    echo

done



echo "Combine final result..."

RES5=$RESDIR/${SID}.norm_${NORM_TYPE}.${BINSIZE}.TopDom_result.juicebox.txt

if [ ! -f "$RES5" ];
then
 
    cd $RESDIR
    
    RESNAME=${SID}.norm_${NORM_TYPE}.${BINSIZE}
    cat ${SID}.*.norm_${NORM_TYPE}.${BINSIZE}.TopDom_result.bed | grep domain > ${RESNAME}.TopDom_TADs.bed
    
    head -n 1 ${SID}.chr1.norm_${NORM_TYPE}.${BINSIZE}.TopDom_result.juicebox.txt > ${RESNAME}.TopDom_result.juicebox.txt
    
    cat ${SID}.*.norm_${NORM_TYPE}.${BINSIZE}.TopDom_result.juicebox.txt | grep -v "comment" >> ${RESNAME}.TopDom_result.juicebox.txt

else
    echo "Result exists, skipping..."
fi

  

