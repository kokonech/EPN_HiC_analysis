# NOTE : this code was fully integrated as an option into InTAD R package
# Check optinion "Integration of chromatin loops"
# https://www.bioconductor.org/packages/release/bioc/html/InTAD.html


# PFA EPD210FH

#INFILE="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.txt"


## RELA BT

#INFILE="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_BT/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.txt"


## RELA EP

#INFILE="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_EP/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.txt"


### PFA enhancer


#INFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.PFA_combined_full.qval_0.05.txt
#GRP=PFA


#INFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.PFA_combined_full.qval_0.1.txt
#GRP=PFA

INFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.PFA_combined_full.pval_0.01.txt
GRP=PFA

#INFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.PFA_combined_min2_specific.qval_0.05.txt
#GRP=PFA


## RELPOS enhancers

INFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.RELA_combined_full.qval_0.05.txt
GRP=RELA

INFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.RELA_combined_full.qval_0.1.txt
GRP=RELA

INFILE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.RELA_combined_full.pval_0.01.txt
GRP=RELA


###########
# Analysis

# bin - bin only
INFILE2=${INFILE/.txt/.enh_contacts.txt}

python findGeneEnhancerConnections.py -i $INFILE -g $GRP

python checkConnectionPresence.py $INFILE2 $GRP

Rscript findLoopCorrelations.R $INFILE2

# upstream/downstream bin allowed
INFILE3=${INFILE/.txt/.enh_contacts_closest.txt}

python findGeneEnhancerConnections.py -i $INFILE -g $GRP --spread

python checkConnectionPresence.py $INFILE3 $GRP

Rscript findLoopCorrelations.R $INFILE3





