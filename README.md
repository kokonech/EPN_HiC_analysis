# Ependymoma tumors HiC data analysis

The project contains analysis code for the ependymoma HiC study

Required environment: Linux, Python, R

Tool details (versions, etc) are provided in launch bash commands, for questions please create an issue or contact repo owner directly. 

## HiC data analysis

### Initial processing

Generate HiC-pro analysis configs + input reads data:
*hiC/prepareForHicProTumors.py*

Generate JuiceBox output (JB views are used for mulitple figures e.g. Figure 1f,2a):
*hiC/convertToJuicebox.sh*

Call TADs via TopDom:
*hiC/runTopDom.sh*

### Loops analysis

Call loops via FitHiC2:
*hiC/loopCalling/runFitHiC.sh*

Find gene-enhancer assoications via loops:
*hiC/loopCalling/connectGenesEnhancersViaLoops.sh*

Inspect expression of genes in loops (Figure 1e):
*hiC/loopCalling/checkExprLevel.R*

Differential loops analysis:
*hiC/loopCalling/diffLoopAnalysis.R*

### Structural variants ###

Run alignment for hicBreakFinder:
*hiC/svAnalysis/mapping_hicBreakFinder.sh*

Call SV with hicBreakFinder:
*hiC/svAnalysis/run_hicBreakFinder.sh*

Generate JuiceBox output:
*hiC/svAnalysis/convertHicBreakFinderToJuiceBox.py*

### CTCF - methylation data analysis

Merge CTCF peaks with differntially methylated regions (Figures 4b,c):
*hyperMethCtcf/compareDMRToCTCF.R*

Filter and annotate gene-enhancer in correlation pairs with CTCF loss in loops:
*hyperMethCtcf/ctcfWithinLoops.py*

Filter and annotate gene-superenhancer pairs with CTCF loss in loops:
*hyperMethCtcf/annotateCtcfLossLoopSE.py*
