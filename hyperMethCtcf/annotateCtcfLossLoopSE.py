import os

# filter gene expr

# This is the main script to find possible candidates for experimental validation from CTCF loss for PFA where loop connects super enhancer (SE) to gene

### Input materials

cell_line_expr = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/rna_sequencing/PFA8936/counts/celllines/EPN_cell_line_counts_full.RPKM.gn.txt"

normalConnections = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/externalData/astroscyte_cerebellum/hicpro_result/fithic2_res/loops5KB_initial/fithicoutput_low_bound_10K_res/FitHiC.spline_pass2.res5000.significances.qval_0.05.txt"

contactsResDir = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/"


######
# PFA


# tumor specific

print "Processing PFA..."


inFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/ctcfLossHypomethSE/PFA_gene_SE_loop_connected.with_motifs.txt"


statsCheck = contactsResDir + "FitHiC_res5000.PFA_combined_full.qval_0.05.txt"
#otherCheck = contactsResDir + "FitHiC_res5000.RELA_combined_full.pval_0.01.txt"
otherCheck = contactsResDir + "FitHiC_res5000.RELA_combined_full.qval_0.05.txt"



cellLine = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.pval_0.01.txt"

pubFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/associateGeneEnhancer/PFA_specific_TE_SE_positive_correlated_genes_in_TAD.txt"

seFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/chip_seq_sequencing/enhancer/QSEA_Enhancer/PFA_specific_SE.txt"

deseqFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/rna_sequencing/DESeq/PFA_DEseq2_contrast.with_names.txt"

pfaRelaFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/binInfo/groupSpecificR2/limma_NG_gmkheidelberg06_ep_subgroups_PFA_vs_RELA_R2.txt"
EXP_STATUS_MAIN="05_pf_epn_a >= 09_st_epn_rela"

r2File = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/binInfo/groupSpecificR2/1527860198_NG_t_pfa_vs_other_epn_NG_tpfa_vs_other_epn.txt"
EXP_STATUS="other < 05_pf_epn_a"

nscFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/binInfo/MDT_CRISPR/CRISPR_PFA_fNSC_MDT.txt"

gbmFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/binInfo/MDT_CRISPR/CRISPR_PFA_vs_GBM_MDT.txt"


# load enhancers 

pubEnh = set()

# group specific
for line in open(pubFile):  
    if line.startswith("peakid"):
        continue
    items = line.strip().split("\t")
    sID = "%s_%s" % (items[0], items[3])
#    #print sID
    pubEnh.add(sID)


# load validated
def loadValidCrispr(inFile):
    res = set()
    for line in open(inFile):
        if line.startswith("Gene"):
            continue
        items = line.strip().split("\t") 
        gName = items[0]
        res.add(gName)

    return res

crisprGbm = loadValidCrispr(gbmFile)
crisprNsc = loadValidCrispr(nscFile)
print "Loaded CRISPR PFA vs GBM",len(crisprGbm)
print "Loaded CRISPR PFA vs fNSC",len(crisprNsc)


#print "chr1:24714449:24717135_GRHL3" in pubEnh

# load RNA-seq result

deGenes = set()
for line in open(deseqFile):
    items = line.strip().split("\t")
    gName = items[1]
    deGenes.add(gName)

print "Loaded RNA-seq DEG:",len(deGenes)

mainGenes = set()
for line in open(pfaRelaFile):
    if line.startswith("#"):
        continue
    items = line.strip().split("\t")
    gName = items[1]
    status = items[-1]
    if EXP_STATUS_MAIN in status:
        #print gName
        mainGenes.add(gName)

print "Loaded PFA vs RELA genes:",len(mainGenes) 



r2Genes = set()
for line in open(r2File):
    if line.startswith("#"):
        continue
    items = line.strip().split("\t")
    gName = items[1]
    status = items[-1]
    if EXP_STATUS in status:
        #print gName
        r2Genes.add(gName)

print "Loaded R2 DEG:", len(r2Genes)



# load other connections

def loadLoops(inFile):

    otherRes = set()
    for line in open(inFile):
        if line.startswith("chr1\tfragment"):
            continue
        items= line.strip().split("\t")
        if items[0]  != items[2]:
            continue
        otherRes.add(  "%s:%s-%s" % (items[0], items[1],items[3]) )

    return otherRes

    



clRes = loadLoops(cellLine)
print "Loaded cell line contacts:", len(clRes)

otherRes = loadLoops(otherCheck)
print "Loaded other group contacts:", len(otherRes)


# Load loop stats

qvals = {}
for line in open(statsCheck):
    if line.startswith("chr1\tfragment"):
            continue
    items= line.strip().split("\t")
    if items[0]  != items[2]:
        continue
    qvals[ "%s:%s-%s" % (items[0], items[1],items[3]) ] = "%s\t%s\t%s" % (items[4],items[5],items[6] )

print "Loaded loop stats" 

#normalRes = loadLoops(normalConnections)
#print "Loaded control normal contacts:", len(normalRes)

    

# load expression

clExpr = {}
clExprHeader = "Expr_EPD210FH\tExpr_BT165\tExpr_EP1NS"

for line in open(cell_line_expr):
    if line.startswith("EPD210FH"):
        continue
    items = line.strip().split("\t")
    clExpr[items[0] ] = "\t".join(items[1:])

print "Loaded expression in cell lines..."

# filter gene expression

genesControl = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/rna_sequencing/DESeq/epn_rna_rpkm_ensmbl.filtered.gn.txt"

filteredGenes = set()

for line in open(genesControl):
    if line.startswith("X11EP21"):
        continue
    items = line.strip().split("\t")
    filteredGenes.add(items[0])

print "Loaded gene filtering:", len(filteredGenes)



# check the stats in arima loops

arcLoopDir = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/CaptureHiC/loops_arcplots"


def loadArcLoops(sId):
    
    inFile = "%s/%s.cis_le_2Mb.bedpe" % (arcLoopDir,sId)

    #print inFile
    assert os.path.exists(inFile)

    res = []
    for line in open(inFile):
        items = line.strip().split("\t")
        #print items
        assert (items[0] == items[3] )
        if (items[1] < items[4]):
            res.append( [ items[0], int(items[1]), int(items[2]),
                            int(items[4]), int(items[5]) ]  )
        else:
            res.append( [ items[0], int(items[4]), int(items[5]),
                            int(items[1]), int(items[2]) ]  )
        
        #break
    #print res
    print "Loaded %d arch loops for %s" % (len(res),sId)
    return res


arcBlock = ["5_Aza_1", "5_Aza_2_1", "5_Aza_2_2", "DMSO_1", "DMSO_2" ]
arcInfo = {}
for sId in arcBlock:
    arcInfo[sId] =  loadArcLoops(sId)

def checkArcLoops(loci, arcLoops):
    ls = loci.split(":")
    lchr = ls[0]
    ls2 = ls[1].split("-")
    lstart = int(ls2[0])
    lend = int(ls2[1])
    #print lchr,lstart,lend

    for l in arcLoops:
        if lchr == l[0]:
            p1 = (lstart > l[1]) & (lstart < l[2])
            p2 = (lend > l[3]) & (lend < l[4])
            if (p1 & p2):
                return True

        
    return False





## compute and write the  result


assert ".txt" in inFile
resName = inFile.replace(".txt", ".annotated.txt")
print "Writing result:", resName

resFile = open(resName , "w")

#header  =("Enhancer\tGene\tCor\tP-val\tContactSupp\tEnhType\tPrevEAG\tDESeq2\tR2\tFoundInOther\tCell_line\tFoundInAstroCb\t%s\n" % clExprHeader)

header  = "\tnum_contacts\tpval\tqval\tR2_PFA_vs_RELA\tDESeq2\tR2_global"
header += "\tLoops_RELA\tLoops_cell_line"
header += "\tCRISPR_vs_fNSC\tCRISPR_vs_GBM\t%s\n" % clExprHeader


total = 0
for line in open(inFile):
    items = line.strip().split("\t")
    if line.startswith("Loop"):
        resFile.write( line.strip() + header )
        continue
    total += 1
 
    #print items
    #  control genes
    gName = items[3]
    #if gName not in filteredGenes:
    #    continue

    # extract connection ID
    cId = "%s:%s-%s" % (items[4],items[5], items[6])
    #print cId

    #print items[-4]
    qval = qvals[cId]
    #eag = 1 if sID in pubEnh else 0
    targ = 1 if gName in mainGenes else 0
    deg = 1 if gName in deGenes else 0
    #fio = 1 if sID in otherStatus else 0
    #cl = 1 if sID in clStatus else 0
    r2 = 1 if gName in r2Genes else 0
    #normal = 1 if sID in normalStatus else 0
    cNsc = 1 if gName in crisprNsc else 0
    cGbm = 1 if gName in crisprGbm else 0
    
    otherStatus = cId in otherRes
    clStatus = cId in clRes

    expr = clExpr[gName]
    for sId in arcInfo:
        arcLoop = checkArcLoops(cId,arcInfo[sId]) 
        if arcLoop:
            print "%s: %s" % (sId, cId)

    line = line.strip()
    line += "\t%s\t%d\t%d\t%d" % (qval,targ,deg,r2)
    line += "\t%d\t%d" % (otherStatus,clStatus)
    line += "\t%d\t%d\t%s\n" % (cNsc, cGbm,expr)
    resFile.write(line)

    #sID =  "%s_%s" % (items[-4],items[-3])
    #val = idCount.get(sID, 0)
    #idCount[sID] =  val + 1
    #resData.add( "\t".join(items[-4:]))

       #print "Example", cId
        #otherStatus.add(sID)

    #if cId in clRes:
    #    clStatus.add(sID)

    #if cId in normalRes:
    #    normalStatus.add(sID)


#print "Processed %d contats" % total
#print "Filtered %d results" % len(resData)
#print resData
#print "Found also in other group:", len(otherStatus)
#print "Found in cell line:", len(clStatus)
#print "Found in cerebellum astrocytes:", len(normalStatus)





























