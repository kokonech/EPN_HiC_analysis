# K.O.: This script generates result for loop conntected gene and ennancers pairs with hyper methylation CTCF loss within. Note that PFA and RELA should be controlled manually.
# Result for RELA was zero.

targetPairs = set()

contactsResDir = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/"

# manual control (also log2FC params below!!!)

# PFA
reference = contactsResDir + "FitHiC_res5000.PFA_combined_full.pval_0.01.enh_contacts_closest.cor_analysis.annotated.txt"
inFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/wgbs/methAnalysis/FitHiC_res5000.PFA_combined_full.pval_0.01.loops_with_CTCF_DMR_in_between.041219.txt"
targGrp = "PFA"


# RELA 
#reference = contactsResDir + "FitHiC_res5000.RELA_combined_full.pval_0.01.enh_contacts_closest.cor_analysis.annotated.txt"
#inFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/wgbs/methAnalysis/FitHiC_res5000.RELA_combined_full.pval_0.01.loops_with_CTCF_DMR_in_between.050320.txt"
#targGrp = "RELA"


for line in open(reference):
   if line.startswith("Enhancer"):
        continue
   items = line.strip().split("\t")
   enh = items[0]
   gene = items[1]
   targetPairs.add("%s|%s" % (enh,gene))
#print targetPairs


# candidates control based on FitHiC EAG main selection

ctrlFile = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/combinedAnalysis/2019_08_27_list_full.txt"

ctrlGenes = set()
for line in open(ctrlFile):
    items = line.strip().split("\t")
    if items[2] != "FitHiC EAG":
        continue
    if items[1] == targGrp:
        ctrlGenes.add(items[0])

print "Control genes:", len(ctrlGenes)

found = set()


resFile = open(inFile.replace(".txt", ".EAG_selected.v3.txt") , "w")
resFile2 = open(inFile.replace(".txt", ".bedpe"), "w")



cId = 1
for line in open(inFile):
    if line.startswith("chr1\t"):
        resFile.write(line)
        continue
    items = line.strip().split("\t")       
    enhBlocks = items[8].split(":")
    enhLoci = enhBlocks[1].split("-")
    enh = "%s:%d-%s" % (enhBlocks[0], int(enhLoci[0]) - 1, enhLoci[1])
    gName = items[9]
    pairId = "%s|%s" % (enh,gName)
    qval = float(items[7])
    ctcfDiff = float(items[20])
    methDiff = float(items[-1])

    #control CTCF loss hypermetylation
    valid = (ctcfDiff < 0) and (methDiff > 0) # PFA
    #valid = (ctcfDiff > 0) and (methDiff < 0) #RELA    

    if (pairId in targetPairs) and valid:
        resFile2.write("%s\t%s\t%s\t" % (items[1], items[2],items[2]) )
        resFile2.write("%s\t%s\t%s\t" % (items[1], items[4],items[4]) )
        resFile2.write("c%d\t%s\t.\t.\t%s\n" % (cId, items[-2], items[-1] ) )
        cId += 1
        if gName in ctrlGenes:
            resFile.write(line)
            #print pairId
            #print "Found  %s:%s-%s" % ( items[1], items[2], items[4] )
        found.add(pairId)

    
print ("Done!")



