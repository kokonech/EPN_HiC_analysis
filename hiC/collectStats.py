import os.path
import sys
import fnmatch
import glob

sFile = open("/omics/odcf/analysis/OE0290_projects/Ependymoma/scripts/hiC/2018_12_Ependymoma_Tumors_HiC_QC_Dixon.txt")

resultDir = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/tumors"

sIds = []

resFile = open( "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/combinedAnalysis/EPN_HiC.QC_report.230222.txt" , "w" )

resFile.write("ID\tSample\tNumLanes\tGroup\t")
resFile.write("ValidInteractions\tTrans\tCis\tCisShort\tCisLarge\n")


def parseQC(path):
    vals = []
    for line in open(path):
        if line.startswith("valid_interaction\t"):
            continue
        items = line.strip().split("\t")
        vals.append(items[1])        
    return vals



for line in sFile:
    if line.startswith("ID") or len(line.strip()) == 0:
        continue
   
    items = line.strip().split()

    mId = items[0]
    sId = items[7]
    sGroup = items[8].split("-")[0]

    numLanes = int(items[6])
    #if numLanes < 2:
    #    continue
    
    print sId,sGroup

    resLine = "%s\t%s\t%d\t%s\t" % (mId, sId, numLanes, sGroup)
    sPath = "%s/%s/hicpro_result/hic_results/data/%s/%s_allValidPairs.mergestat" % (resultDir, sId, sId, sId )
    assert os.path.exists(sPath)

    vals = parseQC(sPath)
    resLine += "\t".join(vals) + "\n"
    
    resFile.write(resLine)


# collect cell lines

sIds = { "RELPOS_BT" : "RELA", "RELPOS_EP" : "RELA", "PFA_EPD210FH" : "PFA", "FFPE" : "PFA", "SC480" : "PFA" }


for sId in sIds:
    print sId
    resLine = "%s\t%s\t2\t%s\t" % (sId, sId, sIds[sId] )
    if sId == "FFPE" or sId == "SC480":
        path = resultDir.replace("tumors","tumor_relapse") + "/%s_results/hicpro_result/hic_results/data/%s/%s_allValidPairs.mergestat" % (sId, sId, sId)
    else:
        path = resultDir.replace("tumors","") + "%s/hicpro_result/hic_results/data/%s/%s_allValidPairs.mergestat" % (sId, sId, sId)
    vals = parseQC(path)
    resLine += "\t".join(vals) + "\n"
    resFile.write(resLine)

sIds2 = { "CHLA4488-5_primary", "CHLA4488-5_relapse",  "case32_relapse", "case32_primary" }

relapseQcDir = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/tumor_relapse/qc"


for sId in sIds2:
    print sId
    resLine = "%s\t%s\t2\tPFA\t" % (sId, sId )
    
    path = "%s/%s_allValidPairs.mergestat" % (relapseQcDir, sId)
    vals = parseQC(path)
    resLine += "\t".join(vals) + "\n"
    resFile.write(resLine)

  
# EPD210FH
#resLine = "EPD210FH\tEPD210FH\t2\tPFA\t"
#print resLine
#path1 = resultDir.replace("tumors","") + "PFA_EPD210FH/hicpro_result/hic_results/data/PFA_EPD210FH/PFA_EPD210FH_allValidPairs.mergestat"
#vals = parseQC(sPath)
#resLine += "\t".join(vals) + "\n"
#resFile.write(resLine)





           
print "Done!"

           

