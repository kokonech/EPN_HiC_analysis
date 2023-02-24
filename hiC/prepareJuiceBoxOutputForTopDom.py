import sys

# example
#inFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_BT/hicpro_result/RELA_BT_chr11_HiC.211117.KRnorm_100Kb.txt"

inFile = sys.argv[1]

vals = inFile.split("/")[-1].split(".")

chrName = vals[1]
binLength = int(vals[-2])

#chrSizeFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/star_index_gencode19/hs37d5_rename.fa.fai"

# 27.04.21 attempt for novel reference
#chrSizeFile="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/EPD210FH_ref/EPD210FH_chr1_chr8_ref_regions.comb_auto.fa.fai"

# 09.08.21 4EP53 RELA
#chrSizeFile="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/c11orf95_RELA_ref/C11orf95_RELA_novel_ref_4EP53.fa.fai"

# 18.08.21 4EP53 RELA
chrSizeFile="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/c11orf95_RELA_ref/C11orf95_RELA_novel_ref_11EP22.fa.fai"


chrLens = {}
for line in open(chrSizeFile):
    items = line.strip().split("\t")
    chrLens[items[0]] = int(items[1])

assert chrName in chrLens
chrLength = chrLens[chrName]

if chrName == "chr11":
    assert chrLength == 135006516

print "Chr name: %s, chr length: %d, bin length: %d" % (chrName, chrLength, binLength)

resFileName = inFile.replace(".txt",  ".TopDom_input.txt")
print "Result file: ", resFileName

resFile = open(resFileName, "w")

locations = set()
valInfo = {}

chrToCheck = resFileName.split("/")[-1].split(".")[1]

for line in open(inFile):
    items = line.strip().split("\t")
    loc1 = int(items[0])
    loc2 = int(items[1])
    #print loc1,loc2
    assert loc1 <= loc2
    locId = "%d-%d" % (loc1,loc2)
    if items[2] == "NaN":
        print "Detcted Nan, skipping..."
        continue
    valInfo[locId] = items[2]
    locations.add(loc1)
    locations.add(loc2)
    
    
finalLoc = max(locations)
numLocs = finalLoc / binLength

print "Loaded locations:"
print "Found locations %d, estimated: %d" % ( numLocs, len(locations) )
print finalLoc

print "Read data. Preparing result..."

start1 = 0
binId1 = 1
bNames = []

writeFirstStr = True

while start1 < chrLength:
    checkLocus1 = start1 + binLength
    start2 = 0
    binId2 = 1
    bVals = []
    while start2 < chrLength:
        checkLocus2 = start2 + binLength
        locOrder = sorted([checkLocus1,checkLocus2])
        locId = "%d-%d" % (locOrder[0], locOrder[1])
        if locId in valInfo:
            value = valInfo[locId]
        else:
            value = "0"
        bName = "HIC_bin%d|hg19|%s:%d-%d" % (binId2, chrName, start2+1,checkLocus2 ) 
        bNames.append(bName)
        #print start1,start2,locId,value
        bVals.append(value)
        start2 += binLength
        binId2 += 1


    #if writeFirstStr:
    #    resFile.write("\t" + "\t".join(bNames) + "\n")
    #    writeFirstStr = False
    #    print "Wrote header.."

    #print bNames[binId1 - 1]
    resFile.write( "%s\t%d\t%d\t%s\n" % (chrName,start1,checkLocus1, "\t".join(bVals) ) )
    binId1 += 1
    start1 += binLength

print "Finished!"









    





