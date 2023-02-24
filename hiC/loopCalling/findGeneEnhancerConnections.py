# gen code TSS

import sys
import argparse

#aconFile = sys.argv[1]
descriptionText = "The script performs association of genes and possible targets from input of contacts."

parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", action="store", required="true", dest="inFile",
        help="Input contacts file.")

parser.add_argument("-g", action="store", required=False, dest="selGroup",
        help="Group for specific enhancers: PFA or RELA")

parser.add_argument("--spread", action="store_true", default=False, dest="spreadMode", help="Option to activate spread check: upstream and downsream of the bin")

args = parser.parse_args()


print "Input contacts file:", args.inFile
conFile = args.inFile
print "Group selected:", args.selGroup
print "Upstream/donswtream spread mode active:", args.spreadMode

geneBinsFile="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/binInfo/genecode19_tss_per_5kBP.txt"

geneBins = {}

for line in open(geneBinsFile):
    items = line.strip().split("\t")
    binId = items[0]
    gName = items[1]

    if binId in geneBins:
        geneBins[binId].append(gName)
    else:
        geneBins[binId] = [gName]

print "Loaded gene bins %d  bins" %  len(geneBins)

# 07.07.18 PFA

status = ""
if args.selGroup == "PFA":
    enhBinsFile="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/binInfo/PFA_specific_peaks_per_5kBP.txt"
elif args.selGroup == "RELA":
    enhBinsFile="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/binInfo/RELA_specific_peaks_per_5kBP.txt"
else:
    print "Using common enhancers"
    status = "common_"
    enhBinsFile="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/binInfo/PFA_RELA_group_specific_enhancer_peaks_per_5kBP.txt"
    

print "Enhancer bins data", enhBinsFile


enhBins = {}
for line in open(enhBinsFile):
    items = line.strip().split("\t")
    binId = items[0]
    enhId = items[1]
    enhBins[binId] = enhId

print "Loaded enhancer %d  bins" %  len(enhBins)

# confident connections

# PFA

#resFile = open( conFile.replace(".txt", ".EPD210FH_enh_contacts_closest.txt"), "w")


# RELA EP

#conFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_EP/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.txt"

#resFile = open( conFile.replace(".txt", ".RELA_EP_enh_contacts_closest.txt"), "w")

# RELA BT

#conFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_BT/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.txt"

   
if args.spreadMode:
    # bin and upstream/downstream closest
    binWindows = [(-1,-1),(-1,0), (-1,1),(0,-1), (0,0),(0,1), (1,-1), (1,0),(1,1)]
else:
    # only the same bin
    binWindows = [(0,0)]
    
 
if len(binWindows) == 1:
    print "Checking only bin-bin connections..."
    resFile = open( conFile.replace(".txt", ".%senh_contacts.txt" % status ), "w")
else:
    print "Checking with upstream/downstream bin..."
    resFile = open( conFile.replace(".txt", ".%senh_contacts_closest.txt" % status), "w")

# NOTE: what about 2 genes per bin?

novelConnections = []

k = 0
LIM = 2500
found = 0

for line in open(conFile):
    line = line.strip()
    if k == 0:
        resFile.write(line + "\ttarget1\ttarget2\n")
        k += 1
        continue
    items = line.split("\t")
    chr1 = items[0]
    mPos1 = int(items[1])

    chr2 = items[2]
    mPos2 = int(items[3])
 
   
    contacts = []
    #print "Start ", mPos1,mPos2

    foundPairs= set()
    for binWin in binWindows:
        mPos1 = int(items[1]) + 2*binWin[0]*LIM
        mPos2 = int(items[3]) + 2*binWin[1]*LIM

        pos1 = "%s:%d-%d" % (chr1, mPos1 - LIM + 1, mPos1 + LIM)
        pos2 = "%s:%d-%d" % (chr2, mPos2 - LIM + 1, mPos2 + LIM)

        #print pos1, pos2
        if pos1 in geneBins:
            genes = geneBins[pos1]
            # print "Found gene ", genes
            if pos2 in enhBins:
                enh = enhBins[pos2]
                pair = "\t%s\t%s\n"  % (",".join(genes),enh)
                if pair not in foundPairs:
                    resFile.write(line + pair )
                    found += 1
                    foundPairs.add(pair)
                #print "Found enhancer ", enh
                #for g in genes:
                #    novelConnections.append( [ enh, g ] )
            

        if pos2 in geneBins:
            genes = geneBins[pos2]
            #print "Found gene ", genes
            if pos1 in enhBins:
                enh = enhBins[pos1]
                #print "Found enhancer", enh
                pair = "\t%s\t%s\n"  % (enh, ",".join(genes))
                if pair not in foundPairs:
                    resFile.write(line + pair )
                    found += 1
                    foundPairs.add(pair)
                #for g in genes:
                #    novelConnections.append( [ enh, g ] )
                #break
        

    k += 1

print "Found %d enhancer-gene connections  out of %d total"  % (found, k)




    







