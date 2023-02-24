import sys

# script check if these genes were identified from previous study


# PFA
#inFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.EPD210FH_enh_contacts_closest.txt"

# RELA BT
#inFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_BT/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.RELA_BT_enh_contacts_closest.txt"

# RELA EP
#inFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/RELPOS_EP/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.RELA_EP_enh_contacts_closest.txt"

assert len(sys.argv) == 3

inFile = sys.argv[1]
grp = sys.argv[2]

print "Input data:", inFile
print "Group:", grp

if (grp != "PFA") and (grp != "RELA"):
    print "ERROR! Incorrect group:", grp
    sys.exit(-1)

# WARNING: multiple genes might be present!
contacts = []
for line in open(inFile):
    items = line.strip().split("\t")
    pos1 = "%s:%s" % (items[0], items[1])
    pos2 = "%s:%s" % (items[2], items[3])
    contacts.append( [ items[-1], items[-2], pos1, pos2 ] )

print "Loaded %d contacts" % len(contacts)

if grp == "PFA":
    # PFA  default
    checkFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/associateGeneEnhancer/PFA_specific_TE_SE_positive_correlated_genes_in_TAD.txt"

    # PFA filtered 
    #checkFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/binInfo/PFA_specific_TE_SE_positive_correlated_genes_in_TAD.filtered.txt"
else:
    assert grp == "RELA"
    # RELA default
    checkFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/associateGeneEnhancer/RELPOS_specific_TE_SE_positive_correlated_genes_in_TAD.txt"

print "Checking status from", checkFile

found = 0
total = 0

resFile =  open(inFile.replace(".txt",".EAG_confirmed.txt"), "w")


for line in open(checkFile):
    if line.startswith("peakid"):
        resFile.write(line.strip() + "\tbin1\tbin2\n")
        continue
    items = line.strip().split("\t")
    ec1 = items[0].split(":")
    ec2 = ec1[1].split("-")
    # adjust
    enh = "%s:%d-%s" % (ec1[0], int(ec2[0]) + 1, ec2[1])
    gene = items[3]
    total += 1
    #print enh,gene
    for c in contacts:
        check1 = (c[0] == enh) and ( gene in c[1].split(",") )
        check2 = (c[1] == enh) and ( gene in c[0].split(",")  )
        if check1 or check2:
            found += 1
            #print "Found!"
            resFile.write(line.strip() + "\t%s\t%s\n" % (c[2],c[3]))
            break       

print "Total found: %d out of %d" % (found,total)


