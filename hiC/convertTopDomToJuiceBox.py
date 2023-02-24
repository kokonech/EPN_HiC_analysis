import sys
import os

# examples

infile = sys.argv[1]

if not os.path.exists(infile):
    print "Error! File not found", infile
    sys.exit(-1)

assert ".TopDom_result.domain" in infile

#format
#r1   x1         x2         chr2   y1         y2         color     comment
#chrX   85000000   89000000   chrX   85000000   89000000   0,255,0   My green region

print "Writing result..."

resfile = open( infile.replace(".domain", ".juicebox.txt"), "w")

resfile.write("r1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment\n")

chrVal = infile.split("/")[-1].split(".")[1]

assert chrVal.startswith("chr")

for line in open(infile):
    items = line.strip().split("\t")
    if line.startswith("#"):
        continue

    if items[5] != "domain":
        continue

    #print items
    idx1 = items[1] 
    idx2 = items[3]

    chrId = items[0]
    start = items[2]  
    end =  items[4]
    tId ="TopDom_%s_%s" % (idx1, idx2)
    newLine = "%s\t%s\t%s\t%s\t%s\t%s\t0,0,255\t%s\n" % ( chrId, start, end, chrId, start,end ,tId )
    resfile.write(newLine)



