import sys
import os

# input file is hicBreakFinder output

infile = sys.argv[1]

if not os.path.exists(infile):
    print "Error! File not found", infile
    sys.exit(-1)

assert ".breaks.txt" in infile

#output format
#r1   x1         x2         chr2   y1         y2         color     comment
#chrX   85000000   89000000   chrX   85000000   89000000   0,255,0   My green region


print "Writing result..."

resfile = open( infile.replace(".txt", ".juicebox.txt"), "w")

resfile.write("r1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment\n")

idx = 1
for line in open(infile):
    items = line.strip().split("\t")
    if line.startswith("#"):
        continue

    #print items
    chr1 = items[1]
    start1 = items[2]  
    end1 =  items[3]

    chr2 = items[5]
    start2 = items[6]
    end2 = items[7]

    tId ="hicBreakFinder_%s_%d" % (items[-1],idx)
    idx += 1
    newLine = "%s\t%s\t%s\t%s\t%s\t%s\t0,255,0\t%s\n" % ( chr1, start1, end1, chr2, start2,end2,tId )
    resfile.write(newLine)



