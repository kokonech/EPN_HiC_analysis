import os.path
import sys
import fnmatch
import glob

dataDir="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/raw_data/HiC/tumors"

dataDir2 = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/raw_data/HiC/tumors_v2/C202SC18122892/raw_data"

sFile = open("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/2018_12_Ependymoma_Tumors_HiC_QC_Dixon.txt")

resultDir = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors"

template = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/config_template_tumor.txt"

print "Preparing input commands..."

scriptFile= open("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/submitHicPro.sh", "w")

# update with novel data
updateIds = [ "RB_94", "RB_95", "RB_96" ]


for line in sFile:
    if line.startswith("ID") or len(line.strip()) == 0:
        continue
   
    items = line.strip().split()
    sId = items[0].replace("-","_")
    sName = items[7]
    sGroup = items[8]
    print sId,sName,sGroup

    sDir = "%s/%s" % (dataDir, sId)

    if not os.path.exists(sDir):
        print "Error! Folder %s does not exist!" % sDir
        sys.exit(-1)

    r1Files = []
    for f in glob.glob("%s/*_1.fq.gz" % sDir):
        r1Files.append(f)

    if sId in updateIds:
        print "Include additional data..."
        sDir2 = "%s/%s" % (dataDir2, sId)
        for f in glob.glob("%s/*_1.fq.gz" % sDir2):
            r1Files.append(f)
    print r1Files

    resDir = resultDir + "/" +  sName
    print resDir
    if not os.path.exists(resDir):
        os.makedirs(resDir)
    
    resDirPath = resDir + "/hicpro_result"       
    if os.path.exists(resDirPath):
        print "HiC analysis was already performed, skipping..."
        continue


    configPath = "%s/config_%s.txt" % (resDir,sName)
    
    resFile = open( configPath  , "w" ) 
    
    for line in open(template):
        if line.startswith("LOGFILE"):
            line = "LOGFILE = hicpro_%s.log\n" % sName
            resFile.write(line)
            continue
        
        if line.startswith("JOB_NAME"):
            line = "JOB_NAME = %s_split\n" % sName
            resFile.write(line)
            continue

        resFile.write(line)

    inputDir = resDir + "/" + "input_data"      
    if not os.path.exists(inputDir):
        os.makedirs(inputDir)

    inSubDir = inputDir + "/" + sName
    if not os.path.exists(inSubDir):
        os.makedirs(inSubDir)

    for f1Path in r1Files:
        assert os.path.exists(f1Path)
        newLoci1 = "%s/%s" % (inSubDir, f1Path.split("/")[-1])
        newLoci1 = newLoci1.replace(".fq.gz",".fastq.gz")
        #print newLoci1
        if not os.path.exists(newLoci1):
            os.symlink(f1Path, newLoci1)
    
        f2Path = f1Path.replace("_1.fq.gz","_2.fq.gz")
        assert os.path.exists(f2Path)
        newLoci2 = "%s/%s" % (inSubDir, f2Path.split("/")[-1])
        newLoci2 = newLoci2.replace(".fq.gz",".fastq.gz")
        
        #print newLoci2
        if not os.path.exists(newLoci2):
            os.symlink(f2Path, newLoci2)

    print "Writing command for",sName
    scriptFile.write("cd %s\n" % resDir )
    scriptFile.write("~/tools/HiCpro/HiC-Pro_2.9.0/bin/HiC-Pro -c %s -i input_data -o hicpro_result -p\n\n" % (configPath.split("/")[-1] ) )
    #sys.exit(-1)
    

print "Finished!"

