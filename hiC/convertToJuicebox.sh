# prepare Juicebox files


#ID=9EP1
#RESDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/9EP1/hicpro_result

#ID=7EP41
#RESDIR=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/tumors/7EP41/hicpro_result

#cd $RESDIR

export TMPDIR=/omics/odcf/analysis/OE0290_projects/Ependymoma/external_data/tmpData
echo $TMPDIR

# initial
#DATADIR=/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/tumors
# latest Capture
DATADIR=/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/CaptureHiC/per_sample


for ID  in `ls ${DATADIR}`
do
    echo $ID

    cd $DATADIR/$ID/hicpro_result
    ls hic_results

    #RESFILE=${ID}_allValidPairs.hic
    RESFILE=${ID}.allValidPairs.hic
    if [ ! -f "$RESFILE" ];
    then
        echo "Prepare JuiceBox file..."
        #~/tools/HiCpro/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i  hic_results/data/${ID}/${ID}_allValidPairs -j ~/tools/juicebox/juicer_tools_0.7.5.jar  -g ~/tools/HiCpro/HiC-Pro_2.9.0/annotation/chrom_hg19.sizes
        ~/tools/HiCpro/HiC-Pro_2.11.4/bin/utils/hicpro2juicebox.sh -i  hic_results/data/${ID}/${ID}.allValidPairs -j ~/tools/juicebox/juicer_tools_1.13.01.jar  -g ~/tools/HiCpro/HiC-Pro_2.9.0/annotation/chrom_hg19.sizes
    else
        echo "Result exists!"
    fi
   
    echo
    #break

done

# other samples

# cmd

#~/tools/HiCpro/HiC-Pro_2.9.0/bin/utils/hicpro2juicebox.sh -i  hic_results/data/${ID}/${ID}_allValidPairs -j ~/tools/juicebox/juicer_tools_0.7.5.jar  -g ~/tools/HiCpro/HiC-Pro_2.9.0/annotation/chrom_hg19.sizes

#cd /icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/astroscyte_cerebellum_rep2/hicpro_result_detailed
ID=AstroCb_rep2
# cmd


# NB file
# cd /icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/SK-N-DZ/hicpro_result
ID=SK_N_DZ
# cmd

#cd /icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/SK-N-DZ/hicpro_result
ID=SK_N_DZ
#

#cd /icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/externalData/SK-N-SH/hicpro_result
ID=SK_N_SH

echo "Done!"
 
