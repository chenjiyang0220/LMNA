
cd `pwd`
# ---------- HOMER analysis example ------------
# findMotifsGenome
for j in LTR LINE SINE ;do
    for i in H3K27ac H3K27me3 H3K4me1 H3K4me3 ;do
        mkdir $j/$i;
        findMotifsGenome.pl ../histone_KOup_TE_family/$j/$i\_$j\_KOupTE.bed mm10 ./$j/$i -size given -len 8,10,12 
    done
done

# annotate peak
for i in "${celltypes[@]}";do
    annotatePeaks.pl "../TE_KOup_enhancer/${i}_KOup_enhancer.bed" mm10 -m ./motif1.motif -mbed ./${i}/annotate_peak/KOup_enhancer_$i\_sites.bed > ./${i}/annotate_peak/KOup_enhancer_$i\_anno.txt
done

