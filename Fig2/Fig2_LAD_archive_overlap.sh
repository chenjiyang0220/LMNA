#!/bin/bash

#### LMNA archive
$base_dir=LMNA

####-----------LADs categories analysis-----------
#### KO BLAD overlap with WT A/BLAD

KO_BLAD=$base_dir/bulkCUT/LAD_calling/results/LMNB1/KO
WT_LAD=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/WT_ABLADs
dest=$base_dir/bulkCUT/LAD_calling/analysis/KO_region/BLAD_switch

for i in merge; do
    bedtools intersect -a $KO_BLAD/$i/BED_files_2states/$i\_LADs.bed -b $WT_LAD/$i\_onlyA.bed >  $dest/$i\_AtoB.bed;
    bedtools intersect -a $KO_BLAD/$i/BED_files_2states/$i\_LADs.bed -b $WT_LAD/$i\_all_nonLADs.bed >  $dest/$i\_NontoB.bed;
    bedtools intersect -a $KO_BLAD/$i/BED_files_2states/$i\_nonLADs.bed -b $WT_LAD/$i\_A_B.bed >  $dest/$i\_ABtoNon.bed;
    bedtools intersect -a $KO_BLAD/$i/BED_files_2states/$i\_nonLADs.bed -b $WT_LAD/$i\_onlyB.bed >  $dest/$i\_BtoNon.bed;
done

#### LMNA WT A/B only LADs & overlap

ALAD=$base_dir/bulkCUT/LAD_calling/results/LMNA/WT
BLAD=$base_dir/bulkCUT/LAD_calling/results/LMNB1/WT
dest=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/WT_ABLADs

for i in merge ;do
    bedtools multiinter -names ALAD BLAD -i $ALAD/$i/BED_files_2states/$i\_LADs.bed $BLAD/$i/BED_files_2states/$i\_LADs.bed > $dest/$i\_WT_ABLADs.bed
done

####----------- WT histone mark region analysis (all antibody) ---------------------
LA=$base_dir/bulkCUT/LAD_calling/results/LMNA/WT
H3K9me2=$base_dir/bulkCUT/KDD_calling/results/H3K9me2/WT
H3K9me3=$base_dir/bulkCUT/KDD_calling/results/H3K9me3/WT
H3K4me1=$base_dir/bulkCUT/peak/WT/H3K4me1
H3K4me3=$base_dir/bulkCUT/peak/WT/H3K4me3
H3K27me3=$base_dir/bulkCUT/peak/WT/H3K27me3
LA=$base_dir/bulkCUT/LAD_calling/results/LMNA/WT
H3K9me2=$base_dir/bulkCUT/KDD_calling/results/H3K9me2/KO
H3K9me3=$base_dir/bulkCUT/KDD_calling/results/H3K9me3/KO

dest=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/active_hist

dest=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/each
for i in AdrenalGland  BoneMarrow  Cerebellum  Eye  Gonad  Kidney    Liver  LymphNode   Muscle   Pituitary  \
SmallIntestine  Spleen   Testicle  Vascular Bladder       BrainStem   Cerebrum    Fat  Heart \
LargeIntestine  Lung   MammaryGland  Pancreas  Skin       SpinalCord      Stomach  Thymus  ;do
    bedtools multiinter -names LA H3K9me2 H3K9me3 H3K4me1 H3K4me3 H3K27me3 \
        -i $LA/$i/BED_files_2states/$i\_LADs.bed \
        $H3K9me2/$i/BED_files_2states/$i\_KDDs.bed $H3K9me3/$i/BED_files_2states/$i\_KDDs.bed \
        $H3K4me1/$i/$i\_peaks_20.broadPeak $H3K4me3/$i/$i\_peaks_20.narrowPeak \
        $H3K27me3/$i/$i\_peaks_20.broadPeak \
        > $dest/$i.each.antibody.bed
done

###------------ WT ALADs overlap with histone mark region ---------------

dest=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/active_hist
LA=$base_dir/bulkCUT/LAD_calling/results/LMNA/WT
for i in AdrenalGland  BoneMarrow  Cerebellum  Eye  Gonad  Kidney    Liver  LymphNode   Muscle   Pituitary  \
SmallIntestine  Spleen   Testicle  Vascular Bladder       BrainStem   Cerebrum    Fat  Heart \
LargeIntestine  Lung   MammaryGland  Pancreas  Skin       SpinalCord      Stomach  Thymus ;do
    bedtools intersect -a $dest/$i.K4.region.bed -b $LA/$i/BED_files_2states/$i\_LADs.bed -wa > $dest/$i.K4.LA.ol.bed
done

## T1/T2 LADs identify(Fig2B 2SH)
dest=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/K9me2_3
LA=$base_dir/bulkCUT/LAD_calling/results/LMNA/WT
for i in AdrenalGland  BoneMarrow  Cerebellum  Eye  Gonad  Kidney    Liver  LymphNode   Muscle   Pituitary  \
SmallIntestine  Spleen   Testicle  Vascular Bladder       BrainStem   Cerebrum    Fat  Heart  \
LargeIntestine  Lung   MammaryGland  Pancreas  Skin       SpinalCord      Stomach  Thymus ;do
    bedtools intersect -a $LA/$i/BED_files_2states/$i\_LADs.bed -b $dest/$i.K9.region.bed  > $dest/$i.K9.LA.bed
done

###-------------C-ALAD analysis---------------

dest=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/CLAD
ALAD=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/ALAD

cd $base_dir/bulkCUT/LAD_calling/analysis/WT_region/ALAD
bed_files=($(find . -maxdepth 1 -name "*_LADs.bed" | sort))
names=()
for file in "${bed_files[@]}"; do
    base_name=$(basename "$file" .bed)
    name_part=${base_name%%_*}  
    names+=("$name_part")
done

names_str=$(IFS=' '; echo "${names[*]}")
echo "Bed files: ${bed_files[@]}"
echo "Names: ${names[@]}"

bedtools multiinter -i "${bed_files[@]}" -names $names_str > $dest/WT_cALAD_multiinter_sub.bed



#####------------ switched BLAD region overlap with K9me2/3 -------------------

WT_ABLADs=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/WT_ABLADs
KO_BLADs=$base_dir/bulkCUT/LAD_calling/results/LMNB1/KO
dest=$base_dir/bulkCUT/LAD_calling/analysis/KO_region/BLAD_switch_new

for tissue in merge ;do
    bedtools intersect -a $WT_ABLADs/$tissue\_AB_merge.bed -b $KO_BLADs/$tissue/BED_files_2states/$tissue\_nonLADs.bed > $dest/$tissue\_LADtoNon.bed
    bedtools intersect -a $WT_ABLADs/$tissue\_all_nonLADs.bed -b $KO_BLADs/$tissue/BED_files_2states/$tissue\_nonLADs.bed > $dest/$tissue\_NontoNon.bed
    bedtools intersect -a $WT_ABLADs/$tissue\_all_nonLADs.bed -b $KO_BLADs/$tissue/BED_files_2states/$tissue\_LADs.bed > $dest/$tissue\_NontoB.bed
    bedtools intersect -a $WT_ABLADs/$tissue\_AB_merge.bed -b $KO_BLADs/$tissue/BED_files_2states/$tissue\_LADs.bed > $dest/$tissue\_LADtoB.bed
done


H3K9me2=$base_dir/bulkCUT/KDD_calling/results/H3K9me2/KO
H3K9me3=$base_dir/bulkCUT/KDD_calling/results/H3K9me3/KO
H3K27me3=$base_dir/bulkCUT/peak/KO/H3K27me3
LMNB1=$base_dir/bulkCUT/LAD_calling/results/LMNB1/KO

switch_region=$base_dir/bulkCUT/LAD_calling/analysis/KO_region/Bswitch_KO_region

for tissue in AdrenalGland  BoneMarrow  Cerebellum  Eye  Gonad  Kidney    Liver  LymphNode   Muscle   Pituitary  \
SmallIntestine  Spleen   Testicle  Vascular Bladder       BrainStem   Cerebrum    Fat  Heart \
LargeIntestine  Lung   MammaryGland  Pancreas  Skin       SpinalCord      Stomach  Thymus ;do
   bedtools intersect -a $dest/$tissue\_LADtoNon.bed -b $H3K9me2/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_LADtoNon_H3K9me2.bed
   bedtools intersect -a $dest/$tissue\_LADtoNon.bed -b $H3K9me3/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_LADtoNon_H3K9me3.bed
   bedtools intersect -a $dest/$tissue\_LADtoNon.bed -b $H3K27me3/$tissue/$tissue\_peaks_20.broadPeak > $switch_region/$tissue\_LADtoNon_H3K27me3.bed
   bedtools intersect -a $dest/$tissue\_LADtoNon.bed -b $LMNB1/$tissue/BED_files_2states/$tissue\_LADs.bed > $switch_region/$tissue\_LADtoNon_LMNB1.bed

   bedtools intersect -a $dest/$tissue\_NontoNon.bed -b $H3K9me2/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_NontoNon_H3K9me2.bed
   bedtools intersect -a $dest/$tissue\_NontoNon.bed -b $H3K9me3/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_NontoNon_H3K9me3.bed
   bedtools intersect -a $dest/$tissue\_NontoNon.bed -b $H3K27me3/$tissue/$tissue\_peaks_20.broadPeak > $switch_region/$tissue\_NontoNon_H3K27me3.bed
   bedtools intersect -a $dest/$tissue\_NontoNon.bed -b $LMNB1/$tissue/BED_files_2states/$tissue\_LADs.bed > $switch_region/$tissue\_NontoNon_LMNB1.bed

   bedtools intersect -a  $dest/$tissue\_NontoB.bed -b $H3K9me2/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_NontoB_H3K9me2.bed
   bedtools intersect -a  $dest/$tissue\_NontoB.bed -b $H3K9me3/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_NontoB_H3K9me3.bed
   bedtools intersect -a  $dest/$tissue\_NontoB.bed -b $H3K27me3/$tissue/$tissue\_peaks_20.broadPeak > $switch_region/$tissue\_NontoB_H3K27me3.bed
   bedtools intersect -a  $dest/$tissue\_NontoB.bed -b $LMNB1/$tissue/BED_files_2states/$tissue\_LADs.bed > $switch_region/$tissue\_NontoB_LMNB1.bed

   bedtools intersect -a  $dest/$tissue\_LADtoB.bed -b $H3K9me2/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_LADtoB_H3K9me2.bed
   bedtools intersect -a  $dest/$tissue\_LADtoB.bed -b $H3K9me3/$tissue/BED_files_2states/$tissue\_KDDs.bed > $switch_region/$tissue\_LADtoB_H3K9me3.bed
   bedtools intersect -a  $dest/$tissue\_LADtoB.bed -b $H3K27me3/$tissue/$tissue\_peaks_20.broadPeak > $switch_region/$tissue\_LADtoB_H3K27me3.bed
   bedtools intersect -a  $dest/$tissue\_LADtoB.bed -b $LMNB1/$tissue/BED_files_2states/$tissue\_LADs.bed > $switch_region/$tissue\_LADtoB_LMNB1.bed

done

#####------------ C/F/nonLAD average CUT&Tag signal using bigWigAverageOverBed(Fig2G)-----------

CLADbed3=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/CLAD/WT_cALD.bed
FLADbed3=$base_dir/bulkCUT/LAD_calling/analysis/WT_region/FLAD
nonLADbed3=$base_dir/bulkCUT/LAD_calling/results/LMNA/WT
output=$base_dir/bulkCUT/analysis/LADs_feature/CLAD/CFLAD_Average

for mark in $marks;do
        bwdir=$base_dir/bulkCUT/bw/$mark/WT
        for tissue in $tissues;do
            bw=$bwdir/$tissue.$mark.log2.10bp.bw
            bigWigAverageOverBed "$bw" "$bed4dir/${tissue}_FLAD.bed" "$output/${mark}_FLAD_${tissue}.output.tab"
            bigWigAverageOverBed "$bw" "$bed4dir/WT_cALD.bed" "$output/${mark}_CLAD_${tissue}.output.tab"
            bigWigAverageOverBed "$bw" "$bed4dir/${tissue}_nonLADs.bed" "$output/${mark}_nonLAD_${tissue}.output.tab"
        done
done




