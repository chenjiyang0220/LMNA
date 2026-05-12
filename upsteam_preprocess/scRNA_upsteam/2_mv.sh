#BSUB -q normal
#BSUB -J mv
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=6]"
#BSUB -n 6

cd /share/home/hanxiaoping/RAWDATA/E200025630/L01

merge=merge
mkdir $merge

for i in LMNA_WT LMNA_KO;do
  # mkdir $merge/$i
  for j in Bladder Pancreas Heart Muscle;do
    mkdir $merge/$i/$j
    for k in {13..16};do
      mkdir $merge/$i/$j/$k
    done
  done
done
mkdir $merge/LMNA_KO/Vascular
mkdir $merge/LMNA_WT/WhiteFat
for i in {13..16};do
  mkdir $merge/LMNA_KO/Vascular/$i
  mkdir $merge/LMNA_WT/WhiteFat/$i
done

for i in {13..16};do
  for j in LMNA_WT LMNA_KO;do
    for k in Bladder Pancreas Heart Muscle;do
      mv process/$i/merge/$j/$k/*.bam $merge/$j/$k/$i
    done
  done

  mv process/$i/merge/LMNA_KO/Vascular/*.bam $merge/LMNA_KO/Vascular/$i
  mv process/$i/merge/LMNA_WT/WhiteFat/*.bam $merge/LMNA_WT/WhiteFat/$i
done
