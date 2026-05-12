#BSUB -q normal
#BSUB -J bulkcut2
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=12]"
#BSUB -n 12

# conda activate macs3

chromSize="/path/mm10_ucsc/mm10_chromSize"

tissue=$(basename `pwd`)
dname=$(dirname `pwd`)
protein=$(basename $dname)
ddname=$(dirname $dname)
type=$(basename $ddname)


IgG=/share/home/hanxiaoping/RAWDATA/LMNA/WXY_bulkcut/final/$type/IgG/$tissue/
sort -k1,1 -k2,2n $tissue.$protein.log2.20kb.bedgraph > $tissue.$protein.log2.20kb.sorted.bedgraph
bedGraphToBigWig $tissue.$protein.log2.20kb.sorted.bedgraph /path/mm10_ucsc/mm10_chromSize.txt $tissue.$protein.log2.20kb.bw 

# shuf

KOdir=/share/home/hanxiaoping/RAWDATA/LMNA/WXY_bulkcut/final/KO
KO=$KOdir/$protein/$tissue/bowtie2.fragments.bed
WT=./bowtie2.fragments.bed

min=`wc -l $KO $WT | awk '{print $1}' | sort -n | sed -n '1p'`
shuf $KO -n $min -o $KOdir/$protein/$tissue/bowtie2.fragments.shuf.bed
shuf $WT -n $min -o bowtie2.fragments.shuf.bed


#---peak---

broad
macs3 callpeak -f BED -t bowtie2.fragments.bed \
	-c $IgG/bowtie2.fragments.bed \
 	-g mm -n $tissue.input --keep-dup all \
 	--nomodel --shift 0 --extsize 200 --broad \
 	--outdir ./peak_input_all/

macs3 callpeak -f BED -t bowtie2.fragments.shuf.bed \
 	-B -g mm -n $tissue --keep-dup 1 --nolambda \
 	--nomodel --shift 0 --extsize 200 --broad \
 	--outdir ./peak_shuf/
cat ./peak_shuf/$tissue\_peaks.broadPeak | awk '$5 > 19{print $0}' > ./peak_shuf/$tissue\_peaks_20.broadPeak

# narrow
macs3 callpeak -f BED -t bowtie2.fragments.shuf.bed \
 	-B -g mm -n $tissue --keep-dup 1 --nolambda \
 	--nomodel --shift 0 --extsize 200  \
 	--outdir ./peak_shuf/

cat ./peak_shuf/$tissue\_peaks.narrowPeak | awk '$5 > 19{print $0}' > ./peak_shuf/$tissue\_peaks_20.narrowPeak

# diff peak

KOdir=/share/home/hanxiaoping/RAWDATA/LMNA/WXY_bulkcut/final/KO
d1=`wc -l $KOdir/$protein/$tissue/bowtie2.fragments.bed | awk '{print $1}' | sort -n | sed -n '1p'`
d2=`wc -l bowtie2.fragments.bed | awk '{print $1}' | sort -n | sed -n '1p'`

cd peak_shuf/
macs3 bdgdiff --t1 $KOdir/$protein/$tissue/peak_shuf/$tissue\_treat_pileup.bdg --c1 $KOdir/$protein/$tissue/peak_shuf/$tissue\_control_lambda.bdg \
--t2 $tissue\_treat_pileup.bdg --c2 $tissue\_control_lambda.bdg  \
-g 75 -l 150 -C 1 --o-prefix diff_ko_vs_wt

# # ---bam compare----

samtools sort -o  bowtie2.sorted.bam  bowtie2.mapped.bam 
samtools index -@ 9 bowtie2.sorted.bam && rm bowtie2.mapped.bam
bamCoverage -b bowtie2.sorted.bam -o $tissue.$protein.bin10.bw \
	--binSize 10 \
	--outFileFormat bigwig \
  --normalizeUsing CPM  

# for LMNB1 LMNA H3K9me3 H3K9me2
bamCompare -b1 bowtie2.sorted.bam -b2 $IgG/bowtie2.sorted.bam \
--binSize 20000 \
--scaleFactorsMethod None \
--normalizeUsing CPM \
--outFileFormat bedgraph \
--operation log2 \
-o $tissue.$protein.log2.20kb.bedgraph

# for H3K4me3 H3K27ac H3K27me3
bamCompare -b1 bowtie2.sorted.bam -b2 $IgG/bowtie2.sorted.bam \
--binSize 10 \
--scaleFactorsMethod None \
--normalizeUsing CPM \
--outFileFormat bigwig \
--operation log2 \
-o $tissue.$protein.log2.10bp.bw