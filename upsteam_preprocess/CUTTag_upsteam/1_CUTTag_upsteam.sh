#BSUB -q normal
#BSUB -J bulkcut
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=12]"
#BSUB -n 12

tmpdir=tmp
mkdir $tmpdir

# trim reverse ME
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -j 0 -m 20 -O 10 \
	-o out1.fastq -p out2.fastq \
	*1.fq.gz *2.fq.gz \
	--info-file=$tmpdir/cut.log > $tmpdir/trim_report.txt

ref=/path/mm10_ucsc/bowtie2/mm10
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant \
 	--phred33 -I 10 -X 700 -p 8 -x ${ref} \
 	-1 out1.fastq -2 out2.fastq -S bowtie2.bam &> $tmpdir/bowtie2.txt && rm $tmpdir/cut.log out1.fastq out2.fastq


# picardCMD="java -jar /path/Drop-seq_tools-2.5.1/3rdParty/picard/picard.jar"
# $picardCMD SortSam I=bowtie2.bam O=bowtie2.sorted.bam SORT_ORDER=coordinate && rm bowtie2.bam
# $picardCMD MarkDuplicates I=bowtie2.sorted.bam O=bowtie2.sorted.rmDup.bam \
#  	REMOVE_DUPLICATES=true METRICS_FILE=picard.rmDup.txt && rm bowtie2.sorted.bam

samtools view -F 0x04 bowtie2.bam | \
	awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' \
	| sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > $tmpdir/fragmentLen.txt

minQualityScore=2
samtools view -H bowtie2.bam > bowtie2.qs$minQualityScore.bam 
samtools view -q $minQualityScore bowtie2.bam >> bowtie2.qs$minQualityScore.bam && rm bowtie2.bam

samtools view -bS -F 0x04 bowtie2.qs$minQualityScore.bam > bowtie2.mapped.bam && rm bowtie2.qs$minQualityScore.bam


bedtools bamtobed -i bowtie2.mapped.bam -bedpe > bowtie2.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' bowtie2.bed > bowtie2.clean.bed && rm bowtie2.bed
cut -f 1,2,6 bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  > bowtie2.fragments.bed && rm bowtie2.clean.bed sorted.name.bam

samtools sort -o bowtie2.sorted.bam bowtie2.mapped.bam && rm bowtie2.mapped.bam
samtools index -@ 9 bowtie2.sorted.bam

#---bedgraph---
tissue=$(basename `pwd`)
dname=$(dirname `pwd`)
protein=$(basename $dname)
ddname=$(dirname $dname)
type=$(basename $ddname)

IgG=/share/home/hanxiaoping/RAWDATA/LMNA/WXY_bulkcut/final/$type/IgG/$tissue/

chromSize="/path/mm10_ucsc/mm10_chromSize.txt"
depth=`wc -l bowtie2.fragments.bed | awk '{print $1}' | sort -n | sed -n '1p'`
cpm=1000000
scaleFactor=$(echo "scale=4; $cpm/$depth"| bc)
bedtools genomecov -i bowtie2.fragments.bed -g $chromSize -bg -scale $scaleFactor   > $tissue\_$protein\_$type.bedgraph
bedGraphToBigWig bowtie2.fragments.bedgraph /path/mm10_ucsc/mm10_chromSize.txt $sample.bw



