#BSUB -q normal
#BSUB -J align_0820
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=9]"
#BSUB -n 9

cutadapt -a CTGTCTCTTATACACA -A CTGTCTCTTATACACA -m 20 -O 12 -o out1.fastq -p out2.fastq *1.fq.gz *2.fq.gz

/share/home/guoguoji/tools/STAR-2.7.10a/bin/Linux_x86_64_static/STAR \
	--genomeDir /share/home/guoguoji/tools/mm10_ucsc/STAR/ \
	--readFilesIn  out1.fastq out2.fastq \
	--readFilesCommand cat \
	--outFileNamePrefix star \
	--outSAMtype BAM SortedByCoordinate \
   --outBAMsortingThreadN 10

sample=$(basename `pwd`)

featureCounts -T 10 -t exon -g gene_id -p \
	-a /share/home/guoguoji/tools/mm10_ucsc/STAR/mm10_ucsc.gtf \
	-o $sample.count.txt starAligned.sortedByCoord.out.bam

TEtranscripts -t starAligned.sortedByCoord.out.bam -c starAligned.sortedByCoord.out.bam --GTF /share/home/guoguoji/tools/TEtranscripts-2.2.3/mm10_ucsc.gtf --TE /share/home/guoguoji/tools/TEtranscripts-2.2.3/mm10_rmsk_TE.gtf --sortByPos
TElocal -b starAligned.sortedByCoord.out.bam --GTF /share/home/guoguoji/tools/TEtranscripts-2.2.3/mm10_ucsc.gtf --TE /share/home/guoguoji/tools/TElocal-1.1.2/mm10_rmsk_TE.gtf.locInd --sortByPos
#rm out1.fastq out2.fastq