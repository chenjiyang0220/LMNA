#BSUB -q normal
#BSUB -J FINAL
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=9]"
#BSUB -n 9

col=$(basename `pwd`) #Pituitary
dropseq_root=/path/Drop-seq_tools-2.5.1

# samtools merge $type.$tissue.all.bam $(find . -name "*.bam")
${dropseq_root}/DigitalExpression -m 8g \
    I=./tag/star_gene_exon_tagged.bam \
    CELL_BARCODE_TAG=XC \
    MOLECULAR_BARCODE_TAG=XM \
    O=$col.dge.txt.gz \
    SUMMARY=$col.dge.summary.txt \
    NUM_CORE_BARCODES=20000 \
    LOCUS_FUNCTION_LIST=INTRONIC \
    TMP_DIR=.

