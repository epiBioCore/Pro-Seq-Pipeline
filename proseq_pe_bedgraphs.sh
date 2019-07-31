#!/bin/bash

###indices
ribo="/local_storage/annotation_db/Mus_musculus/UCSC/mm10/Sequence/RiboHisatIndex/mm10_ribo"
genome="/local_storage/annotation_db/Mus_musculus/UCSC/mm10/Sequence/newHisat/genome"
spikein="/local_storage/annotation_db/Homo_sapiens/UCSC/hg19/Sequence/newHisat/genome"
chrInfo="mm10_chrominfo.txt"

##directories
fastq=Raw_fastq
trim=Trimmed_fastq
align=Alignment
bigwigs=BigWigs


if [ ! -d $bigwigs ]
then

        mkdir $bigwigs

fi



if [ ! -d $trim ]
then

	mkdir $trim
fi

if [ ! -d $align ]
then

        mkdir $align
fi



#for i in `ls $fastq/*fq.gz | cut -f2 -d "/" | cut -f1,2,3,4 -d '_' | uniq`
#do
#
#echo $i 
#
#read1=($fastq/${i}*1.fq.gz)
#
#read2=($fastq/${i}*2.fq.gz)
#
#
#
#java -jar /usr/share/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 14 -phred33 -trimlog $trim/${i}_log.txt \
#$read1 $read2 \
#$trim/${i}_R1_paired_at.fq.gz $trim/${i}_R1_Unpaired_at.fq.gz \
#$trim/${i}_R2_paired_at.fq.gz $trim/${i}_R2_Unpaired_at.fq.gz \
#ILLUMINACLIP:proseq_adapters.fa:2:30:10 MINLEN:10 2> $trim/${i}_trimStats.txt
#
#
#echo $i trim done!
#done
#
##done
#
#for f in $(ls $trim/*_paired_at.fastq.gz | cut -f2 -d "/" | cut -f1,2,3,4 -d '_' | uniq)
#do 
#
#	#echo $f| tee -a rib_alignstat alignstat spikein_stat
#
#	## align to spikein
#	hisat2 -p 10 -x $spikein --no-spliced-alignment --phred33 --un-conc-gz ${align}/${f}_%_no_spikein.fq.gz -1 ${trim}/${f}_R1_paired_at.fq.gz -2 ${trim}/${f}_R2_paired_at.fq.gz 2>> ${align}/${f}_spikein_alignstat.txt > /dev/null
#	##remove ribosomal
#
#	hisat2 -p 10 -x $ribo --no-spliced-alignment --phred33 --un-conc-gz ${align}/${f}_%_norib.fq.gz -1 ${align}/${f}_1_no_spikein.fq.gz -2 ${align}/${f}_1_no_spikein.fq.gz 2>> ${align}/${f}_rib_alignstat.txt > /dev/null
#
#	##align to genome
#
#	hisat2 -p 10 -x $genome --no-spliced-alignment --phred33 -1 ${align}/${f}_1_norib.fq.gz -2 ${align}/${f}_2_norib.fq.gz 2>> ${align}/${f}_alignstat.txt | grep -E 'NH:i:1$|^@'| grep -E 'XM:i:[012]|^@' | samtools view -b - > ${align}/${f}.bam
#	#in next command greps are for eliminating reads aligned to many places (NH) and/or with #more than 2 mismatches (XM):
#
#	samtools sort -@ 6 -m 1G -o ${align}/${f}.bam ${align}/${f}.bam
#
#	samtools index ${align}/${f}.bam
#
#	
#
#done
#



if [ ! -d $bigwigs ]
then

	mkdir $bigwigs

fi


for i in ${align}/*.bam
do

	sample=$(basename $i .bam)
	echo converting $sample to bed files	
        bamToBed -i $i -bedpe -mate1 | awk 'BEGIN {OFS="\t"} ($5 > 20) {print $0}' > tmp.bed

	echo splitting  $sample into 5prime and 3prime samples
	### make bed files for 5' and '3
	awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$9}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$9}' tmp.bed |\
	awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6=="+"?"-":"+"}' | gzip  > ${sample}_3prime.bed.gz 
       ## switching the strand of the read


	awk 'BEGIN{OFS="\t"} ($9 == "+") {print $4,$6-1,$6,$7,$8,$9}; ($9 == "-") {print $4,$5,$5+1,$7,$8,$9}' tmp.bed |\
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6=="+"?"-":"+"}' | gzip  > ${sample}_5prime.bed.gz 
	## switching the strand of the read

	
	##need to normalize by per million, this will change to the number of reads aligning to genebodies
	N=$(wc -l tmp.bed | awk '{print $1}')
	norm=$(echo "scale=6; 1000000/$N"|bc)

	echo converting to bedGraph
	### make bedgraphs for each strand
	genomeCoverageBed -i ${sample}_3prime.bed.gz -scale $norm -g $chrInfo -bg -strand + > ${sample}_3prime_plus.bedGraph
	genomeCoverageBed -i ${sample}_3prime.bed.gz -scale $norm -g $chrInfo -bg -strand - > ${sample}_3prime_minus_tmp.bedGraph

	## for minus strand bedgraphs "invert" coverage
	awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' ${sample}_3prime_minus_tmp.bedGraph  > ${sample}_3prime_minus.bedGraph
	rm ${sample}_3prime_minus_tmp.bedGraph	

	genomeCoverageBed -i ${sample}_5prime.bed.gz -scale $norm -g $chrInfo -bg -strand + > ${sample}_5prime_plus.bedGraph
        genomeCoverageBed -i ${sample}_5prime.bed.gz -scale $norm -g $chrInfo -bg -strand - > ${sample}_5prime_minus_tmp.bedGraph
	awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' ${sample}_5prime_minus_tmp.bedGraph  > ${sample}_5prime_minus.bedGraph
	rm ${sample}_5prime_minus_tmp.bedGraph 	

	echo sorting
	##sort
	LC_COLLATE=C sort -k1,1 -k2,2n ${sample}_3prime_plus.bedGraph > ${sample}_3prime_plus_sorted.bedGraph
	LC_COLLATE=C sort -k1,1 -k2,2n ${sample}_3prime_minus.bedGraph > ${sample}_3prime_minus_sorted.bedGraph

	LC_COLLATE=C sort -k1,1 -k2,2n ${sample}_5prime_plus.bedGraph > ${sample}_5prime_plus_sorted.bedGraph
        LC_COLLATE=C sort -k1,1 -k2,2n ${sample}_5prime_minus.bedGraph > ${sample}_5prime_minus_sorted.bedGraph

	echo converting to BigWig
	##bigWigs
	bedGraphToBigWig ${sample}_3prime_plus_sorted.bedGraph $chrInfo $bigwigs/${sample}_3prime_plus.bw
	bedGraphToBigWig ${sample}_3prime_minus_sorted.bedGraph $chrInfo $bigwigs/${sample}_3prime_minus.bw
	
	bedGraphToBigWig ${sample}_5prime_plus_sorted.bedGraph $chrInfo $bigwigs/${sample}_5prime_plus.bw
	bedGraphToBigWig ${sample}_5prime_minus_sorted.bedGraph $chrInfo $bigwigs/${sample}_5prime_minus.bw



 
	echo Done!
done

