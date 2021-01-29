#! /bin/bash

module load anaconda/2.7
cd /home/jiangyao/project/env_python2
source ./bin/activate

path1="/home/jiangyao/project/simulation/APA_simu/"
path2="/home/jiangyao/project/RNASeqReadSimulator-master/src/"
path3="/home/jiangyao/project/RNASeqReadSimulator-master/demo/input/"
path4="${path1}otherfiles/"

depth=500000

for i in $(seq 372 1000);do

	python ${path2}gensimreads.py -e ${path1}expr${i}_1.txt -n $depth -o ${path4}res1.bed ${path1}ref.bed12
    cat ${path4}res1.bed | python ${path2}getseqfrombed.py - ${path3}reference.fa | ${path2}splitfasta.py -o ${path4}test1_fa

	python ${path2}gensimreads.py -e ${path1}expr${i}_2.txt -n $depth -o ${path4}res2.bed ${path1}ref.bed12
    cat ${path4}res2.bed | python ${path2}getseqfrombed.py - ${path3}reference.fa | ${path2}splitfasta.py -o ${path4}test2_fa

    /home/jiangyao/project/forindex/STAR-2.7.1a/bin/Linux_x86_64/STAR --runThreadN 16 --runMode alignReads --genomeDir /home/jiangyao/project/forindex --readFilesIn ${path4}test1_fa_1.fa ${path4}test1_fa_2.fa --outFileNamePrefix ${path4}/ --outReadsUnmapped Fastx
    samtools view -b ${path4}Aligned.out.sam > ${path4}condition${i}_1.bam
    samtools sort ${path4}condition${i}_1.bam > ${path4}condition${i}_1.sorted.bam
    samtools index ${path4}condition${i}_1.sorted.bam

    /home/jiangyao/project/forindex/STAR-2.7.1a/bin/Linux_x86_64/STAR --runThreadN 16 --runMode alignReads --genomeDir /home/jiangyao/project/forindex --readFilesIn ${path4}test2_fa_1.fa ${path4}test2_fa_2.fa --outFileNamePrefix ${path4}/ --outReadsUnmapped Fastx
    samtools view -b ${path4}Aligned.out.sam > ${path4}condition${i}_2.bam
    samtools sort ${path4}condition${i}_2.bam > ${path4}condition${i}_2.sorted.bam
    samtools index ${path4}condition${i}_2.sorted.bam

    bamCoverage -bs 5 -b ${path4}condition${i}_1.sorted.bam -o ${path1}condition${i}_1.bw
    bamCoverage -bs 5 -b ${path4}condition${i}_2.sorted.bam -o ${path1}condition${i}_2.bw
done
