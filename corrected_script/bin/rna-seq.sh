file1=$1
file2=$2
pre=`basename $file1 |cut -d '_' -f 1-4`;


#if read length < 45bp, segment length = half of read length, segment-mismatches = 0 or 1
mkdir $pre
#gunzip -f -c $file1 >$pre.1.fq
#gunzip -f -c $file2 >$pre.2.fq

tophat --bowtie1 -o $pre -G osativa7_newpol4ab.gtf --no-coverage-search -p 16 --library-type fr-firststrand  /qiyijunlab/00.processed_data/05.zhangbo/xl/geo_analysis/ref/ricev7/ricev7 $file1 $file2

#rm -rf $pre.1.fq
#rm -rf $pre.2.fq
