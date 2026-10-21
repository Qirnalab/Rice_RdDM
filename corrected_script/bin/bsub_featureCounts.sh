i=$1
m=$2
n=$3

name=`basename $i|cut -d '_' -f 1-3`

featureCounts -O -P -s 2 -p -d 10 -D 2000 -C -T 4 -a osativa7_newpol4ab.gtf -o $name $i $m $n
