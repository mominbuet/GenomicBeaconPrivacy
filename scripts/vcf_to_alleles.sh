#!/bin/bash

extension=`echo $1| awk -F"." '{print $NF}'`

echo $extension

case $extension in

	vcf) exe="cat"
		;;
	gz) exe="zcat"
		;;
	bz2) exe="bzcat"
		;;
	*) echo "Unknown extension"
		kill -SIGQUIT $2
		;;
esac

if [[ $3 == "nopass" ]] ; then 
${exe} $1 | grep -v "#" | awk '{split($10,gt,":");if( index(gt[1],"1")>0 && index(gt[1],"0")>0 && length($5)==1 && length($4)==1){print $1,$2,$5}}' | sed 's/chr//g' > $2 ; 
else
${exe} $1 | grep -v "#" | grep "PASS" | awk '{split($10,gt,":");if( index(gt[1],"1")>0 && index(gt[1],"0")>0 && length($5)==1 && length($4)==1){print $1,$2,$5}}' | sed 's/chr//g' > $2 ; 
fi
