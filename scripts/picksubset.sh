#!/bin/bash

infile=$1
nqueries=$2
outfile=$3

cat $infile | shuf | head -n $nqueries > $outfile

