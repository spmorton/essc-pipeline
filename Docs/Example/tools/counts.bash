#!/bin/bash
totals=0
while read i; do
  #x=$(find $i -iname *.pot | wc)
  printf "%30s `find $i -iname *.pot | wc -l`\n" $i
  #set $thisRun
  #printf -v x '%d\n' $1 2>/dev/null
  #totals=$((totals+x))
done <seqs.txt


