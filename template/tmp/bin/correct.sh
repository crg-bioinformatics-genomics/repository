

#***************************************************************************
#                             catRAPID 
#                             -------------------
#    begin                : Nov 2010
#    copyright            : (C) 2010 by G.G.Tartaglia
#    author               : : Tartaglia
#    date                 : :2010-11-25
#    id                   : caRAPID protein-RNA interactions
#    email                : gian.tartaglia@crg.es
#**************************************************************************/
#/***************************************************************************
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU General Public License as published by  *
# *   the Free Software Foundation; either version 2 of the License, or     *
# *   (at your option) any later version.                                   *
# *                                                                         *
# ***************************************************************************/
#!/bin/bash
#bin_path=`which start.modified.sh`
#bin_path=`dirname $bin_path`
params_path="./params"


# normalizes
cp ./results/mat.$1.$2.tr.txt   ./results/norm.$1.$2.tr.tmp
 
# min-max scaling 
mi=`awk '{for(i=1;i<=NF;i++){print $i}}'   ./results/norm.$1.$2.tr.tmp  | sort -n -k 1 -r | tail -1`
ma=`awk '{for(i=1;i<=NF;i++){print $i}}'   ./results/norm.$1.$2.tr.tmp  | sort -n -k 1    | tail -1`

# calculates the ip
cat          ./tmp/score.txt      > ./tmp/analysis.txt
V=`cat tmp/analysis.txt`

#  calculates the dp
sh ./bin/score.3.sh
awk '(NR==1){for(i=2;i<=NF;i++){a[i]=$i}} (NR>1){for(i=2;i<=NF;i++){s=s+$1*a[i]*$i}} END{print s}' ./tmp/rlc.tmp > ./tmp/score.2.txt
spec=`sh ./bin/processer.sh tmp/score.2.txt | awk '{print $3}'`

# takes the scores
sh ./bin/statistics.sh $1 $2 > ./tmp/disc.tmp
fi1=`cat ./tmp/disc.tmp | awk '('$spec'-$2<=0.7)&&($2!~/nan/)&&($2!~/inf/){print $2}      ('$spec'-$2>0.7)||($2~/nan/)||($2~/inf/){print '$spec'}'`
fi2=`cat ./tmp/disc.tmp | awk '      ($3>'$fi1')&&($3!~/nan/)&&($3!~/inf/){print -$3}         ($3<='$fi1')||($3~/nan/)||($3~/inf/){print '$fi1'}' `
fi3=`cat ./tmp/disc.tmp | awk '                   ($1!~/nan/)&&($1!~/inf/){print int($1*100)}               ($1~/nan/)||($1~/inf/){print int(100*'$fi1')}' `
#changes

cat ./results/norm.$1.$2.tr.tmp | awk '{for(i=1;i<=NF;i++){ u=(($i-"'$mi'")/("'$ma'"-"'$mi'")-0.5)*2*3; if(u<0){s=u*'$fi1'}  if(u>=0){s=u*'$fi2'}   printf "%f\t", s}  printf "\n"}' > ./results/norm.$1.$2.tr.txt

# poor frustration of sequences
awk '{s=0; d=0;  for(i=2;i<=NF;i++){s=s+$i; a[i]=$i} for(i=2;i<=NF;i++){d=d+(a[i]-s/(NF-1))^2;}  print sqrt(d/(NF-1))}' ./database/protein.dat | awk '{s=s+$1} END{print s/NR}' > ./tmp/varp.tmp
awk '{s=0; d=0;  for(i=2;i<=NF;i++){s=s+$i; a[i]=$i} for(i=2;i<=NF;i++){d=d+(a[i]-s/(NF-1))^2;}  print sqrt(d/(NF-1))}' ./database/rna.dat     | awk '{s=s+$1} END{print s/NR}' > ./tmp/varr.tmp


# plots the maps
sh bin/plot.4.tr.gp $1.$2 results eps $fi3 $3
