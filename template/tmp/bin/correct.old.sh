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
params_path=./params


# normalizes
as=`awk '{for(i=1;i<=NF;i++){s=s+sqrt(($i)^2)^2; k++}} END{print s}' ./data/ff.$1-$2.mat.txt`
#
an=`awk '{for(i=1;i<=NF;i++){s=s+$i; k++}} END{print k}' ./data/ff.$1-$2.mat.txt `
bn=`awk '{for(i=1;i<=NF;i++){s=s+$i; k++}} END{print k}' ./results/mat.$1.$2.tr.txt`
#
bs=`awk '{for(i=1;i<=NF;i++){s=s+sqrt(($i)^2)^('$an'/'$bn'); k++}} END{print s}' ./results/mat.$1.$2.tr.txt `

# parceval's theoreme
awk '{for(i=1;i<=NF;i++){s=$i/((('$bs')^2)/('$as')^2); s=(s-0.0000148)/0.000069; t=(s+2.5)/(4.29+2.83); u=(t-0.5)*2; printf "%f\t", s} printf "\n"}' ./results/mat.$1.$2.tr.txt  > ./results/norm.$1.$2.tr.tmp
 
# min-max scaling 
mi=`awk '{for(i=1;i<=NF;i++){print $i}}'   ./results/norm.$1.$2.tr.tmp  | sort -n -k 1 -r | tail -1`
ma=`awk '{for(i=1;i<=NF;i++){print $i}}'   ./results/norm.$1.$2.tr.tmp  | sort -n -k 1    | tail -1`

se=`echo $mi $ma | awk 'BEGIN{s=0} ($1<-3)||($2>3){s=1} END{print s}'`

# calculates the score
cat          ./tmp/score.txt      > ./tmp/analysis.txt
V=`cat tmp/analysis.txt`
cat $params_path/negatives.txt    >> ./tmp/analysis.txt
P=`awk '(NR==1){s=$1} (NR>1){if($1<s){k++}} END{print k/(NR-1)}' tmp/analysis.txt`
echo $P > ./tmp/disc.tmp

# overloaded
if [ $se == 1 ] ; then
fi=`cat ./tmp/disc.tmp`
awk '{for(i=1;i<=NF;i++){s=3*'$fi'^3*(($i-"'$mi'")/("'$ma'"-"'$mi'")-0.5)*2;printf "%f\t", s}  printf "\n"}' ./results/norm.$1.$2.tr.tmp > ./results/norm.$1.$2.tr.txt
fi

# normally loaded
if [ $se == 0 ] ; then
fi=`cat ./tmp/disc.tmp`
cat ./results/norm.$1.$2.tr.tmp > ./results/norm.$1.$2.tr.txt
fi

# background
awk '{s=0; d=0;  for(i=2;i<=NF;i++){s=s+$i; a[i]=$i} for(i=2;i<=NF;i++){d=d+(a[i]-s/(NF-1))^2;}  print sqrt(d/(NF-1))}' ./database/protein.dat | awk '{s=s+$1} END{print s/NR}' > ./tmp/varp.tmp
awk '{s=0; d=0;  for(i=2;i<=NF;i++){s=s+$i; a[i]=$i} for(i=2;i<=NF;i++){d=d+(a[i]-s/(NF-1))^2;}  print sqrt(d/(NF-1))}' ./database/rna.dat     | awk '{s=s+$1} END{print s/NR}' > ./tmp/varr.tmp
se2=`paste ./tmp/varp.tmp ./tmp/varr.tmp | awk 'BEGIN{s=0} ($1<0.005)||($2<0.02){s=1} END{print s}'`

if [ $se2 == 1 ]; then
V=0;
fi


# plots the distribution
sh  plot.dist.gp $V


if [ $se2 == 1 ]; then
awk '{for(i=1;i<=NF;i++){printf "%4.2f\t", 0} printf "\n"}' ./results/norm.$1.$2.tr.tmp > ./results/norm.$1.$2.tr.txt
fi 

# plots the maps
sh plot.4.tr.gp $1.$2 results eps
