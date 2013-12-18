#!/bin/sh
export GDFONTPATH=$HOME/Library/Fonts:/Library/Fonts:/System/Library/Fonts

cp $2/norm.$1.tr.txt  $2/n.$1.tr.txt

# dp
sh ./bin/score.3.sh
awk '(NR==1){for(i=2;i<=NF;i++){a[i]=$i}} (NR>1){for(i=2;i<=NF;i++){s=s+$1*a[i]*$i}} END{print s}' ./tmp/rlc.tmp > ./tmp/score.2.txt
spec=`sh ./bin/processer.sh tmp/score.2.txt | awk '{print int($3*100)}'`
inter=$4;



# changes
h=`echo $inter $spec | awk '($1>50)&&($2>50){s=1} ($1<=50)&&($2>50){s=2} ($1>50)&&($2<=50){s=3} ($1<=50)&&($2<=50){s=4} END{print s}'`

if [ $h == 1 ]; then sentence1="Interaction Propensity  = $inter%";   sentence2="Discriminative  Power  = $spec%" ; fi
if [ $h == 2 ]; then sentence1="Interaction Propensity <= 50%";       sentence2="Discriminative  Power  = $spec%" ; fi
if [ $h == 3 ]; then sentence1="Interaction Propensity  = $inter%";   sentence2="Discriminative  Power <= 50%";     fi
if [ $h == 4 ]; then sentence1="Interaction Propensity <  33%";       sentence2="Discriminative  Power <  33%"   ;  fi
# changes

#
pl=`wc      ./results/n.$1.tr.txt | awk '{print $1}'`
rl=`head -1 ./results/n.$1.tr.txt | awk '{print NF}'`
#
gnuplot << EOF
#
set cbrange [-3.95:3.95]
set terminal jpeg large enhanced  size 1024,768
set output "./eps/matrix.jpg"
set label "Name: \"$5\"" at graph  0.0, graph  1.15
set label "$sentence1" at graph  0.0, graph  1.10
set label "$sentence2" at graph  0.0, graph  1.05
set multiplot
set size 1,1
set origin 0.1,0.15
set size 0.8,0.8
set xrange [1:$rl]
set yrange [1:$pl]
set border 0 lw 3
set pm3d map 
set palette model RGB defined (0 "blue", 100 "white", 110 "white", 200 "red") 
set ylabel "Protein Residue Index"
set xlabel "RNA Nucleotide Index"
set xtics nomirror 
set ytics nomirror 
set xtics out
set ytics out

splot "results/n.$1.tr.txt" matrix notitle


EOF
