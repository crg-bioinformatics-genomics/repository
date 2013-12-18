for g in `cat $1`; do 
awk '('$g'>$1){p++} END{printf "%4.2f\t%4.2f\t", '$g', p/NR}' params/negatives.txt;  
awk '('$g'<$1){p++} END{printf "%4.2f\n", p/NR}'              params/positives.txt;
done | awk '{f=0.33; s=($2*(1-$3)*f)/($3*(1-$2)*(1-f)+$2*(1-$3)*f); print $1,s, $2^0.5, 1-$3}'
