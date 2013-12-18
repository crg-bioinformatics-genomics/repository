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

# takes sequences with certain lengths (eg 50 750 50 3000)

# selects the lengths
more ./database/protein.txt | awk '(length($2)>'$1')&&(length($2)<'$2')' | sort -k 1 | awk '{a[NR]=$0} END{for(i=1;i<=NR;i++){if(a[i]!=a[i+1]){print a[i]}}}'| awk '($2!~/U/)||($2!~/X/)||($2!~/Z/)||($2!~/B/)' >     ./database/protein.$1.$2.txt
more ./database/rna.txt     | awk '(length($2)>'$3')&&(length($2)<'$4')' | sort -k 1 | awk '{a[NR]=$0} END{for(i=1;i<=NR;i++){if(a[i]!=a[i+1]){print a[i]}}}'| awk '($2!~/X/)||($2!~/Z/)||($2!~/B/)'            >     ./database/rna.$3.$4.txt

# finds the interactions
for j in `cat ./database/rna.$3.$4.txt      | awk '{print $1}'`;  do grep $j     ./database/interactions.rna.txt; done > ./tmp/sel.1.tmp; 
for j in `cat ./database/protein.$1.$2.txt  | awk '{print $1}'`;  do grep $j     ./tmp/sel.1.tmp; done > ./tmp/sel.2.tmp
sort -k 1 ./tmp/sel.2.tmp | awk '{a[NR]=$0} END{for(i=1;i<=NR;i++){if(a[i]!=a[i+1]){print a[i]}}}' > ./database/interactions.rna.$1.$2.txt

