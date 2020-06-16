FNAME=$1
cat $FNAME | grep POS | cut -f 2 -d "[" | cut -f1 -d "]" | sed 's/,/ /g' > pos.data
cat $FNAME | grep VEL | cut -f 2 -d "[" | cut -f1 -d "]" | sed 's/,/ /g' >vel.data
cat $FNAME | grep ENERGY | cut -f 2 -d ":" > energy.data
