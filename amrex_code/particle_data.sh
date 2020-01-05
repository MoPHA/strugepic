#cat data.out | grep POS | cut -f 2 -d"[" | sed 's/,0]$//g' | sed 's/,/ /g' > data.data
cat data.out | grep POS | cut -f 2 -d "[" | cut -f1 -d "]" | sed 's/,/ /g' > data.data
cat data.out | grep VEL | cut -f 2 -d "[" | cut -f1 -d "]" | sed 's/,/ /g' >vel_data.data
cat data.out | grep ENERGY | cut -f 2 -d ":" > energy.data
cat data.out | grep MOMENTUM | cut -f 2 -d ":" > momentum.data
