cat data.out | grep POS | cut -f 2 -d"[" | sed 's/,0]$//g' | sed 's/,/ /g' > data.data
