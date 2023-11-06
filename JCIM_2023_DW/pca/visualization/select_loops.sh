num=$(ls | grep cluster | wc -l)

for ((j=0; j<${num}; j++));do
    cd cluster$j
    python ~/python/select_dcd.py -s Representative.pdb -t within1sigma.dcd -n 200
    cd ..
done
