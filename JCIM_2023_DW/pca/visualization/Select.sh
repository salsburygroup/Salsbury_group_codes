#for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops;do

for i in 170sLoop 180sLoop 220sLoop 60sLoop gammaLoop helix1 regulatoryLoops connection catalyticPocket 30sLoop;do
    cd ${i}
    ../select_loops.sh
    cd ..
done

for i in catalyticTriad exositeI exositeII;do
    cd ${i}
    ../select_residues.sh
    cd ..
done
