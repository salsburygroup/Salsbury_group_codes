#for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1;do
for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 regulatoryLoops connection 30sLoop;do
    cd ${i}
    ../visualize.sh
    cd ..
done
