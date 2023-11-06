for i in 170sLoop 220sLoop 60sLoop catalyticPocket exositeI exositeII gammaLoop helix1 helix2 connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/IMWKRescaled/clusters
    ../../../visualize.sh
    cd ../../..
done


for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/HDBSCAN/clusters
    ../../../visualize.sh
    cd ../../..
done
