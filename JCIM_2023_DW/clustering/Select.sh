for i in 170sLoop 180sLoop 220sLoop 60sLoop gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/HDBSCAN/clusters
    ../../../select_loops.sh
    cd ../../..
done

for i in 170sLoop 220sLoop 60sLoop gammaLoop helix1 helix2 connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/IMWKRescaled/clusters
    ../../../select_loops.sh
    cd ../../..
done



for i in catalyticPocket catalyticTriad exositeI exositeII;do
    cd protein_stride100_${i}/HDBSCAN/clusters
    ../../../select_residues.sh
    cd ../../..
done

for i in catalyticPocket exositeI exositeII;do
    cd protein_stride100_${i}/IMWKRescaled/clusters
    ../../../select_residues.sh
    cd ../../..
done
