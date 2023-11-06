for i in 170sLoop 220sLoop 60sLoop catalyticPocket exositeI exositeII gammaLoop helix1 helix2 connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/IMWKRescaled/clusters
    ../../../connect_image.sh
    cp image.png ../../../image_${i}_AH.png
    cd ../../..
done


for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/HDBSCAN/clusters
    ../../../connect_image.sh
    cp image.png ../../../image_${i}_HD.png
    cd ../../..
done
