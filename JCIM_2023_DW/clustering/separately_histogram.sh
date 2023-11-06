for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1 30sLoop;do
    cp protein_stride100_${i}/IMWKRescaled/separately_histogram.png separately_histogram/separately_histogram_${i}_AH.png
    cp protein_stride100_${i}/HDBSCAN/separately_histogram.png separately_histogram/separately_histogram_${i}_HD.png
done
