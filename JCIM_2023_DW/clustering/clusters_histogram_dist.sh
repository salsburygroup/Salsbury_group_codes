for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1 30sLoop;do
#for i in 60sLoop;do
    cd protein_stride100_${i}/IMWKRescaled
    python ../../clusters_histogram_dist.py -title $i -d timeseries.txt -m IMWKRescaled -sp 3 -tm 1ns
    cd ..
    cd HDBSCAN
    python ../../clusters_histogram_dist.py -title $i -d timeseries.txt -m HDBSCAN -sp 3 -tm 1ns
    cd ../..
done
