if [ ! -d Plot_timeseries ]; then
    mkdir Plot_timeseries
fi

for i in 170sLoop 220sLoop 60sLoop catalyticPocket exositeI exositeII gammaLoop helix1 helix2 connection betaSheet1 30sLoop;do
    cp protein_stride100_${i}/IMWKRescaled/Plot_timeseries.png Plot_timeseries/Plot_timeseries_${i}_AH.png
done

for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1 30sLoop;do
    cp protein_stride100_${i}/HDBSCAN/Plot_timeseries.png Plot_timeseries/Plot_timeseries_${i}_HD.png
done
