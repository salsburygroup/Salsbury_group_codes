for i in catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 regulatoryLoops connection 30sLoop;do
    cd ${i}
    ../visualize.sh
    cd ..
done
