#for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops;do
#for i in 170sLoop 180sLoop 220sLoop 60sLoop gammaLoop helix1 helix2 connection betaSheet1 catalyticPocket;do
for i in 170sLoop 180sLoop 220sLoop 60sLoop gammaLoop helix1 connection 30sLoop;do
    cd ${i}
    ../image_crop_loops.sh
    cd ..
done

for i in regulatoryLoops;do
    cd ${i}
    ../image_crop_RegulatoryLoops.sh
    cd ..
done

for i in catalyticPocket;do
    cd ${i}
    ../image_crop_catalyticPocket.sh
    cd ..
done

for i in catalyticTriad exositeI exositeII;do
    cd ${i}
    ../image_crop_residues.sh
    cd ..
done
