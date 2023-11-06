for i in catalyticTriad exositeI exositeII;do
    cd ${i}
    ../image_crop_residues.sh
    cd ..
done
