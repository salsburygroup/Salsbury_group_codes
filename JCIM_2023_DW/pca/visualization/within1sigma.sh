path=/deac/salsburyGrp/wud18/md/TM/pca/visualization
#for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops betaSheet1 connection;do
for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 regulatoryLoops connection 30sLoop;do
    cd ${i}
    num=$(ls | grep cluster | wc -l)
    for ((j=0; j<${num}; j++));do
	cd cluster$j
	python /home/wud18/python/within1sigma.py -s ${path}/${i}/cluster${j}/rep.pdb -t ${path}/${i}/cluster${j}/all.dcd -o1 ${path}/${i}/cluster${j}/Representative.pdb -o2 ${path}/${i}/cluster${j}/within1sigma.dcd
	cd ..
    done
    cd ..
done
