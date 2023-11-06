path=/deac/salsburyGrp/wud18/md/TM/clustering
for i in 170sLoop 220sLoop 60sLoop catalyticPocket exositeI exositeII gammaLoop helix1 helix2 connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/IMWKRescaled/clusters
    num=$(ls | grep cluster | wc -l)
    for ((j=0; j<${num}; j++));do
	cd cluster$j
	python /home/wud18/python/within1sigma.py -s ${path}/protein_stride100_${i}/IMWKRescaled/clusters/cluster${j}/rep.pdb -t ${path}/protein_stride100_${i}/IMWKRescaled/clusters/cluster${j}/all.dcd -o1 ${path}/protein_stride100_${i}/IMWKRescaled/clusters/cluster${j}/Representative.pdb -o2 ${path}/protein_stride100_${i}/IMWKRescaled/clusters/cluster${j}/within1sigma.dcd
	cd ..
    done
    cd ../../..
done

for i in 170sLoop 180sLoop 220sLoop 60sLoop catalyticPocket catalyticTriad exositeI exositeII gammaLoop helix1 helix2 regulatoryLoops connection betaSheet1 30sLoop;do
    cd protein_stride100_${i}/HDBSCAN/clusters
    num=$(ls | grep cluster | wc -l)
    for ((j=0; j<${num}; j++));do
	cd cluster$j
	python /home/wud18/python/within1sigma.py -s ${path}/protein_stride100_${i}/HDBSCAN/clusters/cluster${j}/rep.pdb -t ${path}/protein_stride100_${i}/HDBSCAN/clusters/cluster${j}/all.dcd -o1 ${path}/protein_stride100_${i}/HDBSCAN/clusters/cluster${j}/Representative.pdb -o2 ${path}/protein_stride100_${i}/HDBSCAN/clusters/cluster${j}/within1sigma.dcd
	cd ..
    done
    cd ../../..
done
