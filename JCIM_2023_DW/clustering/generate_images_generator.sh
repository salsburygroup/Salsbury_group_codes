cp generate_images_loops.vmd protein_stride100_gammaLoop/IMWKRescaled/clusters/generate_images.vmd
#cp generate_images_regulatoryLoops.vmd protein_stride100_regulatoryLoops/IMWKRescaled/clusters/generate_images.vmd
cp generate_images_catalyticPockets.vmd protein_stride100_catalyticPocket/IMWKRescaled/clusters/generate_images.vmd
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 82 to 94' && mv generate_images.vmd protein_stride100_60sLoop/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 55 to 62' && mv generate_images.vmd protein_stride100_30sLoop/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 158 to 166' && mv generate_images.vmd protein_stride100_helix1/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 205 to 212' && mv generate_images.vmd protein_stride100_helix2/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 204 to 219' && mv generate_images.vmd protein_stride100_170sLoop/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 262 to 274' && mv generate_images.vmd protein_stride100_220sLoop/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 79 135 241' && mv generate_images.vmd protein_stride100_catalyticTriad/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 57 98 104 106 109 142 143' && mv generate_images.vmd protein_stride100_exositeI/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 125 134 281 284 288' && mv generate_images.vmd protein_stride100_exositeII/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 225 to 239' && mv generate_images.vmd protein_stride100_180sLoop/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 167 to 170' && mv generate_images.vmd protein_stride100_connection/IMWKRescaled/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 171 to 176' && mv generate_images.vmd protein_stride100_betaSheet1/IMWKRescaled/clusters

cp generate_images_loops.vmd protein_stride100_gammaLoop/HDBSCAN/clusters/generate_images.vmd
cp generate_images_regulatoryLoops.vmd protein_stride100_regulatoryLoops/HDBSCAN/clusters/generate_images.vmd
cp generate_images_catalyticPockets.vmd protein_stride100_catalyticPocket/HDBSCAN/clusters/generate_images.vmd
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 82 to 94' && mv generate_images.vmd protein_stride100_60sLoop/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 55 to 62' && mv generate_images.vmd protein_stride100_30sLoop/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 158 to 166' && mv generate_images.vmd protein_stride100_helix1/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 205 to 212' && mv generate_images.vmd protein_stride100_helix2/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 204 to 219' && mv generate_images.vmd protein_stride100_170sLoop/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 262 to 274' && mv generate_images.vmd protein_stride100_220sLoop/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 79 135 241' && mv generate_images.vmd protein_stride100_catalyticTriad/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 57 98 104 106 109 142 143' && mv generate_images.vmd protein_stride100_exositeI/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 125 134 281 284 288' && mv generate_images.vmd protein_stride100_exositeII/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 225 to 239' && mv generate_images.vmd protein_stride100_180sLoop/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 167 to 170' && mv generate_images.vmd protein_stride100_connection/HDBSCAN/clusters
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 171 to 176' && mv generate_images.vmd protein_stride100_betaSheet1/HDBSCAN/clusters
