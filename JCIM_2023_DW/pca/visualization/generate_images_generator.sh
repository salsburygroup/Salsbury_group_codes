cp generate_images_loops.vmd gammaLoop/generate_images.vmd
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 82 to 94' && mv generate_images.vmd 60sLoop
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 158 to 166' && mv generate_images.vmd helix1
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 205 to 212' && mv generate_images.vmd helix2
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 204 to 219' && mv generate_images.vmd 170sLoop
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 262 to 274' && mv generate_images.vmd 220sLoop
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 79 135 241' && mv generate_images.vmd catalyticTriad
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 57 98 104 106 109 142 143' && mv generate_images.vmd exositeI
python ~/python/Replace.py -fi generate_images_residues.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 125 134 281 284 288' && mv generate_images.vmd exositeII
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 225 to 239' && mv generate_images.vmd 180sLoop
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 167 to 170' && mv generate_images.vmd connection
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 171 to 176' && mv generate_images.vmd betaSheet1
python ~/python/Replace.py -fi generate_images_loops.vmd -fo generate_images.vmd -ci 'resid 182 to 190' -co 'resid 55 to 62' && mv generate_images.vmd 30sLoop
