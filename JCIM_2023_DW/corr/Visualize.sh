python ~/python/Connect_image.py -im thrombin/corr.png TM56/corr_TM56_thrombin.png TM456/corr_TM456_thrombin.png TM456_TM56/corr_TM456_TM56.png -r 2 -c 2 -o corr.png
python ~/python/Connect_image.py -im thrombin/corr.png TM56/corr_TM56_thrombin.png TM456/corr_TM456_thrombin.png -r 1 -c 3 -o corr2.png

python ~/python/Connect_image.py -im visualize/TM56/TM56_thrombin_0.5_ig1.png visualize/TM456/TM456_thrombin_0.5_ig1.png visualize/TM456_TM56/TM456_TM56_0.5_ig1.png -r 1 -c 3 -o corr_visualization0.5_TM456_TM56.png

python ~/python/Connect_image.py -im visualize/TM56/TM56_thrombin_0.5_ig1.png visualize/TM456/TM456_thrombin_0.5_ig1.png visualize/white.png -r 1 -c 3 -o corr_visualization0.5.png
