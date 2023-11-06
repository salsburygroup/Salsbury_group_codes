num=$(ls | grep cluster | wc -l)
if [ $num -gt 20 ];then
    num=20
fi

# 1 figure
if [ $num == 1 ] ; then
    mv image0_label.png image.png
fi

# 2 figures
if [ $num == 2 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png -r 1 -c 2 -o image.png
fi

# 3 figures
if [ $num == 3 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png  -r 1 -c 3 -o image.png
fi

# 4 figures
if [ $num == 4 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png -r 2 -c 2 -o image.png
fi

# 5 figures
if [ $num == 5 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png white_label.png -r 3 -c 2 -o image.png
fi

# 6 figures
if [ $num == 6 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png -r 3 -c 2 -o image.png
fi

# 7 figures
if [ $num == 7 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png white_label.png -r 4 -c 2 -o image.png
fi

# 8 figures
if [ $num == 8 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png -r 4 -c 2 -o image.png
fi

# 9 figures
if [ $num == 9 ] ; then
    python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png -r 3 -c 3 -o image.png
fi

# 10 figures
if [ $num == 10 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png white_label.png white_label.png -r 4 -c 3 -o image.png

fi

# 11 figures
if [ $num == 11 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png white_label.png -r 4 -c 3 -o image.png

fi

# 12 figures
if [ $num == 12 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png -r 4 -c 3 -o image.png

fi

# 13 figures
if [ $num == 13 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png white_label.png white_label.png -r 5 -c 3 -o image.png

fi

# 14 figures
if [ $num == 14 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png image13_label.png white_label.png -r 5 -c 3 -o image.png

fi

# 15 figures
if [ $num == 15 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png image13_label.png image14_label.png -r 5 -c 3 -o image.png

fi

# 16 figures
if [ $num == 16 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png image13_label.png image14_label.png image15_label.png -r 4 -c 4 -o image.png

fi

# 17 figures
if [ $num == 17 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png image13_label.png image14_label.png image15_label.png image16_label.png white_label.png -r 6 -c 3 -o image.png

fi

# 18 figures
if [ $num == 18 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png image13_label.png image14_label.png image15_label.png image16_label.png image17_label.png -r 6 -c 3 -o image.png

fi

# 19 figures
if [ $num == 19 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png image13_label.png image14_label.png image15_label.png image16_label.png image17_label.png image18_label.png white_label.png -r 5 -c 4 -o image.png

fi

# 20 figures
if [ $num == 20 ] ; then

python ~/python/Connect_image.py -im image0_label.png image1_label.png image2_label.png image3_label.png image4_label.png image5_label.png image6_label.png image7_label.png image8_label.png image9_label.png image10_label.png image11_label.png image12_label.png image13_label.png image14_label.png image15_label.png image16_label.png image17_label.png image18_label.png image19_label.png -r 5 -c 4 -o image.png

fi

rm *_label.png
