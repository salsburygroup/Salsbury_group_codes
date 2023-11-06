num=$(ls | grep cluster | wc -l)

for ((j=0; j<${num}; j++));do
    line=$((${j}+1))p
    dir=$(du -sh cluster* | sort -hr | sed 's/cluster/\ncluster/g' | sed '/K/d' | sed '/M/d' | sed '/G/d' | sed -n $line)
    mv ${dir} clusterN${j}
done

for ((j=0; j<${num}; j++));do
    mv clusterN${j} cluster${j}
done
