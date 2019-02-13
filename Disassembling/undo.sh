#! /bin/bash

for i in `seq 0 268`
do
	mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd
done
