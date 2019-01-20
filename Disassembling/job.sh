#! /bin/bash


AR=(

"0" "1" "2" "3" "4" "5" "6" "7" "8" "11" "12" "13" "14" "18" "19" "20" "21" "24" "25" "26" "27" "29" "32" "33" "34" "37" "38" "39" "40" "41" "42" "43" "44" "47" "48" "50" "54" "55" "56" "57" "59" "62" "66" "74" "80" "88" "91" "158" "187" "223" "224" "226"  

)

for i in ${AR[*]};
do
        mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd 
	./Simp ../data/pc/voxels_seg_$i.pcd 0.01f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd
        mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
done

mv ../data/pc/done/voxels_diff.pcd ../data/pc/voxels_diff.pcd 
./Simp ../data/pc/voxels_diff.pcd 0.005f ../data/mesh_simp/voxels_diff.vtk ../data/mesh_simp/voxels_diff.obj ../data/pc_simp/voxels_diff.pcd
mv ../data/pc/voxels_diff.pcd ../data/pc/done/voxels_diff.pcd


for i in `seq 0 239`
do
	doIt=1
	for j in ${AR[*]};
	do
		if [ "$i" -eq "$j" ]
		then 
			doIt=0
		fi
	done
	if [ "$doIt" -eq 1 ]
	then
	        mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd 
		./Simp ../data/pc/voxels_seg_$i.pcd 0.001f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd
        	mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
	fi
done
