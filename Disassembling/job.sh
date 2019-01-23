#! /bin/bash


AR=(

"0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "16" "17" "18" "19" "20" "21" "23" "24" "25" "26" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "40" "41" "42" "43" "44" "45" "46" "48" "49" "51" "52" "53" "54" "56" "57" "58" "59" "61" "62" "63" "65" "70" "71" "73" "74" "77" "80" "93" "101" "105" "179" "181" "188" 

)

for i in ${AR[*]};
do
        mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd 
	./Simp ../data/pc/voxels_seg_$i.pcd 0.05f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd
        mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
done

mv ../data/pc/done/voxels_diff.pcd ../data/pc/voxels_diff.pcd 
./Simp ../data/pc/voxels_diff.pcd 0.005f ../data/mesh_simp/voxels_diff.vtk ../data/mesh_simp/voxels_diff.obj ../data/pc_simp/voxels_diff.pcd
mv ../data/pc/voxels_diff.pcd ../data/pc/done/voxels_diff.pcd


for i in `seq 0 262`
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
		./Simp ../data/pc/voxels_seg_$i.pcd 0.01f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd
        	mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
	fi
done
