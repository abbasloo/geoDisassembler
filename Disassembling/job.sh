#! /bin/bash


AR=(

"0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "42" "43" "45" "46" "47" "48" "49" "51" "52" "53" "55" "56" "62" "63" "65" "66" "68" "70" "72" "73" "75" "76" "78" "80" "84" "88" "89" "108" "113" "115" "184" "185" "187" "192" "199" "222"

)

for i in ${AR[*]};
do
        mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd 
	./Simp ../data/pc/voxels_seg_$i.pcd 0.02f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd 1 200
        mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
done

mv ../data/pc/done/voxels_diff.pcd ../data/pc/voxels_diff.pcd 
./Simp ../data/pc/voxels_diff.pcd 0.006f ../data/mesh_simp/voxels_diff.vtk ../data/mesh_simp/voxels_diff.obj ../data/pc_simp/voxels_diff.pcd 0 
mv ../data/pc/voxels_diff.pcd ../data/pc/done/voxels_diff.pcd


for i in `seq 0 268`
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
		./Simp ../data/pc/voxels_seg_$i.pcd 0.02f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd 1 200
        	mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
	fi
done
