#! /bin/bash


for i in `seq 0 268`
do
	mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd
done

mv ../data/pc/done/voxels_diff.pcd ../data/pc/voxels_diff.pcd

AR=(

"3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "17" "18" "19" "20" "21" "22" "24" "25" "26" "27" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "42" "45" "46" "47" "48" "49" "51" "52" "53" "55" "56" "62" "63" "65" "66" "68" "70" "72" "73" "76" "78" "80" "84" "88" "89" "108" "113" "115" "184" "185" "187" "192" "199" "222" 

)

for i in ${AR[*]};
do
        mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd 
	./Simp ../data/pc/voxels_seg_$i.pcd 0.02f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd 1 800
        mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
done


AR=(

"0" "2" "43" "75" "246" "1"

)

for i in ${AR[*]};
do
        mv ../data/pc/done/voxels_seg_$i.pcd ../data/pc/voxels_seg_$i.pcd 
	./Simp ../data/pc/voxels_seg_$i.pcd 0.02f ../data/mesh_simp/voxels_seg_$i.vtk ../data/mesh_simp/voxels_seg_$i.obj ../data/pc_simp/voxels_seg_$i.pcd 1 200
        mv ../data/pc/voxels_seg_$i.pcd ../data/pc/done/voxels_seg_$i.pcd
done
