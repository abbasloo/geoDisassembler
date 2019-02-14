#! /bin/bash


#for i in `seq 0 268`
#do
#	mv ../data/pc/done/Hall3_seg_$i.pcd ../data/pc/Hall3_seg_$i.pcd
#done

#mv ../data/pc/done/Hall3_diff.pcd ../data/pc/Hall3_diff.pcd

AR=(

"0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "16" "17" "18" "19" "20" "21" "23" "24" "25" "26" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "40" "41" "42" "43" "44" "45" "46" "48" "49" "51" "52" "53" "54" "56" "57" "58" "59" "61" "62" "63" "65" "70" "71" "73" "74" "77" "80" "93" "101" "105" "179" "181" "188" 

)

for i in ${AR[*]};
do
        mv ../data/pc/done/Hall3_seg_$i.pcd ../data/pc/Hall3_seg_$i.pcd 
	./Simp ../data/pc/Hall3_seg_$i.pcd 0.02f ../data/mesh_simp/Hall3_seg_$i.vtk ../data/mesh_simp/Hall3_seg_$i.obj ../data/pc_simp/Hall3_seg_$i.pcd 0
        mv ../data/pc/Hall3_seg_$i.pcd ../data/pc/done/Hall3_seg_$i.pcd
done

AR=(

"1" "3" "4" "5" "9" "14" "15" "16" "23" "25" "27" "28" "31" "32" "33" "35" "38" "44" "46" "49" "50" "51" "53" "55" "56" "61" "66" "72" "78" "79" "83" "91" "113" "130" 

)

for i in ${AR[*]};
do
	mv ../data/pc/done/diff_seg_$i.pcd ../data/pc/diff_seg_$i.pcd 
	./Simp ../data/pc/diff_seg_$i.pcd 0.006f ../data/mesh_simp/diff_seg_$i.vtk ../data/mesh_simp/diff_seg_$i.obj ../data/pc_simp/diff_seg_$i.pcd 0
	mv ../data/pc/diff_seg_$i.pcd ../data/pc/done/diff_seg_$i.pcd
done


