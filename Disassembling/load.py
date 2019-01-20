import csv
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D 


##############################################################################
# descriptors

Feat = []
with open('../Segmentation/Feat.csv', 'rb') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=';')
     for row in spamreader:
         Feat_ = []
	 for i in range(len(row)):
	     Feat_ = np.append(Feat_, row[i].split(', '))
	 Feat__ = []
	 for i in range(len(Feat_)):
		Feat__ = np.append(Feat__, Feat_[i].split('('))
	 Feat___ = []
	 for i in range(len(Feat__)):
		Feat___ = np.append(Feat___, Feat__[i].split(')'))
         Feat.append(Feat___)

ID = 30
#plt.figure(np.str(ID))
this = Feat[ID]
l = len(this)
for i in range(2, l-1):
	this[i] = np.float(this[i])
this = this[2:l-1]
#plt.plot(this)
#plt.show()

##############################################################################
# connectivities

con = [[]]
with open('../Segmentation/SegCon.csv', 'rb') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=';')
     for row in spamreader:
	 con_ = []
	 for i in range(len(row)):
	     if row[i] != '':
	     	con_ = np.append(con_, np.int(row[i]))
         con.append(con_)	 

ID = 100
seg = ''
for i in range(len(con[ID])):
	if i == 0:
		seg += ' ../data/mesh/voxels_seg_' + np.str(np.int(con[ID][i])) + '.vtk'
	if i != 0:
		seg += ' ../data/pc/voxels_seg_' + np.str(np.int(con[ID][i])) + '.pcd'

#os.system('pcl_viewer ' + seg)

#############################################################################
# floors and roofs

ar0 = []
ar1 = []
ar2 = []
ar3 = []
ar4 = []
ar5 = []
ar6 = []
ar7 = []
ar8 = []
with open('../Segmentation/Seg.csv', 'rb') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=';')
     for row in spamreader:
	ar0 = np.append(ar0, row[0])
	ar1 = np.append(ar1, np.float(row[1]))
	ar2 = np.append(ar2, np.float(row[2]))
	ar3 = np.append(ar3, np.float(row[3]))	
	ar4 = np.append(ar4, np.float(row[4]))	
	ar5 = np.append(ar5, np.float(row[5]))	
	ar6 = np.append(ar6, np.float(row[6]))	
	ar7 = np.append(ar7, np.float(row[7]))
	ar8 = np.append(ar8, np.float(row[8]))	

#plt.figure(); plt.hist(ar1, bins=100)
#plt.figure(); plt.hist(ar2, bins=100)
#plt.figure(); plt.hist(ar3, bins=100)
#plt.figure(); plt.hist(ar7, bins=100)

"""
idx = ar8 > 500
idx = np.where(idx)[0]
fig = plt.figure(); 
ax = fig.add_subplot(111, projection='3d');
for i in range(len(idx)):
	n = con[idx[i]]
        for j in range(len(n)):
		if ar8[np.int(n[j])] > 500:
			ax.plot([ar1[idx[i]], ar1[np.int(n[j])]], [ar2[idx[i]], ar2[np.int(n[j])]], [ar3[idx[i]], ar3[np.int(n[j])]], color='red');
	ax.text(ar1[idx[i]], ar2[idx[i]], ar3[idx[i]], np.str(ar0[idx[i]]), color='blue')
plt.axis('off');
plt.show();
"""

idx = ar7 < 1e-4
idx_ = ar8 > 1000

idx_4 = ar4 < 1e-2
idx_4 = idx_4*1
idx_5 = ar5 < 1e-2
idx_5 = idx_5*1
idx_6 = ar6 < 1e-2
idx_6 = idx_6*1
idx_456 = idx_4 + idx_5 + idx_6
idx_456 = idx_456 <= 1  

idx = idx*idx_
idx = np.where(idx)
a = ar0[idx]

idx = ar2 > 1.0
idx_ = ar8 > 500

idx_4 = ar4 < 1e-2
idx_4 = idx_4*1
idx_5 = ar5 < 1e-2
idx_5 = idx_5*1
idx_6 = ar6 < 1e-2
idx_6 = idx_6*1
idx_456 = idx_4 + idx_5 + idx_6
idx_456 = idx_456 <= 1 

idx = idx*idx_*idx_456
idx = np.where(idx)
b = ar0[idx]

seg = ''
for i in range(len(a)):
	seg += ' ../data/pc/voxels_seg_' + np.str(a[i]) + '.pcd'
for i in range(len(b)):
	seg += ' ../data/pc/voxels_seg_' + np.str(b[i]) + '.pcd'
#os.system('pcl_viewer ' + seg)

##############################################################################
# exteriors

idx_4 = ar4 < 1e-3
idx_4 = idx_4*1
idx_5 = ar5 < 1e-3
idx_5 = idx_5*1
idx_6 = ar6 < 1e-3
idx_6 = idx_6*1
idx_456 = idx_4 + idx_5 + idx_6
idx_456 = idx_456 <= 1 

con = []
with open('../Segmentation/Ext.csv', 'rb') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=';')
     for row in spamreader:
	 con = np.append(con, np.int(row[0]))

seg = ''
d = ''
for i in range(len(con)):
	doIt = True
        #for j in range(len(a)):
	#	if np.int(a[j]) == np.int(con[i]):
	#		doIt = False
        if (ar8[np.int(con[i])] > 100 and idx_456[np.int(con[i])]):
		if doIt:
	        	seg += ' ../data/pc/voxels_seg_' + np.str(np.int(con[i])) + '.pcd'
                        d += '"' + np.str(np.int(con[i])) + '" '
print d
os.system('pcl_viewer ' + seg)

##############################################################################
# anything that is not in exteriors, interior walls

idx_4 = ar4 < 1e-4
idx_4 = idx_4*1
idx_5 = ar5 < 1e-4
idx_5 = idx_5*1
idx_6 = ar6 < 1e-4
idx_6 = idx_6*1
idx_456 = idx_4 + idx_5 + idx_6
idx_456 = idx_456 <= 1

seg = ''
for i in range(len(idx_456)):
        doIt = True
        for j in range(len(con)):
		if i == np.int(con[j]) and idx_456[i]:
			doIt = False
        if (ar8[i] > 10 and  doIt):
	        seg += ' ../data/pc/voxels_seg_' + np.str(i) + '.pcd'
#os.system('pcl_viewer ' + seg)

##############################################################################
# pillars

idx_4 = ar4 < 5e-5
idx_4 = idx_4*1
idx_5 = ar5 < 5e-5
idx_5 = idx_5*1
idx_6 = ar6 < 5e-5
idx_6 = idx_6*1
idx_456 = idx_4 + idx_5 + idx_6
idx_456 = idx_456 <= 1

seg = ''
for i in range(len(idx_456)):
        if (~idx_456[i]):
	        seg += ' ../data/pc/voxels_seg_' + np.str(i) + '.pcd'
#os.system('pcl_viewer ' + seg)

##############################################################################
# layout

idx_4 = ar4 < 1e-3
idx_4 = idx_4*1
idx_5 = ar5 < 1e-3
idx_5 = idx_5*1
idx_6 = ar6 < 1e-3
idx_6 = idx_6*1
idx_456 = idx_4 + idx_5 + idx_6
idx_456 = idx_456 <= 1 

con = []
with open('../Segmentation/Ext.csv', 'rb') as csvfile:
     spamreader = csv.reader(csvfile, delimiter=';')
     for row in spamreader:
	 con = np.append(con, np.int(row[0]))

seg = ''
for i in range(len(con)):
	doIt = True
        for j in range(len(a)):
		if np.int(a[j]) == np.int(con[i]):
			doIt = False
        for j in range(len(b)):
		if np.int(b[j]) == np.int(con[i]):
			doIt = False
        if (ar8[np.int(con[i])] > 100 and idx_456[np.int(con[i])]):
		if doIt:
	        	seg += ' ../data/pc/voxels_seg_' + np.str(np.int(con[i])) + '.pcd'
#os.system('pcl_viewer ' + seg)

