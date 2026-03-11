import json
import io
import os
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from numpy import linalg as la
from PIL import Image, ImageDraw, ImageSequence, ImageFont

with open('plot-dbh.json','r') as indata:
	data = json.load(indata)
	print(data.keys())

#print(data['hop'][0])
def pmd_spaghetti_plt(cutoff,meas_1,meas_2,ax,para,hop):
	x,y,x2,y2=[],[],[],[]
	n = 1
	print('length para=')
	print(len(para))
	for i in para:
		for j in range(len(i)): 
			x.append(i[j][meas_1])
			y.append(i[j][meas_2])
		ax.plot(x,y,marker='o',markersize=0.,linewidth=0.5,alpha=0.3)
		x=list()
		y=list()
		n = n + 1 

	hop_l = len(hop)
	for k in range(len(para)):
		u = 1
		if len(hop[k]) > 0: 
			n = -u
			j = hop[k][n]
			j = j - 1
			if j < cutoff:
				print(para[k][j]) 
				x2.append(para[k][j][meas_1])
				y2.append(para[k][j][meas_2])
			else:
				u = u + 1
		else:
			continue
	ax.scatter(x2,y2,color='black',marker='o',s=3,zorder=3)

fig = plt.figure()
ax = fig.add_subplot(111)


#CN-CN
fig = plt.figure(figsize=(6, 6))  # Adjusting figure size to be square
ax = fig.add_subplot(111, aspect='equal')  # Setting aspect ratio to equal
pmd_spaghetti_plt(2001,1,2,ax,data['data'],data['hop'])
ax.set_xlim(1.2,5.3)
ax.set_ylim(1.2,5.3)
plt.xticks(fontsize=15, weight='bold')
plt.yticks(fontsize=15, weight='bold')
ax.set_xlabel(r'C-N (Å)',color='blue',fontsize=14,labelpad=8)
ax.set_ylabel(r'C-N (Å)',color='green',fontsize=14,labelpad=8)
plt.savefig('DBH-CN-CN.png',dpi=1200,bbox_inches='tight')
plt.show()


#CC-Dihedral 
#pmd_spaghetti_plt(2001,4,0,ax,data['data'],data['hop'])
#ax.set_xlim(0,180)
#ax.set_ylim(1.2,2.7)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('DBH-CC-Dihedral.png',dpi=1400)
#ax.set_xlabel(r'Dihedral H-C-C-C (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'C-C (Å)',color='red',fontsize=14,labelpad=8)

##parameters: B 3 4 B 3 6 B 4 7 B 6 7 D 15 2 3 5 D 14 1 4 5 D 4 1 2 3

#C3-N D 15 2 3 5 
#pmd_spaghetti_plt(2001,4,1,ax,data['data'],data['hop'])
#ax.set_xlim(0,180)
#ax.set_ylim(1.2,2.7)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('C3N-Dihedral-Inversion.png',dpi=1400)
#ax.set_xlabel(r'Dihedral 15 2 3 5 (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'Dihedral 14 1 4 5 (°)',color='red',fontsize=14,labelpad=8)

#C4-N D 14 1 4 5
#pmd_spaghetti_plt(2001,5,2,ax,data['data'],data['hop'])
#ax.set_xlim(0,185)
#ax.set_ylim(1.2,5.3)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('C4N-Dihedral-Inversion.png',dpi=1400)
#ax.set_xlabel(r'Dihedral 15 2 3 5 (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'Dihedral 14 1 4 5 (°)',color='red',fontsize=14,labelpad=8)

#D 15 2 3 5 D 14 1 4 5
#pmd_spaghetti_plt(2001,4,5,ax,data['data'],data['hop'])
#ax.set_xlim(0,180)
#ax.set_ylim(1.2,2.7)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('DBH-Dihedral-Inversion.png',dpi=1400)
#ax.set_xlabel(r'Dihedral 15 2 3 5 (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'Dihedral 14 1 4 5 (°)',color='red',fontsize=14,labelpad=8)

#D 15 2 3 5 D 4 1 2 3
#pmd_spaghetti_plt(2001,4,6,ax,data['data'],data['hop'])
#ax.set_xlim(0,180)
#ax.set_ylim(0,75)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('DBH-Dihedral-Inversion.png',dpi=1400)
#ax.set_xlabel(r'Dihedral 15 2 3 5 (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'Dihedral 4 1 2 3 (°)',color='red',fontsize=14,labelpad=8)

#CC-D 4 1 2 3
#pmd_spaghetti_plt(2001,6,0,ax,data['data'],data['hop'])
#ax.set_xlim(0,180)
#ax.set_ylim(1.2,2.7)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('DBH-CC-Dihedral.png',dpi=1400)
#ax.set_xlabel(r'Dihedral H-C-C-C (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'C-C (Å)',color='red',fontsize=14,labelpad=8)

