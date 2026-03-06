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
import seaborn as sns
from scipy.stats import gaussian_kde

with open('plot-dbh.json','r') as indata:
	data = json.load(indata)
	print(data.keys())

#print(data['hop'][0])
def pmd_spaghetti_plt(cutoff,meas_1,meas_2,ax,para,hop):
	x,y,x2,y2=[],[],[],[]
	n = 1
	for i in para:
		for j in range(len(i)): 
			x.append(i[j][meas_1])
			y.append(i[j][meas_2])
		#ax.plot(x,y,marker='o',markersize=0.,linewidth=0.5,alpha=0.3)
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

	classification = []
	for i in range(len(x2)):
		if x2[i]>110 and y2[i]>110:
			#ax.scatter(x2[i],y2[i],color='#460F62',marker='o',s=2,zorder=3)
			classification.append("Inverted-HP")
		elif x2[i]>110 or y2[i]>110:
			#ax.scatter(x2[i],y2[i],color='#1E9A8A',marker='o',s=2,zorder=3)
			classification.append("Partial-Inverted-HP")
		else:
			#ax.scatter(x2[i],y2[i],color='#E6D222',marker='o',s=2,zorder=3)
			classification.append("Retained-HP")
	# Count and print the number of each species
	from collections import Counter
	species_counts = Counter(classification)
	for species, count in species_counts.items():
		print(f"{species}: {count}")
	
	results_df=pd.DataFrame(list(zip(x2,y2,classification)),columns =['D 15 2 3 5','D 14 1 4 5','Type'])
	results_df.to_csv('classify_hopping.csv',index=False)

	#ax.scatter(x2,y2,color='black',marker='o',s=3,zorder=3)

	# Compute KDE Density
	xy = np.vstack([x2, y2])  
	kde = gaussian_kde(xy)(xy)  # Evaluate density at each point

	# Normalize KDE values for colormap
	kde_normalized = (kde - kde.min()) / (kde.max() - kde.min())					

	# Scatter plot with density-based coloring
	sc = ax.scatter(x2, y2, c=kde_normalized, cmap='coolwarm', marker='o', s=8, zorder=3)

	# Colorbar
	cbar = plt.colorbar(sc, ax=ax, label="Density (KDE)")



fig = plt.figure()
ax = fig.add_subplot(111)


#CN-CN
#fig = plt.figure(figsize=(6, 6))  # Adjusting figure size to be square
#ax = fig.add_subplot(111, aspect='equal')  # Setting aspect ratio to equal
#pmd_spaghetti_plt(2001,1,2,ax,data['data'],data['hop'])
#ax.set_xlim(1.2,5.3)
#ax.set_ylim(1.2,5.3)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('DBH-CN-CN.png',dpi=1200)
#ax.set_xlabel(r'C-N (Å)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'C-N (Å)',color='red',fontsize=14,labelpad=8)


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

#D 15 2 3 5 D 14 1 4 5
pmd_spaghetti_plt(2001,4,5,ax,data['data'],data['hop'])
ax.set_xlim(40,185)
ax.set_ylim(40,185)
plt.xticks(fontsize=15, weight='bold')
plt.yticks(fontsize=15, weight='bold')
ax.set_xlabel(r'Dihedral H-C-C-C (°)',color='blue',fontsize=14,labelpad=8)
ax.set_ylabel(r'Dihedral H-C-C-C (°)',color='green',fontsize=14,labelpad=8)
plt.savefig('DBH-Dihedral-Inversion.png',dpi=1400,bbox_inches='tight')


#D 15 2 3 5 D 4 1 2 3
#pmd_spaghetti_plt(2001,4,5,ax,data['data'],data['hop'])
#ax.set_xlim(40,185)
#ax.set_ylim(40,185)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('DBH-Dihedral-Inversion.png',dpi=1400)
#ax.set_xlabel(r'Dihedral 15 2 3 5 (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'Dihedral 4 1 2 3 (°)',color='red',fontsize=14,labelpad=8)
#plt.savefig('DBH-Dihedral-Inversion.png',dpi=1400)

#CC-D 4 1 2 3
#pmd_spaghetti_plt(2001,6,0,ax,data['data'],data['hop'])
#ax.set_xlim(0,180)
#ax.set_ylim(1.2,2.7)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('DBH-CC-Dihedral.png',dpi=1400)
#ax.set_xlabel(r'Dihedral H-C-C-C (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'C-C (Å)',color='red',fontsize=14,labelpad=8)

#CN-Dihedral-4
#pmd_spaghetti_plt(2001,2,5,ax,data['data'],data['hop'])
#ax.set_xlim(0,180)
#ax.set_ylim(1.2,2.7)
#plt.xticks(fontsize=15, weight='bold')
#plt.yticks(fontsize=15, weight='bold')
#plt.savefig('CN-Dihedral-4.png',dpi=1400)
#ax.set_xlabel(r'CN-Dihedral-4 (°)',color='blue',fontsize=14,labelpad=8)
#ax.set_ylabel(r'C-C (Å)',color='red',fontsize=14,labelpad=8)

##parameters: B 3 4 B 3 6 B 4 7 B 6 7 D 15 2 3 5 D 14 1 4 5 D 4 1 2 3
