#Script adapted by Leticia A. Gomes to determine geometric parameters from xyz files 
import sys,os
import pandas as pd
import numpy as np
from numpy import linalg as la
import sys

def BND(xyz,var):
    ## This function calculate distance
    ## a<->b

    var=[int(x) for x in var]
    a,b=var[0:2]

    v1=xyz[a-1]
    v2=xyz[b-1]
    r=la.norm(v1-v2)

    return r

def AGL(xyz,var):
    ## This function calculate angle
    ## a<-b->c

    var=[int(x) for x in var]
    a,b,c=var[0:3]

    r1=np.array(xyz[a-1])
    r2=np.array(xyz[b-1])
    r3=np.array(xyz[c-1])
    v1=r1-r2
    v2=r3-r2
    v1=v1/la.norm(v1)
    v2=v2/la.norm(v2)
    cosa=np.dot(v1,v2)
    alpha=np.arccos(cosa)*57.2958

    return alpha

def Deh(xyz,para):
    ## This function calculate dihedral distance
    ##   n1    n2
    ##    |    |
    ## a<-b-><-c->d

    r1=xyz[para[0]-1]
    r2=xyz[para[1]-1]
    r3=xyz[para[2]-1]
    r4=xyz[para[3]-1]
    v1=r1-r2
    v2=r3-r2
    v3=r2-r3
    v4=r4-r3
    n1=np.cross(v1,v2)
    n2=np.cross(v3,v4)
    n1=n1/la.norm(n1)
    n2=n2/la.norm(n2)
    cosb=np.dot(n1,n2)
    value=np.arccos(cosb)*57.2958

    return value

def DHD2(xyz,var):
    ## This function calculate dihedral angle involving dummpy center
    ##   n1    n2
    ##    |    |
    ## a,b<-c-><-d->e,f

    a,b,c,d,e,f=var
    r1=np.array(xyz[a-1])
    r2=np.array(xyz[b-1])
    r3=np.array(xyz[c-1])
    r4=np.array(xyz[d-1])
    r5=np.array(xyz[e-1])
    r6=np.array(xyz[f-1])
    v1=(r1+r2)/2-r3
    v2=r4-r3
    v3=r3-r4
    v4=(r5+r6)/2-r4
    n1=np.cross(v1,v2)
    n2=np.cross(v3,v4)
    n1=n1/la.norm(n1)
    n2=n2/la.norm(n2)
    cosb=np.dot(n1,n2)
    beta=np.arccos(cosb)*57.2958
    axis=np.cross(n1,n2)
    pick=np.argmax(np.abs(axis))
    sign=np.sign(axis[pick]/v2[pick])  # find the projection with largest magnitude (non-zero), then just compare it to avoid 0/0
    if sign == -1:
        beta=360-beta

    return beta

## read xyz file
input=open(sys.argv[1]).read().splitlines()
coord_dict={}
for n, line in enumerate(input):
    if 'traj' in line:
        natom=int(input[n-1])       # number of atom is at previous line
        ntraj=int(line.split()[1])  # trajectory index is the second string in the comment line
        time=int(line.split()[3])   # time step is the forth string in the	comment	line
        coord=np.array([x.split() for x in input[n+1:int(n+1+natom)]])[:,1:4].astype(float)  # extract and format  coordinate into numpy array
        coord_dict[ntraj]=[time,coord] # store data into dictionary for later usage
key=[]
CN1=[]
CN2=[]
N2=[]
CC=[]
D=[]
CC1=[]
CC2=[]
CC3=[]
CC4=[]
CC5=[]
CH1=[]
CH2=[]
CH3=[]
CH4=[]
CH5=[]
CH6=[]
CH7=[]
CH8=[]

# loop over all structures
for ntraj,snapshot in coord_dict.items():
    time,coord=snapshot
    key.append(ntraj)
    CN1.append(BND(coord,[3,6]))
    CN2.append(BND(coord,[4,7]))
    N2.append(BND(coord,[7,6]))
    CC.append(BND(coord,[4,3]))
    D.append(Deh(coord,[15,2,3,5]))
    CC1.append(BND(coord,[1,4]))
    CC2.append(BND(coord,[1,2]))
    CC3.append(BND(coord,[2,3]))
    CC4.append(BND(coord,[3,5]))
    CC5.append(BND(coord,[5,4]))
    CH1.append(BND(coord,[4,10]))
    CH2.append(BND(coord,[1,13]))
    CH3.append(BND(coord,[1,14]))
    CH4.append(BND(coord,[2,12]))
    CH5.append(BND(coord,[2,15]))
    CH6.append(BND(coord,[3,11]))
    CH7.append(BND(coord,[5,8]))
    CH8.append(BND(coord,[5,9]))

# Classify species
classification = []
for i in range(len(key)):
    if CN1[i] < 2.0 and CN2[i] < 2.0:  # Fixed indexing issue
        classification.append("Reactant")
    elif CN1[i] < 2.0 or CN2[i] < 2.0:
        classification.append("1-Bond")
    elif CC1[i] > 2.0 or CC2[i] > 2.0 or CC3[i] > 2.0 or CC4[i] > 2.0 or CC5[i] > 2.0:
        classification.append("Unphysical CC")
    elif CH1[i] > 2.0 or CH2[i] > 2.0 or CH3[i] > 2.0 or CH4[i] >2.0 or CH5[i] > 2.0 or CH6[i] > 2.0 or CH7[i] > 2.0 or CH8[i] > 2.0:
        classification.append("Unphysical CH")
    elif CC[i] > 2.0:
        classification.append("Diradical")
    elif D[i] < 110.0:
        classification.append("Retained")
    else:
        classification.append("Inverted")

results_df=pd.DataFrame(list(zip(key,CN1,CN2,N2,CC,D,classification)),columns =['Trajectory','C3-N6 A','C4-N7 A','N6-N7 A','C3-C4 A','D15-2-3-5','Specie'])   
results_df.to_csv('classify_product.csv',index=False) 

# Count and print the number of each species
from collections import Counter
species_counts = Counter(classification)
for species, count in species_counts.items():
    print(f"{species}: {count}")
