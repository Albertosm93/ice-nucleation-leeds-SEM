# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:15:52 2023

@author: py15asm
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
import pandas as pd
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

###functions###
def classification_barbados(bins,initiate=False,finaldata=False):
    
    #labels
    topic_labels = ['Other','Carbonaceous','S rich','Na rich','Metal rich','Ca rich','Al-Si rich','Si only']
    
    if initiate==True:
        #number of categories
        cath=8
        
        #initiate final_data array
        finaldata=np.zeros((cath,len(bins)-1))
        
    if initiate==False:
        finaldata[0,j]=float(types.count('Cl')+types.count('Cl+K')+types.count('No oxigen')+types.count('All')+types.count('Mainly O')+types.count('No Classification'))/len(types)
                
        finaldata[1,j]=float(types.count('Oxygen only')+types.count('Oxigen+Cu(trace)')+types.count('Oxygen +K')+types.count('P rich')+types.count( 'False Si '))/len(types)   
        
        finaldata[2,j]=float(types.count('Sulfate aerosol'))/len(types)  
    
        finaldata[3,j]=float(types.count('Aged Sea Aerosol')+types.count('Sea aerosol (Cl)')+types.count('Sea aerosol'))/len(types) 
        
        finaldata[4,j]=float(types.count('Metallic Fe')+types.count('Metallic Ti')+types.count('Metallic Cu')+types.count('Metallic Mn')+types.count('Metallic Al')+types.count('False Cr')+types.count('Metallic Zn')+types.count('Metallic Pb'))/len(types)
        
        finaldata[5,j]=float(types.count('Ca rich+NaCl')+types.count('Ca rich')+types.count('Ca pure')+types.count('Gypsium'))/len(types)
        
        finaldata[6,j]=float(types.count('Silica mixtures +NaCl')+types.count('Aluminosillicates+NaCl')+types.count('Sillica mixtures ')+types.count('Aluminosillicates'))/len(types)
        
        finaldata[7,j]=float(types.count('Sillica'))/len(types)
        
    return finaldata, topic_labels

def classification_iceland(bins,initiate=False,finaldata=False):
    
    #labels
    topic_labels = ['Other','Carbonaceous','S rich','Cl rich','Na rich','Metal rich','Ca rich','Al-Si rich','Si only','Si rich']
    
    if initiate==True:
        #number of categories
        cath=10
        
        #initiate final_data array
        finaldata=np.zeros((cath,len(bins)-1))
        
        return finaldata, topic_labels
        
    if initiate==False:
        finaldata[0,j]=float(types.count('No oxigen')+types.count('All')+types.count('Mainly O')+types.count('No Classification'))/len(types)
                
        finaldata[1,j]=float(types.count('Oxygen only')+types.count('Oxigen+Cu(trace)')+types.count('Oxygen +K')+types.count('P rich')+types.count( 'False Si '))/len(types)   
        
        finaldata[2,j]=float(types.count('Sulfate aerosol'))/len(types) 
        
        finaldata[3,j]=float(types.count('Cl')+types.count('Cl+K'))/len(types) 

        finaldata[4,j]=float(types.count('Aged Sea Aerosol')+types.count('Sea aerosol (Cl)')+types.count('Sea aerosol'))/len(types) 
        
        finaldata[5,j]=float(types.count('Metallic Fe')+types.count('Metallic Ti')+types.count('Metallic Cu')+types.count('Metallic Mg')+types.count('Metallic Mn')+types.count('Metallic Al')+types.count('False Cr')+types.count('Metallic Zn')+types.count('Metallic Pb'))/len(types) 
        
        finaldata[6,j]=float(types.count('Ca rich')+types.count('Ca rich+NaCl')+types.count('Ca pure')+types.count('Gypsium'))/len(types)
        
        finaldata[7,j]=float(types.count('Aluminosillicates')+types.count('Aluminosillicates+NaCl'))/len(types) 
        
        finaldata[8,j]=float(types.count('Sillica'))/len(types)
        
        finaldata[9,j]=float(types.count('Sillica mixtures ')+types.count('Silica mixtures +NaCl'))/len(types)
        
        
        return finaldata


###inputs###

nbins_submicron=3 
nbins_supermicron=7
lower_lim=0.3 #in um
upper_lim=20 #in um

count_position=0.5

fontsize=12

hatches=True

hatches_style=1 # could be 1 or 2

colors=['green','lightgrey','yellow','cornflowerblue','firebrick','sandybrown','orange','khaki']

classification_function=classification_iceland

path=r'C:\Users\py15asm\Documents\GitHub\ice-nucleation-leeds-SEM\data\chemical_processed\examples\\171002_1624_1640_432L_OP_UP_chemical_processed.csv'
plt.figure(figsize=(6,3),dpi=300)


#colors iceland
colors=[]
K=10
for k in range(K):
    color = plt.cm.Paired((K-k)/float(K), 1)
    colors.append(color)
        





        
    
    
    
#hatches
if hatches==True:
    plt.rcParams['hatch.linewidth'] = 0.3
if hatches==False:
    plt.rcParams['hatch.linewidth'] = 0
    
patterns = [ "\\" , "|" , "/","\\" , "|" , "/","\\" , "|" , "/","\\" , "|" , "/"  ]
if hatches_style==2:
    patterns = [ "+" , "x" , "o","+" , "x" , "o","+" , "x" , "o","+" , "x" , "o"]

#bins
bins=np.zeros(nbins_submicron+nbins_supermicron+1)
bins[:(nbins_submicron+1)]=np.logspace(np.log10(lower_lim),np.log10(1),nbins_submicron+1) 
bins[(nbins_submicron+1):]=np.logspace(np.log10(1),np.log10(upper_lim),nbins_supermicron+1)[1:]

bins=np.around(bins,1)


#pandas data read
data=pd.read_csv(path)
ECDm=data["ECD (um)"].values
Class=data["Class"].values

totals=np.zeros(len(bins)-1)

finaldata,topic_labels=classification_function(bins,initiate=True)
for j in range(0,len(bins)-1):
    types=[]
    for i in range(0,len(ECDm)):
        if float(ECDm[i])>=bins[j] and float(ECDm[i])<bins[j+1] and not (Class[i]=='Metallic Cr'):
            totals[j]+=1
            types.append(Class[i])
        
    if not types==[]:   
        finaldata=classification_function(bins,initiate=False,finaldata=finaldata)
                

        
        
        
                
xnames=np.linspace(0,len(bins)-1,len(bins)-1)+0.5
xticks=np.zeros(len(bins)-1)
  

K, N = finaldata.shape
ind = np.arange(N)
width = 1.08

plt.bar(xnames, finaldata[0,:])

plots = []
height_cumulative = np.zeros(N)


for k in range(K):
    color = plt.cm.Paired((K-k)/float(K), 1)
    if k == 0:
        p = plt.bar(xnames, finaldata[k,:], width, color=colors[k],linewidth=0.4,hatch=patterns[k])

    else:
        p = plt.bar(xnames, finaldata[k,:], width, bottom=height_cumulative,color=colors[k],linewidth=0.3,edgecolor='k',hatch=patterns[k])
    height_cumulative += finaldata[k,:]
    plots.append(p)   
    
plt.bar(xnames, np.ones(finaldata.shape[1]), width,linewidth=2,edgecolor='k',fill=False)

plt.ylim((0, 1))
plt.ylabel('proportion',fontsize=fontsize+4)
plt.xlabel('diameter ($\mathrm{\mu m}$)',fontsize=fontsize+4)


total=(bins)

    
finaldatabulk=finaldata*totals  
errfinaldatabulk=np.sqrt(finaldatabulk)
errfinaldata=errfinaldatabulk/totals


for i in range(0,len(totals)):
    totals[i]=str(totals[i])
  
for i in range(0,N):
    a=str(totals[i])
    plt.annotate(a[:-2],((xnames[i]-count_position),0.95),fontsize=fontsize-2)
plt.annotate('1',(15,0.05))

    

    
    
xnames2=np.linspace(0,len(bins),len(bins))

plt.xticks(xnames2,bins,fontsize=fontsize)

print(np.sum(finaldata,axis=0)) 
print('this should be all 1s')
plt.yticks(np.linspace(0,1,5),fontsize=fontsize)

plt.tick_params(axis='x',which='both',bottom=False,top=False) 
plt.tick_params(axis='y',which='both',bottom=False,top=False) 
#p1.tick_params(labelsize=13)

lgd=plt.legend([p[0] for p in plots[::-1]], topic_labels[::-1],bbox_to_anchor=(1.05, 1), loc=2,borderaxespad=0.,fontsize=fontsize+2)
#plt.text(13.5, 0.02, 'N= '+str(np.sum(totals))[:-2], fontsize=16)



