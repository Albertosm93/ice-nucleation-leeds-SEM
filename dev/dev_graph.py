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
def classification_barbados(bins,j=False,classes=False,initiate=False,finaldata=False):
    
    #labels
    topic_labels = ['Other','Carbonaceous','S rich','Na rich','Metal rich','Ca rich','Al-Si rich','Si only']
    
    #colors
    colors = ['green','lightgrey','yellow','cornflowerblue','firebrick','sandybrown','orange','khaki']
    
    #all categories in algorithm below
    all_categories = np.array(['No oxigen','No oxigen','All','Mainly O','Oxygen only','Oxigen+Cu(trace)','Oxygen +K','P rich' 'No Classification','False Si ','Sulfate aerosol','Metallic Fe','Metallic Ti','Metallic Cu','Metallic Mn','Metallic Al','Metallic Cr','False Cr','Metallic Zn','Metallic Pb','Aged Sea Aerosol','Sea aerosol (Cl)','Sea aerosol','Cl','Cl+K','Ca rich','Ca rich+NaCl','Ca pure','Gypsium','Aluminosillicates','Aluminosillicates+NaCl','Sillica','Sillica mixtures ','Silica mixtures +NaCl','No Classification'])

    if initiate==True:
        #number of categories
        cath=len(topic_labels)
        
        #initiate final_data array
        finaldata=np.zeros((cath,len(bins)-1))
        
        return finaldata, topic_labels, colors, all_categories
    
        
    if initiate==False:
        finaldata[0,j]=float(classes.count('Cl')+classes.count('Cl+K')+classes.count('No oxigen')+classes.count('All')+classes.count('Mainly O')+classes.count('No Classification'))/len(classes)
                
        finaldata[1,j]=float(classes.count('Oxygen only')+classes.count('Oxigen+Cu(trace)')+classes.count('Oxygen +K')+classes.count('P rich')+classes.count( 'False Si '))/len(classes)   
        
        finaldata[2,j]=float(classes.count('Sulfate aerosol'))/len(classes)  
    
        finaldata[3,j]=float(classes.count('Aged Sea Aerosol')+classes.count('Sea aerosol (Cl)')+classes.count('Sea aerosol'))/len(classes) 
        
        finaldata[4,j]=float(classes.count('Metallic Fe')+classes.count('Metallic Ti')+classes.count('Metallic Cu')+classes.count('Metallic Mn')+classes.count('Metallic Al')+classes.count('False Cr')+classes.count('Metallic Zn')+classes.count('Metallic Pb')+classes.count('Metallic Cr'))/len(classes)
        
        finaldata[5,j]=float(classes.count('Ca rich+NaCl')+classes.count('Ca rich')+classes.count('Ca pure')+classes.count('Gypsium'))/len(classes)
        
        finaldata[6,j]=float(classes.count('Silica mixtures +NaCl')+classes.count('Aluminosillicates+NaCl')+classes.count('Sillica mixtures ')+classes.count('Aluminosillicates'))/len(classes)
        
        finaldata[7,j]=float(classes.count('Sillica'))/len(classes)
        
        return finaldata

def classification_iceland(bins,j=False,classes=False,initiate=False,finaldata=False):
    
    #labels
    topic_labels = ['Other','Carbonaceous','S rich','Cl rich','Na rich','Metal rich','Ca rich','Al-Si rich','Si only','Si rich']
    
    #all categories in algorithm below
    all_categories=np.array(['No oxigen','No oxigen','All','Mainly O','Oxygen only','Oxigen+Cu(trace)','Oxygen +K','P rich' 'No Classification','False Si ','Sulfate aerosol','Metallic Fe','Metallic Ti','Metallic Cu','Metallic Mn','Metallic Al','Metallic Cr','False Cr','Metallic Zn','Metallic Pb','Aged Sea Aerosol','Sea aerosol (Cl)','Sea aerosol','Cl','Cl+K','Ca rich','Ca rich+NaCl','Ca pure','Gypsium','Aluminosillicates','Aluminosillicates+NaCl','Sillica','Sillica mixtures ','Silica mixtures +NaCl','No Classification'])

    #colors
    colors=[]
    K=10
    for k in range(K):
        color = plt.cm.Paired((K-k)/float(K), 1)
        colors.append(color)
        
    if initiate==True:
        #number of categories
        cath=len(topic_labels)
        
        #initiate final_data array
        finaldata=np.zeros((cath,len(bins)-1))
        
        return finaldata, topic_labels, colors, all_categories
        
    if initiate==False:
        finaldata[0,j]=float(classes.count('No oxigen')+classes.count('All')+classes.count('Mainly O')+classes.count('No Classification'))/len(classes)
                
        finaldata[1,j]=float(classes.count('Oxygen only')+classes.count('Oxigen+Cu(trace)')+classes.count('Oxygen +K')+classes.count('P rich')+classes.count( 'False Si '))/len(classes)   
        
        finaldata[2,j]=float(classes.count('Sulfate aerosol'))/len(classes) 
        
        finaldata[3,j]=float(classes.count('Cl')+classes.count('Cl+K'))/len(classes) 

        finaldata[4,j]=float(classes.count('Aged Sea Aerosol')+classes.count('Sea aerosol (Cl)')+classes.count('Sea aerosol'))/len(classes) 
        
        finaldata[5,j]=float(classes.count('Metallic Fe')+classes.count('Metallic Ti')+classes.count('Metallic Cu')+classes.count('Metallic Mg')+classes.count('Metallic Mn')+classes.count('Metallic Al')+classes.count('False Cr')+classes.count('Metallic Zn')+classes.count('Metallic Pb')+classes.count('Metallic Cr'))/len(classes) 
        
        finaldata[6,j]=float(classes.count('Ca rich')+classes.count('Ca rich+NaCl')+classes.count('Ca pure')+classes.count('Gypsium'))/len(classes)
        
        finaldata[7,j]=float(classes.count('Aluminosillicates')+classes.count('Aluminosillicates+NaCl'))/len(classes) 
        
        finaldata[8,j]=float(classes.count('Sillica'))/len(classes)
        
        finaldata[9,j]=float(classes.count('Sillica mixtures ')+classes.count('Silica mixtures +NaCl'))/len(classes)
        
        
        return finaldata

def classification_alaska(bins,j=False,classes=False,initiate=False,finaldata=False):
    
    #labels
    topic_labels = ['Other','Carbonaceous','S rich','Cl rich','Na rich','Metal rich','Ca rich','Al-Si rich','Si only','Si rich']
    
    #all categories in algorithm below
    all_categories=np.array(['No oxigen','No oxigen','All','Mainly O','Oxygen only','Oxigen+Cu(trace)','Oxygen +K','P rich' 'No Classification','False Si ','Sulfate aerosol','Metallic Fe','Metallic Ti','Metallic Cu','Metallic Mn','Metallic Al','Metallic Cr','False Cr','Metallic Zn','Metallic Pb','Aged Sea Aerosol','Sea aerosol (Cl)','Sea aerosol','Cl','Cl+K','Ca rich','Ca rich+NaCl','Ca pure','Gypsium','Aluminosillicates','Aluminosillicates+NaCl','Sillica','Sillica mixtures ','Silica mixtures +NaCl','No Classification'])

    #colors
    colors=['green','lightgrey','yellow','palevioletred','cornflowerblue','firebrick','sandybrown','orange','tan','khaki']
        
    if initiate==True:
        #number of categories
        cath=len(topic_labels)
        
        #initiate final_data array
        finaldata=np.zeros((cath,len(bins)-1))
        
        return finaldata, topic_labels, colors, all_categories
        
    if initiate==False:
        finaldata[0,j]=float(classes.count('No oxigen')+classes.count('All')+classes.count('Mainly O')+classes.count('No Classification'))/len(classes)
                
        finaldata[1,j]=float(classes.count('Oxygen only')+classes.count('Oxigen+Cu(trace)')+classes.count('Oxygen +K')+classes.count('P rich')+classes.count( 'False Si '))/len(classes)   
        
        finaldata[2,j]=float(classes.count('Sulfate aerosol'))/len(classes) 
        
        finaldata[3,j]=float(classes.count('Cl')+classes.count('Cl+K'))/len(classes) 

        finaldata[4,j]=float(classes.count('Aged Sea Aerosol')+classes.count('Sea aerosol (Cl)')+classes.count('Sea aerosol'))/len(classes) 
        
        finaldata[5,j]=float(classes.count('Metallic Fe')+classes.count('Metallic Ti')+classes.count('Metallic Cu')+classes.count('Metallic Mn')+classes.count('Metallic Mg')+classes.count('Metallic Al')+classes.count('False Cr')+classes.count('Metallic Zn')+classes.count('Metallic Pb'))/len(classes)   #+classes.count('Metallic Cr')
        
        finaldata[6,j]=float(classes.count('Ca rich')+classes.count('Ca rich+NaCl')+classes.count('Ca pure')+classes.count('Gypsium'))/len(classes)
        
        finaldata[7,j]=float(classes.count('Aluminosillicates')+classes.count('Aluminosillicates+NaCl'))/len(classes) 
        
        finaldata[8,j]=float(classes.count('Sillica'))/len(classes)
        
        finaldata[9,j]=float(classes.count('Sillica mixtures ')+classes.count('Silica mixtures +NaCl'))/len(classes)
        
        
        return finaldata


def classification_standard(bins,j=False,classes=False,initiate=False,finaldata=False):
    
    #labels
    topic_labels = ['Other','Carbonaceous','S rich','Cl rich','Na rich','Metal rich','Ca rich','Al-Si rich','Si only']
    
    #all categories in algorithm below
    all_categories=np.array(['No oxigen','No oxigen','All','Mainly O','Oxygen only','Oxigen+Cu(trace)','Oxygen +K','P rich' 'No Classification','False Si ','Sulfate aerosol','Metallic Fe','Metallic Ti','Metallic Cu','Metallic Mn','Metallic Al','Metallic Cr','False Cr','Metallic Zn','Metallic Pb','Aged Sea Aerosol','Sea aerosol (Cl)','Sea aerosol','Cl','Cl+K','Ca rich','Ca rich+NaCl','Ca pure','Gypsium','Aluminosillicates','Aluminosillicates+NaCl','Sillica','Sillica mixtures ','Silica mixtures +NaCl','No Classification'])

    #colors
    colors=['green','lightgrey','yellow','palevioletred','cornflowerblue','firebrick','sandybrown','orange','khaki']
        
    if initiate==True:
        #number of categories
        cath=len(topic_labels)
        
        #initiate final_data array
        finaldata=np.zeros((cath,len(bins)-1))
        
        return finaldata, topic_labels, colors, all_categories
        
    if initiate==False:
        finaldata[0,j]=float(classes.count('No oxigen')+classes.count('All')+classes.count('Mainly O')+classes.count('No Classification'))/len(classes)
                
        finaldata[1,j]=float(classes.count('Oxygen only')+classes.count('Oxigen+Cu(trace)')+classes.count('Oxygen +K')+classes.count('P rich')+classes.count( 'False Si '))/len(classes)   
        
        finaldata[2,j]=float(classes.count('Sulfate aerosol'))/len(classes) 
        
        finaldata[3,j]=float(classes.count('Cl')+classes.count('Cl+K'))/len(classes) 

        finaldata[4,j]=float(classes.count('Aged Sea Aerosol')+classes.count('Sea aerosol (Cl)')+classes.count('Sea aerosol'))/len(classes) 
        
        finaldata[5,j]=float(classes.count('Metallic Fe')+classes.count('Metallic Ti')+classes.count('Metallic Cu')+classes.count('Metallic Mn')+classes.count('Metallic Mg')+classes.count('Metallic Al')+classes.count('False Cr')+classes.count('Metallic Zn')+classes.count('Metallic Pb'))/len(classes)   #+classes.count('Metallic Cr')
        
        finaldata[6,j]=float(classes.count('Ca rich')+classes.count('Ca rich+NaCl')+classes.count('Ca pure')+classes.count('Gypsium'))/len(classes)
        
        finaldata[7,j]=float(classes.count('Aluminosillicates')+classes.count('Aluminosillicates+NaCl')+classes.count('Sillica mixtures ')+classes.count('Silica mixtures +NaCl'))/len(classes) 
        
        finaldata[8,j]=float(classes.count('Sillica'))/len(classes)
        
        
        return finaldata

        




    
    
def plot_SEM_composition(nbins_submicron,nbins_supermicron,lower_lim,upper_lim,mag_threshold,count_position,fontsize,hatches,hatches_style,classification_function,path):
    
    ### Hatches ###
    if hatches==True:
        plt.rcParams['hatch.linewidth'] = 0.3
    if hatches==False:
        plt.rcParams['hatch.linewidth'] = 0
        
    patterns = [ "\\" , "|" , "/","\\" , "|" , "/","\\" , "|" , "/","\\" , "|" , "/"  ]
    
    if hatches_style==2:
        patterns = [ "+" , "x" , "o","+" , "x" , "o","+" , "x" , "o","+" , "x" , "o"]
    
    
    ### Bins ###
    bins=np.zeros(nbins_submicron+nbins_supermicron+1) #define bins
    bins[:(nbins_submicron+1)]=np.logspace(np.log10(lower_lim),np.log10(mag_threshold),nbins_submicron+1) #high magniication
    bins[(nbins_submicron+1):]=np.logspace(np.log10(mag_threshold),np.log10(upper_lim),nbins_supermicron+1)[1:] #low magnification
    
    #round bins
    bins=np.around(bins,1)
    
    
    ### Read and filter data ###
    
    data=pd.read_csv(path)
    data=data[data['Count']>0] #filter no counts
    
    #extract ECD and class
    ECDm=data["ECD (um)"].values
    Class=data["Class"].values
    Count=data["Count"].values
    
    
    ### Count classes ###
    
    #total number of included particles per size bin
    totals=np.zeros(len(bins)-1)
    
    #initiate the classification function to get class and ECD array
    finaldata,topic_labels, colors, all_categories = classification_function(bins,initiate=True)
    
    #loop for each size bin and append 
    for j in range(0,len(bins)-1):
        
        # list of particle classes in the size bin
        classes=[]
        
        # fill the list of particle classes in this size bin
        for i in range(0,len(ECDm)):
            if float(ECDm[i])>=bins[j] and float(ECDm[i])<bins[j+1] and not (Class[i]=='Metallic Cr'):
                totals[j]+=1
                classes.append(Class[i])
            
        # if the list of particles is not empty, call classification function to put these particles in the final classified particles array
        if not classes==[]:   
            finaldata=classification_function(bins,j,classes,initiate=False,finaldata=finaldata)
          
            
    K, N = finaldata.shape #K number of categories, N number of size bins
    
    
    ### Testing ###
    
    #calculating total proportion of particles in each size bin
    final_proportion=np.sum(finaldata,axis=0).round(5)
    
    #number of size bins where the proportion is either 1 or 0 (no particles have been missing in the classification)
    full=np.array([np.isin(final_proportion,[0,1])==True]).sum()
    
    #if full is not the same as the number of size bins, some particles have not been classified 
    if full!=N:
        print("Some particles have not been classified: ")
        
        for index, item in enumerate(Class[~np.isin(Class,all_categories)]):
            print("Class: "+item)
            print("ECD: "+str(ECDm[~np.isin(Class,all_categories)][index]))
        
    
    
    ### Plotting ###       
        
    #basic settings and initiating counters
    xnames=np.linspace(0,len(bins)-1,len(bins)-1)+0.5  #bin positions
    K, N = finaldata.shape #K number of categories, N number of size bins
    width = 1.08 #bar width
    plots = [] #list with the plots
    height_cumulative = np.zeros(N) #cumulative height of the bars initiated at 0
    
    # iterate to plot each the bar categories
    for k in range(K):
        #first bar
        if k == 0:
            p = plt.bar(xnames, finaldata[k,:], width, color=colors[k],linewidth=0.4,hatch=patterns[k])
        
        #rest of the bars
        else:
            p = plt.bar(xnames, finaldata[k,:], width, bottom=height_cumulative,color=colors[k],linewidth=0.3,edgecolor='k',hatch=patterns[k])
        
        #add last height to cumulative heights
        height_cumulative += finaldata[k,:]
        
        #append plot
        plots.append(p)   
        
    #plot bar edges (empty bars with thick edges)
    plt.bar(xnames, np.ones(finaldata.shape[1]), width,linewidth=2,edgecolor='k',fill=False)
    
    #add legend
    lgd=plt.legend([p[0] for p in plots[::-1]], topic_labels[::-1],bbox_to_anchor=(1.05, 1), loc=2,borderaxespad=0.,fontsize=fontsize+2)
    
    #plot settings
    plt.ylim((0, 1))
    plt.ylabel('proportion',fontsize=fontsize+4)
    plt.xlabel('diameter ($\mathrm{\mu m}$)',fontsize=fontsize+4)
    plt.yticks(np.linspace(0,1,5),fontsize=fontsize)
    plt.tick_params(axis='x',which='both',bottom=False,top=False) 
    plt.tick_params(axis='y',which='both',bottom=False,top=False) 
    
    #round total number of particles and make it str
    totals=totals.round(0).astype(str)
      
    # add the number of particles per bin on top of the bar
    for index,a in enumerate(totals):
        plt.annotate(a[:-2],((xnames[index]-count_position),0.95),fontsize=fontsize-2)
    
    # define bin label position
    xnames_label=np.linspace(0,len(bins),len(bins))
    
    #plot bin labels
    plt.xticks(xnames_label,bins,fontsize=fontsize)
    
    
    return finaldata
    
    
    
    
    
###inputs###
nbins_submicron=3 
nbins_supermicron=7
lower_lim=0.3 #in um
upper_lim=10 #in um
mag_threshold=0.8
count_position=0.5
fontsize=12
hatches=True
hatches_style=1 # could be 1 or 2

classification_function=classification_standard

path=r'C:\Users\py15asm\Documents\GitHub\ice-nucleation-leeds-SEM\data\chemical_processed\examples\\171002_1624_1640_432L_OP_UP_chemical_processed.csv'
#path=r'C:\Users\py15asm\Documents\GitHub\ice-nucleation-leeds-SEM\data\\chemical_processed\180316_FAAM_2045_2126_0.51V_chem_processed.csv'

plt.figure(figsize=(6,3),dpi=300)

plot_SEM_composition(nbins_submicron,nbins_supermicron,lower_lim,upper_lim,mag_threshold,count_position,fontsize,hatches,hatches_style,classification_function,path)  
    
