import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
import pandas as pd
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

"""
ADD DESCRIPTION
"""
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