import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def SEM_csv_to_df(path,skip_empty=True):
    """
    Reads SEM compositional data from a csv via np. 
    It slices the csv containig data from several SEM images into each image.
    Converts the data of each image into df and appends them all together.
    Extracts only the necessary columns
    
    Notes:
    The initial csv cannot contain ae, micrometer or circle (in direction) symbols. 
    You need to convert all the micrometers into um and remove the circles. 
    """
    
    #read data as np array
    data=np.genfromtxt("data/chemical_raw/"+path,dtype=str,delimiter=',')
    
    #get column with indexes of each data section
    index=data[:,0]
    id_index=[int(index) for index,value in enumerate(index) if value=="Id"]
    
    #iterate over the indexes of each data section, convert it to df and append them all together 
    for index,value in enumerate(id_index[:]):
        
        # slice each data section from main data
        if index!=len(id_index)-1:
            subdata=data[id_index[index]:id_index[index+1]-4,:]
        if index==len(id_index)-1:
            subdata=data[id_index[index]:-4,:]
            
        #convert to df and remove duplicates
        df=pd.DataFrame(subdata[1:,:],columns=subdata[0,:])
        df=df.loc[:,~df.columns.duplicated()].reset_index(drop=True)
        
        #concatenate all dfs
        if index==0:
            final_df=df
        if index>0:
            final_df=pd.concat([final_df,df],ignore_index=True)
            
    ### format and filtering ###
    
    #convert Nan and empty spaces into 0
    final_df[pd.isna(final_df)]=0
    final_df[final_df.values=='']=0
            
    #keep only interest columns
    df_SEM=extract_columns(final_df,skip_empty)
    
    #save to csv within outputs folder
    df_SEM.to_csv("data/chemical_processed/"+path[:-4]+"_processed.csv",index=False)  

    return df_SEM
    
def extract_columns(df_SEM,skip_empty=True):
    """
    Given a df with SEM compositional output, it filters the Id, class, ECD, count and weight percentage columns
    """
    
    #get the name of ECD label (as this sometimes varies from export to export)
    try:    
        ECD = [s for s in df_SEM.columns if "ECD" in s][0]
    except IndexError:
        print("Error: .csv does not contain ECD column")
        pass
    
    #basic columns to extract
    basic_labels=['Id', 'Class', ECD,'Count']
    
    #get list of elements present in the df
    element_list=[x for x in df_SEM.columns if ("%" in x) and not ("?" in x)]
    
    #merge basic columns and element composition columns
    extract=basic_labels+element_list
    
    #extract columns
    df_SEM=df_SEM[extract]
    
    #rename ECD column to correct format
    df_SEM=df_SEM.rename(columns={ECD:"ECD (um)"})

    #convert data into numeric
    df_SEM[element_list]=df_SEM[element_list].astype(float)
    df_SEM[['ECD (um)','Count']]=df_SEM[['ECD (um)','Count']].astype(float)
    
    #remove particles with 0 in category
    if skip_empty==True:
        df_SEM=df_SEM.drop(df_SEM[df_SEM['Class']==0].index) 
    
    return df_SEM



def plot_SEM_composition(
       path,
       nbins_submicron,                        
       nbins_supermicron,
       lower_lim,
       upper_lim,
       mag_threshold,
       count_position,
       fontsize,
       hatches,
       hatches_style,
       classification_function,
       legend_plot=True
):
    """
    Plots size-resolved chemical composition of an SEM sample
    
    Inputs:
       - path: path to the csv file with the particles
       - nbins_submicron: number of bins for high magnification analysis                       
       - nbins_supermicron: number of bins for low magnification analysis
       - lower_lim: lower limit of the high magnification analysis
       - upper_lim: upper limit of the analysis
       - mag_threshold: lower limit of the low magnification analysis
       - count_position: position of the number of particles in each bin 9set to 0.5 for large plots)
       - fontsize: fontsize of the plot. set to 12 for large plots
       - hatches: can be True or False if hatches are required
       - hatches_style: can be 1 or 2 for different styles
       - classification_function: set to be a previously imported classification function
           - classification_alaska,
           - classification_barbados,
           - classification_iceland,
           - classification_standard
       - legend_plot: can be True or False if legend is needed or not
       
    Returns:
     - plot
     - finaldata: array with dimensions of number of bins x number of categories.
         Containins the fraction of particles in each size bin in each category.
    """
    
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
        
    ### print number of analysed particles ###
    print("N = "+str(np.sum(totals)))
    
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
    plt.bar(xnames, np.ones(finaldata.shape[1]), width,linewidth=fontsize/6.0,edgecolor='k',fill=False)
    
    #add legend
    if legend_plot==True:
        lgd=plt.legend(
            [p[0] for p in plots[::-1]], 
            topic_labels[::-1],
            bbox_to_anchor=(1.05, 1), 
            loc=2,borderaxespad=0.,
            fontsize=fontsize+2
        )
    
    #plot settings
    plt.ylim((0, 1))
    plt.ylabel('proportion',fontsize=fontsize+2)
    plt.xlabel('diameter ($\mathrm{\mu m}$)',fontsize=fontsize+2)
    plt.yticks(np.linspace(0,1,5),fontsize=fontsize)
    plt.tick_params(axis='x',which='both',bottom=False,top=False) 
    plt.tick_params(axis='y',which='both',bottom=False,top=False) 
    
    #round total number of particles and make it str
    totals=totals.round(0).astype(str)
      
    # add the number of particles per bin on top of the bar
    for index,a in enumerate(totals):
        plt.annotate(a[:-2],((xnames[index]-count_position),0.94),fontsize=fontsize-2)
    
    # define bin label position
    xnames_label=np.linspace(0,len(bins),len(bins))
    
    #plot bin labels
    plt.xticks(xnames_label,bins,fontsize=fontsize)
    
    
    
    return finaldata
