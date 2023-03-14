import pandas as pd
import numpy as np



def SEM_csv_to_df(path):
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
    data=np.genfromtxt("data/"+path,dtype=str,delimiter=',')
    
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
    df_SEM=extract_columns(final_df)
    
    #save to csv within outputs folder
    df_SEM.to_csv("outputs/"+path[:-4]+"_processed.csv",index=False)  

    return df_SEM
    
def extract_columns(df_SEM):
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
    df_SEM=df_SEM.drop(df_SEM[df_SEM['Class']==0].index)
    
    

    
    return df_SEM
