# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:13:46 2023

@author: py15asm
"""



import os, os.path
#os.environ["PROJ_LIB"] = "C:\Users\py15asm\AppData\Local\Continuum\\anaconda2\\Library\share" #importat to add this before importing Basemap
from scipy.io.idl import readsav
import glob
import numpy as np
import pandas as pd
from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
from netCDF4 import Dataset
from scipy.optimize import curve_fit
import matplotlib
matplotlib.style.use('classic')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

#%%

counter_tot_final=0

#if ymax=0, it doesn;t limit the y axis. If it is above 0, that will be the limit of the y axis
def SEM_analysys(path__lowmag,path__highmag,lowlim,medlim,uplim,air,area__highmag,area__lowmag,path,pathcloud,pathcsv,pathpcasp,pathcdp,path_chem,title,ymax,letter,suptitle,legend='yes'):
    counter_tot=0
    filt=11*10**8 #11cm2 to um2
    
    ## bins    
    #morph
    lower=np.log10(lowlim)
    middle=np.log10(medlim)
    upper=np.log10(uplim)
    nbins_highmag=5
    nbins_lowmag=10
    bins_highmag=np.logspace(lower,middle,nbins_highmag)
    bins_lowmag=np.logspace(middle,upper,nbins_lowmag)
    
    #composition
    BinsC=np.logspace(np.log10(0.3),np.log10(20),11)        
    BinsC[:4]=np.logspace(np.log10(lowlim),np.log10(1),4) 
    BinsC[3:]=np.logspace(np.log10(1),np.log10(10),8)
    BinsC=BinsC.round(1)
    
    #make morph bins from composition bins
    bins_highmag=BinsC[:4]
    bins_lowmag=BinsC[3:]
    nbins_highmag=len(bins_highmag)
    nbins_lowmag=len(bins_lowmag)
    
    #importing _lowmag means low mag and _highmag high mag
    data__lowmag=np.genfromtxt(path__lowmag,delimiter=',',skip_header=0,dtype=np.float64)
    data__highmag=np.genfromtxt(path__highmag,delimiter=',',skip_header=0,dtype=np.float64)
    
    data__lowmag_s=np.genfromtxt(path__lowmag,delimiter=',',skip_header=0,dtype=str)
    data__highmag_s=np.genfromtxt(path__highmag,delimiter=',',skip_header=0,dtype=str)
    
    #arranging data
    ECD__lowmag=[]
    ECD__highmag=[]
    
    
    counter=0
    for i in range(0,len(data__lowmag)): 
        if not math.isnan(data__lowmag[i,0]) and not data__lowmag_s[i,1]=='Metallic Cr':
            ECD__lowmag.append(data__lowmag[i,6])
            counter=counter+1
            counter_tot=counter_tot+1
            
    counter=0
    for i in range(0,len(data__highmag)):
        if not math.isnan(data__highmag[i,0]) and not data__highmag_s[i,1]=='Metallic Cr':
            ECD__highmag.append(data__highmag[i,6])
            ECD__lowmag.append(data__highmag[i,6])
            counter=counter+1
            counter_tot=counter_tot+1
    

    
    bincentre_highmag=np.zeros(len(bins_highmag)-1)
    bincentre_lowmag=np.zeros(len(bins_lowmag)-1)
    dlogDp_highmag=np.zeros(len(bins_highmag)-1)
    dlogDp_lowmag=np.zeros(len(bins_lowmag)-1)
    
    for i in range(0,len(bins_highmag)-1):
        bincentre_highmag[i]=(bins_highmag[i]+bins_highmag[i+1])/2
        dlogDp_highmag[i]=math.log10(bins_highmag[i+1])-math.log10(bins_highmag[i])
        
    for i in range(0,len(bins_lowmag)-1):
        bincentre_lowmag[i]=(bins_lowmag[i]+bins_lowmag[i+1])/2
        dlogDp_lowmag[i]=math.log10(bins_lowmag[i+1])-math.log10(bins_lowmag[i])
    
    hist__highmag,edges=np.histogram(ECD__highmag,bins=bins_highmag)
    hist__lowmag,edges=np.histogram(ECD__lowmag,bins=bins_lowmag)
    
    errhist__highmag=np.sqrt(hist__highmag)
    errhist__lowmag=np.sqrt(hist__lowmag)
    
    #calculations Nstands for dndlogdp and A dadlogdp
    Nhist__highmag=hist__highmag*filt/(area__highmag*air*dlogDp_highmag)
    Nhist__lowmag=hist__lowmag*filt/(area__lowmag*air*dlogDp_lowmag)
    
    errNhist__highmag=errhist__highmag*filt/(area__highmag*air*dlogDp_highmag)
    errNhist__lowmag=errhist__lowmag*filt/(area__lowmag*air*dlogDp_lowmag)
    
    Ahist__highmag=Nhist__highmag*math.pi*bincentre_highmag**2
    Ahist__lowmag=Nhist__lowmag*math.pi*bincentre_lowmag**2
    
    errAhist__highmag=errNhist__highmag*math.pi*bincentre_highmag**2
    errAhist__lowmag=errNhist__lowmag*math.pi*bincentre_lowmag**2
    
    N=np.concatenate((Nhist__highmag,Nhist__lowmag))
    
    
    errN=np.concatenate((errNhist__highmag,errNhist__lowmag))
    
    
    A=np.concatenate((Ahist__highmag,Ahist__lowmag))
    
    errA=np.concatenate((errAhist__highmag,errAhist__lowmag))
    
    
    Bins=np.concatenate((bincentre_highmag,bincentre_lowmag))
    Binsupedge=np.concatenate((bins_highmag[:-1],bins_lowmag[:]))
    dlogDp=np.concatenate((dlogDp_highmag,dlogDp_lowmag))
    
    
    #getting CEDA data
    if pcaspflag=='yes':
        
        data=Dataset(path)
        datacloud=Dataset(pathcloud)
        datacsv=np.genfromtxt(pathcsv,delimiter=',',dtype=str)
        
        
        
        nstops=(len(datacsv[:,0])-3)/2
        
        sh0=0
        sm0=0
        ss0=0
        eh0=0
        em0=0
        es0=0
        
        sh1=0
        sm1=0
        ss1=0
        eh1=0
        em1=0
        es1=0
        
        sh2=0
        sm2=0
        ss2=0
        eh2=0
        em2=0
        es2=0
        
        sh3=0
        sm3=0
        ss3=0
        eh3=0
        em3=0
        es3=0
        
        if nstops==0:
            sh0=int(str(datacsv[1,0])[:2])
            sm0=int(str(datacsv[1,0])[2:4])
            ss0=int(str(datacsv[1,0])[4:6])
            eh0=int(str(datacsv[2,0])[:2])
            em0=int(str(datacsv[2,0])[2:4])
            es0=int(str(datacsv[2,0])[4:6])
            
        if nstops==1:
            sh0=int(str(datacsv[1,0])[:2])
            sm0=int(str(datacsv[1,0])[2:4])
            ss0=int(str(datacsv[1,0])[4:6])
            eh0=int(str(datacsv[2,0])[:2])
            em0=int(str(datacsv[2,0])[2:4])
            es0=int(str(datacsv[2,0])[4:6])
            sh1=int(str(datacsv[3,0])[:2])
            sm1=int(str(datacsv[3,0])[2:4])
            ss1=int(str(datacsv[3,0])[4:6])
            eh1=int(str(datacsv[4,0])[:2])
            em1=int(str(datacsv[4,0])[2:4])
            es1=int(str(datacsv[4,0])[4:6])
            
        if nstops==2:
            sh0=int(str(datacsv[1,0])[:2])
            sm0=int(str(datacsv[1,0])[2:4])
            ss0=int(str(datacsv[1,0])[4:6])
            eh0=int(str(datacsv[2,0])[:2])
            em0=int(str(datacsv[2,0])[2:4])
            es0=int(str(datacsv[2,0])[4:6])
            sh1=int(str(datacsv[3,0])[:2])
            sm1=int(str(datacsv[3,0])[2:4])
            ss1=int(str(datacsv[3,0])[4:6])
            eh1=int(str(datacsv[4,0])[:2])
            em1=int(str(datacsv[4,0])[2:4])
            es1=int(str(datacsv[4,0])[4:6])
            sh2=int(str(datacsv[5,0])[:2])
            sm2=int(str(datacsv[5,0])[2:4])
            ss2=int(str(datacsv[5,0])[4:6])
            eh2=int(str(datacsv[6,0])[:2])
            em2=int(str(datacsv[6,0])[2:4])
            es2=int(str(datacsv[6,0])[4:6])
            
        if nstops==3:
            sh0=int(str(datacsv[1,0])[:2])
            sm0=int(str(datacsv[1,0])[2:4])
            ss0=int(str(datacsv[1,0])[4:6])
            eh0=int(str(datacsv[2,0])[:2])
            em0=int(str(datacsv[2,0])[2:4])
            es0=int(str(datacsv[2,0])[4:6])
            sh1=int(str(datacsv[3,0])[:2])
            sm1=int(str(datacsv[3,0])[2:4])
            ss1=int(str(datacsv[3,0])[4:6])
            eh1=int(str(datacsv[4,0])[:2])
            em1=int(str(datacsv[4,0])[2:4])
            es1=int(str(datacsv[4,0])[4:6])
            sh2=int(str(datacsv[5,0])[:2])
            sm2=int(str(datacsv[5,0])[2:4])
            ss2=int(str(datacsv[5,0])[4:6])
            eh2=int(str(datacsv[6,0])[:2])
            em2=int(str(datacsv[6,0])[2:4])
            es2=int(str(datacsv[6,0])[4:6])
            sh3=int(str(datacsv[7,0])[:2])
            sm3=int(str(datacsv[7,0])[2:4])
            ss3=int(str(datacsv[7,0])[4:6])
            eh3=int(str(datacsv[8,0])[:2])
            em3=int(str(datacsv[8,0])[2:4])
            es3=int(str(datacsv[8,0])[4:6])
        
        datapcasp=np.genfromtxt(pathpcasp,delimiter=',',skip_header=3)
        datacdp=np.genfromtxt(pathcdp,delimiter=',',skip_header=3)
        timecloud=np.array(datacloud.variables['Time'][:])
        PCASPtot=np.array(datacloud.variables['PCAS2CON'][:])
        CDPtot=np.array(datacloud.variables['CDP_CONC'][:])
        itime=int(timecloud[0])
        start0=sh0*3600+sm0*60+ss0
        start1=sh1*3600+sm1*60+ss1
        start2=sh2*3600+sm2*60+ss2
        start3=sh3*3600+sm3*60+ss3
        end0=eh0*3600+em0*60+es0
        end1=eh1*3600+em1*60+es1
        end2=eh2*3600+em2*60+es2
        end3=eh3*3600+em3*60+es3
        
        #obtaining the bincentres and the logarithmic bin widths from the pcasp cdp
        Bincenterpcasp=datapcasp[6,:][1:-1]
        errBincenterpcasp=datapcasp[7,:][1:-1]
        Bincentercdp=datacdp[6,:][1:-1]
        errBincentercdp=datacdp[7,:][1:-1]
        dlogDppcasp=datapcasp[12,:][1:-1]
        errdlogDppcasp=datapcasp[13,:][1:-1]
        dlogDpcdp=datacdp[12,:][1:-1]
        errdlogDpcdp=datacdp[13,:][1:-1]
        
        #calculating flows to calculate total number of counts, 
        #the pcasp flow rate is in cm3 s-1, the air speed in ms-1 and the area in mm2.We convert these two last to cm
        pcaspflow=datacloud.variables['PCAS2_FL'][:]
        CDPsamplingarea=0.252/100
        cdpflow=data.variables['TAS_RVSM'][:]*100*CDPsamplingarea 
        cdpflow[cdpflow==-9999.000]=np.nan
        pcaspflow[pcaspflow==-9999.000]=np.nan
        
        
        
          
        #getting pcasp data
        PCASP=np.zeros((len(timecloud),30))
        for i in range(1,31):
            number=str(i)
            if len(number)==1:
                number='0'+number
            PCASP[:,i-1]=datacloud.variables['PCAS2_'+str(number)][:]
        
        timeint=np.zeros(len(timecloud))
        for i in range(0,len(timecloud)):
            timeint[i]=int(timecloud[i])
        #arrangeing the data deppending the number of stops. It does it for both the dn and the flow
        PCASP[PCASP==-9999.000]=np.nan
        pcaspflow[pcaspflow==-9999.0]=np.nan
        
        if nstops==0:        
            interestingrangepcasp=np.array((PCASP[(start0-itime):(end0-itime),:]))
            interestingrangepcaspflow=np.array((pcaspflow[(start0-itime):(end0-itime)]))
            interestingrangepcasptot=np.array((PCASPtot[(start0-itime):(end0-itime)]))
            
        if nstops==1:        
            interestingrangepcasp=np.concatenate((PCASP[(start0-itime):(start1-itime),:],PCASP[(end1-itime):(end0-itime),:]))
            interestingrangepcaspflow=np.concatenate((pcaspflow[(start0-itime):(start1-itime)],pcaspflow[(end1-itime):(end0-itime)]))
            interestingrangepcasptot=np.concatenate((PCASPtot[(start0-itime):(start1-itime)],PCASPtot[(end1-itime):(end0-itime)]))
            
        if nstops==2:        
            interestingrangepcasp=np.concatenate((PCASP[(start0-itime):(start1-itime),:],PCASP[(end1-itime):(start2-itime),:],PCASP[(end2-itime):(end0-itime),:]))
            interestingrangepcaspflow=np.concatenate((pcaspflow[(start0-itime):(start1-itime)],pcaspflow[(end1-itime):(start2-itime)],pcaspflow[(end2-itime):(end0-itime)]))
            interestingrangepcasptot=np.concatenate((PCASPtot[(start0-itime):(start1-itime)],PCASPtot[(end1-itime):(start2-itime)],PCASPtot[(end2-itime):(end0-itime)]))
    
        interestingrangepcaspflow[interestingrangepcaspflow==-9999.000]=np.nan
        dnpcasp=np.nanmean(interestingrangepcasp,axis=0)
        PCASPflow=np.nansum(interestingrangepcaspflow,axis=0)       
        
        #this piece of code averages the bins 5,6, and 15,16 and excludes the last
        #gain stages PCASP the gain stages are from 5 to 6 and from 15 to 16
        fstage=4
        sstage=14
        newdnpcasp=np.zeros(len(dnpcasp)-3)
        newdnpcasp[:fstage]=dnpcasp[:fstage]
        newdnpcasp[fstage]=(dnpcasp[fstage]+dnpcasp[fstage+1])
        newdnpcasp[fstage+1:sstage-1]=dnpcasp[fstage+2:sstage]
        newdnpcasp[sstage-1]=(dnpcasp[sstage]+dnpcasp[sstage+1])
        newdnpcasp[sstage:]=dnpcasp[sstage+2:-1]
        
        newdlogDppcasp=np.zeros(len(dlogDppcasp)-3)
        newdlogDppcasp[:fstage]=dlogDppcasp[:fstage]
        newdlogDppcasp[fstage]=(dlogDppcasp[fstage]+dlogDppcasp[fstage+1])
        newdlogDppcasp[fstage+1:sstage-1]=dlogDppcasp[fstage+2:sstage]
        newdlogDppcasp[sstage-1]=(dlogDppcasp[sstage]+dlogDppcasp[sstage+1])
        newdlogDppcasp[sstage:]=dlogDppcasp[sstage+2:-1]
        
        newerrdlogDppcasp=np.zeros(len(errdlogDppcasp)-3)
        newerrdlogDppcasp[:fstage]=errdlogDppcasp[:fstage]
        newerrdlogDppcasp[fstage]=(errdlogDppcasp[fstage]+errdlogDppcasp[fstage+1])
        newerrdlogDppcasp[fstage+1:sstage-1]=errdlogDppcasp[fstage+2:sstage]
        newerrdlogDppcasp[sstage-1]=(errdlogDppcasp[sstage]+errdlogDppcasp[sstage+1])
        newerrdlogDppcasp[sstage:]=errdlogDppcasp[sstage+2:-1]
        
        newBincenterpcasp=np.zeros(len(Bincenterpcasp)-3)
        newBincenterpcasp[:fstage]=Bincenterpcasp[:fstage]
        newBincenterpcasp[fstage]=(Bincenterpcasp[fstage]+Bincenterpcasp[fstage+1])/2
        newBincenterpcasp[fstage+1:sstage-1]=Bincenterpcasp[fstage+2:sstage]
        newBincenterpcasp[sstage-1]=(Bincenterpcasp[sstage]+Bincenterpcasp[sstage+1])/2
        newBincenterpcasp[sstage:]=Bincenterpcasp[sstage+2:-1]
        
        newerrBincenterpcasp=np.zeros(len(errBincenterpcasp)-3)
        newerrBincenterpcasp[:fstage]=errBincenterpcasp[:fstage]
        newerrBincenterpcasp[fstage]=(errBincenterpcasp[fstage]+errBincenterpcasp[fstage+1])
        newerrBincenterpcasp[fstage+1:sstage-1]=errBincenterpcasp[fstage+2:sstage]
        newerrBincenterpcasp[sstage-1]=(errBincenterpcasp[sstage]+errBincenterpcasp[sstage+1])
        newerrBincenterpcasp[sstage:]=errBincenterpcasp[sstage+2:-1]
        
        del dnpcasp
        dlogDppcasp=0
        errdlogDppcasp=0
        Bincenterpcasp=0
        errBincenterpcasp=0
        
        dnpcasp=newdnpcasp
        dlogDppcasp=newdlogDppcasp
        errdlogDppcasp=newerrdlogDppcasp
        Bincenterpcasp=newBincenterpcasp
        errBincenterpcasp=newerrBincenterpcasp
        
        
        
        nbin=len(Bincenterpcasp)
        dndlogdppcasp=dnpcasp/dlogDppcasp
        
        
        
        #getting cdpdata
        
        CDP=np.zeros((len(timecloud),30))
        for i in range(1,31):
            number=str(i)
            if len(number)==1:
                number='0'+number
            
            CDP[:,i-1]=datacloud.variables['CDP_'+str(number)][:]
        
        timeint=np.zeros(len(timecloud))
        for i in range(0,len(timecloud)):
            timeint[i]=int(timecloud[i])
        #arrangeing the data deppending the number of stops   
        CDP[CDP==-9999.000]=np.nan
        cdpflow[cdpflow==-9999.000]=np.nan
    
        if nstops==0:        
            interestingrangecdp=CDP[(start0-itime):(end0-itime),:]
            interestingrangecdpflow=cdpflow[(start0-itime):(end0-itime)]
            interestingrangecdptot=CDPtot[(start0-itime):(end0-itime)]
                                    
            
        if nstops==1:        
            interestingrangecdp=np.concatenate((CDP[(start0-itime):(start1-itime),:],CDP[(end1-itime):(end0-itime),:]))
            interestingrangecdpflow=np.concatenate((cdpflow[(start0-itime):(start1-itime)],cdpflow[(end1-itime):(end0-itime)]))
            interestingrangecdptot=np.concatenate((CDPtot[(start0-itime):(start1-itime)],CDPtot[(end1-itime):(end0-itime)]))
            
        if nstops==2:        
            interestingrangecdp=np.concatenate((CDP[(start0-itime):(start1-itime),:],CDP[(end1-itime):(start2-itime),:],CDP[(end2-itime):(end0-itime),:]))
            interestingrangecdpflow=np.concatenate((cdpflow[(start0-itime):(start1-itime)],cdpflow[(end1-itime):(start2-itime)],cdpflow[(end2-itime):(end0-itime)]))
            interestingrangecdptot=np.concatenate((CDPtot[(start0-itime):(start1-itime)],CDPtot[(end1-itime):(start2-itime)],CDPtot[(end2-itime):(end0-itime)]))
        
        
        filtered=0
        for i in range(0,len(interestingrangecdptot)):
            if interestingrangecdptot[i]>0.3 or interestingrangepcasptot[i]>4500:
                interestingrangepcasp[i,:]=np.nan
                interestingrangecdp[i,:]=np.nan
                interestingrangepcaspflow[i]=np.nan
                interestingrangecdpflow[i]=np.nan
                filtered+=1
                interestingrangepcasptot[i]=np.nan
                interestingrangecdptot[i]=np.nan
           
               
        print('percentage of filtered data ')      
        print(float(100*filtered/len(interestingrangepcasp)))
        
        
        interestingrangecdpflow[interestingrangecdpflow==-9999.000]=np.nan
        dncdp=np.nanmean(interestingrangecdp,axis=0)
        CDPflow=np.nansum(interestingrangecdpflow,axis=0) 
        
        nbin=len(Bincentercdp)
        
        dndlogdpcdp=dncdp/dlogDpcdp  #dont forget to exclude the first bin
        #calcualting areas
        dadlogdppcasp=dndlogdppcasp*np.pi*Bincenterpcasp**2
        dadlogdpcdp=dndlogdpcdp*np.pi*Bincentercdp**2
        
        #calculating error bar, including the error in dn, the error in D, the error in dlogdp and the discretisation error. See phil rosemberg document. 
        dnpcasp[dnpcasp=='nan']=0
        dncdp[dncdp=='nan']=0
        errdnpcasp=np.sqrt(PCASPflow*dnpcasp)+(1/np.sqrt(12))
        errdncdp=np.sqrt(CDPflow*dncdp)+(1/np.sqrt(12))
        
        errdndlogdpcdp=np.sqrt(((dndlogdpcdp*errdncdp/(CDPflow*dncdp))**2)+((dndlogdpcdp*errdlogDpcdp/dlogDpcdp)**2))
        errdndlogdppcasp=np.sqrt(((dndlogdppcasp*errdnpcasp/(PCASPflow*dnpcasp))**2)+((dndlogdppcasp*errdlogDppcasp/dlogDppcasp)**2))
        
        errdadlogdpcasp=np.sqrt(((dadlogdppcasp*errdnpcasp/(PCASPflow*dnpcasp))**2)+((2*dadlogdppcasp*errBincenterpcasp/Bincenterpcasp)**2)+(dadlogdppcasp*errdlogDppcasp/dlogDppcasp)**2)
        errdadlogdcdp=np.sqrt(((dadlogdpcdp*errdncdp/(CDPflow*dncdp))**2)+((2*dadlogdpcdp*errBincentercdp/Bincentercdp)**2)+(dadlogdpcdp*errdlogDpcdp/dlogDpcdp)**2)
        
        
        sizedist=np.concatenate((dadlogdppcasp[1:-1],dadlogdpcdp[1:]))
        bins=np.concatenate((Bincenterpcasp[1:-1],Bincentercdp[1:]))
        errorsda=np.concatenate((errdadlogdpcasp[1:-1],errdadlogdcdp[1:]))
        
         
        
        #pressure and temp corrections
        time=data.variables['Time'][:]
        Temp=data.variables['TAT_ND_R'][:]
        Pres=data.variables['PS_RVSM'][:]
        Palt=data.variables['PALT_RVS'][:]
        Ralt=data.variables['HGT_RADR'][:]
        
        Presstd=1000
        Tempstd=273.17
        '''
        fig = plt.figure(figsize=(7, 8))
        ax1 = fig.add_subplot(311)
        plt.plot(time,Temp)
        plt.axvspan(start0,end0,alpha=0.5)
        pylab.ylim([250,320])
        ax2 = fig.add_subplot(312)
        plt.plot(time,Pres)
        plt.axvspan(start0,end0,alpha=0.5)
        ax3 = fig.add_subplot(313)
        plt.plot(time,Palt,color='b', label='Pressure')
        plt.plot(time,Ralt,color='r', label='Radar')
        plt.axvspan(start0,end0,alpha=0.5)
        #plt.legend(loc=1)
        pylab.ylim([0,400])
        '''
        
        Temp[Temp==-9999.000]=np.nan
        interestingrangeTemp=Temp[(start0-itime):(end0-itime)]
        Tempavg=np.nanmean(interestingrangeTemp) 
        
        Pres[Pres==-9999.000]=np.nan
        interestingrangePres=Pres[(start0-itime):(end0-itime)]
        Presavg=np.nanmean(interestingrangePres)     
          
        Correction=(Presavg/Presstd)*(Tempstd/Tempavg)
        
        dndlogdppcasp=dndlogdppcasp/Correction
        #errnumberpcasp=errnumberpcasp/Correction
        dndlogdpcdp=dndlogdpcdp/Correction
        #errnumbercdp=errnumbercdp/Correction
        
        dadlogdppcasp=dadlogdppcasp/Correction
        errdadlogdpcasp=errdadlogdpcasp/Correction
        dadlogdpcdp=dadlogdpcdp/Correction
        errdadlogdcdp=errdadlogdcdp/Correction
        
        errdndlogdpcasp=errdadlogdpcasp/Correction
        errdndlogdcdp=errdadlogdcdp/Correction
    
    
    
    #plotting 
    
    plt.figure(figsize=(12,4),dpi=300)
    plt.subplots_adjust(wspace=0.15, hspace=0)
    plt.suptitle(suptitle,fontsize=18,y=1.02)
    p1=plt.subplot(1,2,1)
    plt.errorbar(Bins,N,color='b',yerr=errN,label='SEM')
    if pcaspflag=='yes':
        plt.errorbar(Bincenterpcasp[:],dndlogdppcasp[:],yerr=errdndlogdppcasp[:],xerr=errBincenterpcasp,color='g',label='Optical Probes',marker='o',markersize=3,linestyle=' ')
        plt.errorbar(Bincentercdp[:],dndlogdpcdp[:],yerr=errdndlogdpcdp[:],xerr=errBincentercdp,color='g',marker='o',markersize=3,linestyle=' ')
    plt.xscale('log')
    plt.yscale('log')
    plt.tick_params(labelsize=14)
    #plt.legend(loc=1)
    pylab.xlim([0.1,50])
    pylab.ylim([10**-4,10**5])
    plt.xlabel('diameter ($\mathrm{\mu m}$)',size=14)
    plt.ylabel('dN/dlogDp (cm$^{-3}$)',size=14)
    f = lambda x,pos: str(x).rstrip('0').rstrip('.')
    p1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(f))
    plt.legend(frameon=False,fontsize=14,loc=1,numpoints=1)
    plt.text(-0.12, 1.1, str(letter), transform=p1.transAxes,size=16, weight='bold')
    '''
    p2=plt.subplot(1,2,2)
    plt.errorbar(Bins,A,color='b',yerr=errA,label='SEM')
    if pcaspflag=='yes':
        plt.errorbar(Bincenterpcasp[:],dadlogdppcasp[:],color='g',yerr=errdadlogdpcasp[:],xerr=errBincenterpcasp,linestyle=' ',marker='o',markersize=3,label='Optical Probes')
        plt.errorbar(Bincentercdp[:],dadlogdpcdp[:],color='g',yerr=errdadlogdcdp[:],xerr=errBincentercdp,linestyle=' ',marker='o',markersize=3)
    plt.xscale('log')
    #plt.yscale('log')
    
    plt.legend(loc=1)
    pylab.xlim([0.2,50])
    plt.gca().set_ylim(bottom=0)
    if ymax>0:
        plt.gca().set_ylim(top=ymax)
    #pylab.ylim([0,250])
    plt.tick_params(labelsize=14)
    plt.xlabel('Diameter ($\mathrm{\mu m}$)',size=14)
    plt.ylabel('dA/dlogDp ($\mu m$$^{2}$ cm$^{-3}$)',size=14)
    f = lambda x,pos: str(x).rstrip('0').rstrip('.')
    p1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(f))
    plt.legend(frameon=False,fontsize=14,loc=1,numpoints=1)
    #plt.suptitle=title
    '''


    path=path_chem
    data=np.genfromtxt(path,delimiter=',',skip_header=1,dtype=str)
    ECDm=[]
    Class=[]
    counter=0
    for i in range(0,len(data)): 
        if not data[i,1]=='' and data[i,1]!='Class' and data[i,9]!='':
            ECDm.append(data[i,6])
            Class.append(data[i,1])
            counter=counter+1
            

    totals=np.zeros(len(BinsC)-1)
    cath=10
    finaldata=np.zeros((cath,len(BinsC)-1))
    for j in range(0,len(BinsC)-1):
        types=[]
        for i in range(0,len(ECDm)):
            if float(ECDm[i])>=BinsC[j] and float(ECDm[i])<BinsC[j+1] and not (Class[i]=='Metallic Cr'):
                totals[j]+=1
                types.append(Class[i])
            
        if not types==[]:       
                    
            finaldata[0,j]=float(types.count('No oxigen')+types.count('All')+types.count('Mainly O')+types.count('No Classification'))/len(types)
                    
            finaldata[1,j]=float(types.count('Oxygen only')+types.count('Oxigen+Cu(trace)')+types.count('Oxygen +K')+types.count('P rich')+types.count( 'False Si '))/len(types)   
            
            finaldata[2,j]=float(types.count('Sulfate aerosol'))/len(types) 
            
            finaldata[3,j]=float(types.count('Cl')+types.count('Cl+K'))/len(types) 
    
            finaldata[4,j]=float(types.count('Aged Sea Aerosol')+types.count('Sea aerosol (Cl)')+types.count('Sea aerosol'))/len(types) 
            
            finaldata[5,j]=float(types.count('Metallic Fe')+types.count('Metallic Ti')+types.count('Metallic Cu')+types.count('Metallic Mn')+types.count('Metallic Mg')+types.count('Metallic Al')+types.count('False Cr')+types.count('Metallic Zn')+types.count('Metallic Pb'))/len(types)   #+types.count('Metallic Cr')
            
            finaldata[6,j]=float(types.count('Ca rich')+types.count('Ca rich+NaCl')+types.count('Ca pure')+types.count('Gypsium'))/len(types)
            
            finaldata[7,j]=float(types.count('Aluminosillicates')+types.count('Aluminosillicates+NaCl'))/len(types) 
            
            finaldata[8,j]=float(types.count('Sillica'))/len(types)
            
            finaldata[9,j]=float(types.count('Sillica mixtures ')+types.count('Silica mixtures +NaCl'))/len(types)
            
            
                    
    xnames=1.05*(np.linspace(0,len(BinsC)-1,len(BinsC)-1)-0.5)
    xticks=np.zeros(len(BinsC)-1)
      
    p3=plt.subplot(1,2,2)
        
    K, N = finaldata.shape
    ind = np.arange(N)
    width = 1.08
    
    plt.bar(xnames, finaldata[0,:])
    
    plots = []
    height_cumulative = np.zeros(N)
    
    patterns = [ "\\" , "|" , "/","\\" , "|" , "/","\\" , "|" , "/","\\" , "|" , "/"  ]
    colors=['green','lightgrey','yellow','palevioletred','cornflowerblue','firebrick','sandybrown','orange','tan','khaki']
    plt.rcParams['hatch.linewidth'] = 0.3
    
    for k in range(K):
        color = plt.cm.Paired((K-k)/float(K), 1)
        color=colors[k]
        if k == 0:
            p = plt.bar(xnames, finaldata[k,:], width, color=color,linewidth=0.4,hatch=patterns[k])
    
        else:
            p = plt.bar(xnames, finaldata[k,:], width, bottom=height_cumulative,color=color,linewidth=0.4,hatch=patterns[k])
        height_cumulative += finaldata[k,:]
        plots.append(p)    
    plt.ylim((0, 1))
    plt.ylabel('proportion',size=14)
    plt.xlabel('diameter ($\mathrm{\mu m}$)',size=14)
    
    #plt.title(title)
    total=(BinsC)
    
        
    finaldatabulk=finaldata*totals  
    errfinaldatabulk=np.sqrt(finaldatabulk)
    errfinaldata=errfinaldatabulk/totals
    
    
    for i in range(0,len(totals)):
        totals[i]=str(totals[i])
      
    for i in range(0,N):
        a=str(totals[i])
        plt.annotate(a[:-2],(xnames[i]-0.475,0.95),fontsize=12)
    plt.annotate('1',(15,0.05))
    
    for i in xnames:
        plt.axvline(x=i-0.575,color='black',linewidth=2)
        
    for i in range(0,len(BinsC)):
        BinsC[i]=round(BinsC[i],1)
        
        
    xnames=np.linspace(0,len(BinsC),len(BinsC))+0.5
    
    plt.axvline(xnames[-1]-1.00,color='black',linewidth=2) #+0.11
    plt.xticks(1.05*xnames-1.55,BinsC)
    
    print(np.sum(finaldata,axis=0))
    print('this should be all 1s')
    plt.yticks(np.linspace(0,1,5))
    
    plt.tick_params(axis='x',which='both',bottom=False,top=False) 
    if legend=='yes':
        topic_labels = ['Other','Carbonaceous','S rich','Cl rich','Na rich','Metal rich','Ca rich','Al-Si rich','Si only','Si rich']
        lgd=plt.legend([p[0] for p in plots[::-1]], topic_labels[::-1],bbox_to_anchor=(1.05, 1), loc=2,borderaxespad=0.,fontsize=13)
    plt.text(13.5, 0.02, 'N= '+str(np.sum(totals))[:-2], fontsize=16)
    #plt.title(title,size=14)
    #pylab.xlim([0.1,13])
    
    pylab.savefig('\\\\envdfs.leeds.ac.uk\\a86\py15asm\IV paper\Codes\SEM_graphs\\' + title+'.png',bbox_inches='tight',dpi = 300)


    plt.clf()

    
    
    #non classified figures
    nonclas=0
    for i in range(0,len(Class)):
        if Class[i]!='No oxigen' and Class[i]!='No oxigen' and Class[i]!='All' and Class[i]!='Mainly O' and Class[i]!='Oxygen only' and Class[i]!='Oxigen+Cu(trace)' and Class[i]!='Oxygen +K' and Class[i]!='P rich' 'No Classification' and Class[i]!='False Si ' and Class[i]!='Sulfate aerosol' and Class[i]!='Metallic Fe' and Class[i]!='Metallic Ti' and Class[i]!='Metallic Cu' and Class[i]!='Metallic Mn' and Class[i]!='Metallic Al' and Class[i]!='False Cr' and Class[i]!='Metallic Zn' and Class[i]!='Metallic Pb' and Class[i]!='Aged Sea Aerosol' and Class[i]!='Sea aerosol (Cl)' and Class[i]!='Sea aerosol' and Class[i]!='Cl' and Class[i]!='Cl+K' and Class[i]!='Ca rich' and Class[i]!='Ca rich+NaCl' and Class[i]!='Ca pure' and Class[i]!='Gypsium' and Class[i]!='Aluminosillicates' and Class[i]!='Aluminosillicates+NaCl' and Class[i]!='Sillica' and Class[i]!='Sillica mixtures ' and Class[i]!='Silica mixtures +NaCl':
                nonclas+=1
                print(Class[i])
    print('non-clasified particles: '+str(nonclas))

    
    
    #dust area
    N=np.concatenate((Nhist__highmag,Nhist__lowmag))
    plt.figure(figsize=(6,4.5),dpi=300)
    dustcompositions=np.sum(finaldata[6:,:],axis=0)
    errdustcompositions=np.sqrt(np.sum(finaldatabulk[6:,:],axis=0))/totals
    distbins=BinsC
    dustA=np.zeros(len(A))
    errdustA=np.zeros(len(A))
    dustN=np.zeros(len(A))

    
    
    for i in range(0,len(A)):
        for j in range(0,len(dustcompositions)):
            if Bins[i]>BinsC[j] and Bins[i]<BinsC[j+1]:
                print(i)
                print(j)
                dustA[i]=float(A[i]*dustcompositions[j])
                dustN[i]=float(N[i]*dustcompositions[j])
                errdustA[i]=np.sqrt(float(dustcompositions[j]*errA[i])**2+float(errdustcompositions[j]*A[i])**2)
                
                distbins=BinsC
                
                
    #NaCl area
    NaClcompositions=finaldata[4]
    errNaClcompositions=np.sqrt(finaldatabulk[4])/totals
    NaClbins=BinsC
    NaClA=np.zeros(len(A))
    NaClN=np.zeros(len(A))
    errNaClA=np.zeros(len(A))

    
    
    for i in range(0,len(A)):
        for j in range(0,len(NaClcompositions)):
            if Bins[i]>BinsC[j] and Bins[i]<BinsC[j+1]:
    
                NaClA[i]=float(A[i]*NaClcompositions[j])
                NaClN[i]=float(N[i]*NaClcompositions[j])
                errNaClA[i]=np.sqrt(float(NaClcompositions[j]*errA[i])**2+float(errNaClcompositions[j]*A[i])**2)
                
                NaClbins=BinsC
                
                
                
                
    plt.errorbar(Bins,A,color='b',yerr=errA,label='total')
    plt.errorbar(Bins,dustA,color='g',yerr=errdustA,label='mineral dust')                
    plt.errorbar(Bins,NaClA,color='r',yerr=errNaClA,label='NaCl')
        
    #plt.title('Barbados 170803')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend(loc=1)
    pylab.xlim([0.1,50])
    plt.gca().set_ylim(bottom=0)
    
    plt.xlabel('size ($\mu m$)')
    plt.ylabel('dA/dlogDp ($\mu m$$^{2}$ cm$^{-3}$)')
    
    plt.legend()
    plt.savefig(title+'dust_area.png',dpi = 300)
    
    #integrate dust 
    UpdustA=dustA+errdustA
    LowdustA=dustA-errdustA
    LowdustA[LowdustA<0]=0
    
    dustA[np.isnan(dustA)]=0
    LowdustA[np.isnan(LowdustA)]=0
    UpdustA[np.isnan(UpdustA)]=0
    A[np.isnan(A)]=0
    
    Iup=np.sum(UpdustA*dlogDp)
    Ilow=np.sum(LowdustA*dlogDp)
    I=np.sum(dustA*dlogDp)
    Itotal=np.sum(A*dlogDp)
    print(Iup)#*10**-6)
    print(I)#*10**-6)
    print(Ilow)#*10**-6)
    print('precetnage of dust = '+str(100.0*I/Itotal))
    IN=np.sum(dustN*dlogDp)
    print('number of dust particles')
    print(IN)
    
    
    #integrate NaCl
    UpNaClA=NaClA+errNaClA
    LowNaClA=NaClA-errNaClA
    LowNaClA[LowNaClA<0]=0
    
    NaClA[np.isnan(NaClA)]=0
    LowNaClA[np.isnan(LowNaClA)]=0
    UpNaClA[np.isnan(UpNaClA)]=0
    A[np.isnan(A)]=0
    
    Iup=np.sum(UpNaClA*dlogDp)
    Ilow=np.sum(LowNaClA*dlogDp)
    I=np.sum(NaClA*dlogDp)
    Itotal=np.sum(A*dlogDp)
    print(Iup)#*10**-6
    print(I)#*10**-6
    print(Ilow)#*10**-6
    print('precetnage of NaCl = '+str(100.0*I/Itotal))
    
    IN=np.sum(NaClN*dlogDp)
    print('number of NaCl particles')
    print(IN)
    
    
    
    return counter_tot



#%%
path__lowmag='\\\\envdfs.leeds.ac.uk\\a86\py15asm\\SEM exports\\180316_FAAM_2045_2126_0.51V\Morph\\180316_FAAM_2045_2126_0.51V_morph_lowmag.csv'
path__highmag='\\\\envdfs.leeds.ac.uk\\a86\py15asm\\SEM exports\\180316_FAAM_2045_2126_0.51V\Morph\\180316_FAAM_2045_2126_0.51V_morph_highmag.csv'

lowlim=0.3
medlim=0.8
uplim=10


air=(2215.0/60)*75*0.51*1000 #2215 secs

area__highmag=104*93+78.5*75+86.2*82+86.7*80+86.5*82.6+85*77

area__lowmag=area__highmag+563*528+686*643+747*758+747*758+747*758

path='\\\\envdfs.leeds.ac.uk\\a86\py15asm\MACSSIMIZE\CEDA data\\C087 180316\\core_faam_20180316_v004_r0_c087_1hz.nc'
pathcloud='\\\\envdfs.leeds.ac.uk\\a86\py15asm\MACSSIMIZE\CEDA data\\C087 180316\\core-cloud-phy_faam_20180316_v501_r0_c087.nc'
pathcsv='\\\\envdfs.leeds.ac.uk\\a86\py15asm\MACSSIMIZE\Data\\180316\\EF600\\180316_FAAM_2044_2126_LOW_OP_565L_Tef0.45\\info.csv'
data=Dataset(path)
datacloud=Dataset(pathcloud)
#datacsv=np.genfromtxt(pathcsv,delimiter=',',dtype=str)
pcaspflag='yes'
pathpcasp='\\\\envdfs.leeds.ac.uk\\a86\py15asm\MACSSIMIZE\CEDA data\Calibrations\PCASP2\MACSSIMIZE PCASP2 180504 calibration, n=1.56.csv'
pathcdp='\\\\envdfs.leeds.ac.uk\\a86\py15asm\MACSSIMIZE\CEDA data\Calibrations\CDP\\master MACSSIMIZE 180226 CDP calibrated diams n=1.56.csv'

path_chem='\\\\envdfs.leeds.ac.uk\\a86\py15asm\\SEM exports\\180316_FAAM_2045_2126_0.51V\Chem\\180316_FAAM_2045_2126_0.51V_chem.csv'
os.chdir('\\\\envdfs.leeds.ac.uk\\a86\py15asm\SEM exports\\180316_FAAM_2045_2126_0.51V')
title='180316 at 2045 to 2146'



#%% repeat from here but using another dataset

### inputs ###
path_highmag=r"C:\Users\py15asm\Documents\GitHub\ice-nucleation-leeds-SEM\data\chemical_processed\180316_FAAM_2045_2126_0.51V_morph_highmag_processed.csv"
path_lowmag=r"C:\Users\py15asm\Documents\GitHub\ice-nucleation-leeds-SEM\data\chemical_processed\180316_FAAM_2045_2126_0.51V_morph_lowmag_processed.csv"
area_highmag=104*93+78.5*75+86.2*82+86.7*80+86.5*82.6+85*77
area_lowmag=area_highmag+563*528+686*643+747*758+747*758+747*758
area_filter=np.pi*(40*1000/2)**2 #calculaton based on 47 mm
area_filter=11*10**8 #value by HP in um
air=(2215.0/60)*75*0.51*1000 #in cm3
lowlim=0.3
medlim=0.8
uplim=10
### importing and basic formatting ###

data_highmag=pd.read_csv(path_highmag)
data_lowmag=pd.read_csv(path_lowmag)

mag_threshold=0.8
#exclude Cr
data_highmag_processed=data_highmag[data_highmag["Class"]!="Metallic Cr"][["Id","ECD (um)"]]
data_lowmag_processed=data_lowmag[data_lowmag["Class"]!="Metallic Cr"][["Id","ECD (um)"]]

data_lowmag_processed=data_lowmag_processed.append(data_highmag_processed)

#data_highmag_processed=data_highmag_processed[data_highmag_processed["ECD (um)"]<mag_threshold]
#data_lowmag_processed=data_lowmag_processed[data_lowmag_processed["ECD (um)"]>mag_threshold]


### bins ###

bins=np.logspace(np.log10(0.3),np.log10(10),11)        
bins[:4]=np.logspace(np.log10(lowlim),np.log10(0.8),4) 
bins[3:]=np.logspace(np.log10(0.8),np.log10(10),8)
bins=bins.round(1)

bins_highmag=bins[:4]
bins_lowmag=bins[3:]
nbins_highmag=len(bins_highmag)
nbins_lowmag=len(bins_lowmag)

bincentre_highmag=np.zeros(len(bins_highmag)-1)
bincentre_lowmag=np.zeros(len(bins_lowmag)-1)
dlogDp_highmag=np.zeros(len(bins_highmag)-1)
dlogDp_lowmag=np.zeros(len(bins_lowmag)-1)

for i in range(0,len(bins_highmag)-1):
    bincentre_highmag[i]=(bins_highmag[i]+bins_highmag[i+1])/2
    dlogDp_highmag[i]=math.log10(bins_highmag[i+1])-math.log10(bins_highmag[i])
    
for i in range(0,len(bins_lowmag)-1):
    bincentre_lowmag[i]=(bins_lowmag[i]+bins_lowmag[i+1])/2
    dlogDp_lowmag[i]=math.log10(bins_lowmag[i+1])-math.log10(bins_lowmag[i])
    
    
### calculations ###
N_highmag,edges=np.histogram(data_highmag_processed['ECD (um)'],bins=bins_highmag)
N_lowmag,edges=np.histogram(data_lowmag_processed['ECD (um)'],bins=bins_lowmag)

N_highmag_err=np.sqrt(N_highmag)
N_lowmag_err=np.sqrt(N_lowmag)

dNdlogDp_highmag=N_highmag*area_filter/(air*area_highmag*dlogDp_highmag)
dNdlogDp_lowmag=N_lowmag*area_filter/(air*area_lowmag*dlogDp_lowmag)

dAdlogDp_highmag=dNdlogDp_highmag*np.pi*bincentre_highmag**2
dAdlogDp_lowmag=dNdlogDp_lowmag*np.pi*bincentre_lowmag**2

#concatenating
dAdlogDp=np.concatenate([dAdlogDp_highmag,dAdlogDp_lowmag])
bincentre=np.concatenate([bincentre_highmag,bincentre_lowmag])
dlogDp=np.concatenate([dlogDp_highmag,dlogDp_lowmag])

#plot
plt.plot(bincentre,dAdlogDp)
plt.xscale('log')

### note from here you need to run a dev chem to have finaldata
dust_frac=finaldata[5:,:].sum(axis=0)
SS_frac=finaldata[3,:]

dust_area=np.sum(dust_frac*dAdlogDp*dlogDp)
area=np.sum(dAdlogDp*dlogDp)
print(dust_area)
print(100*dust_area/area)

SS_area=np.sum(SS_frac*dAdlogDp*dlogDp)
print(SS_area)
print(100*SS_area/area)

