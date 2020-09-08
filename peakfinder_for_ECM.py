#import basic libraries for plotting, data structures and signal processing
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas
from scipy.signal import find_peaks, argrelmin
from scipy.interpolate import CubicSpline
import imageio
import os
import fnmatch
import datetime

#%%
strain_key=pandas.DataFrame({('x87a', 'WT', 'NID-1::wSc'),
                             ('x78', 'mec-4(u253)', 'NID-1::wSc'),
                             ('GN965', 'WT', 'LAM-2::mNG'),
                             ('x100','mec-1(e1526)','LAM-2::mNG')}, columns=['Strain','Allele','Label'])

def strainkey(strain):
    allele = strain_key['Allele'][(strain_key['Strain']==strain)]
    label = strain_key['Label'][(strain_key['Strain']==strain)]
    return allele, label

#%%
fpath = 'C:/Users/alaka/Work_GDrive/ECM manuscript/raw data/Laminin and nidogen puncta distribution/straightened/'
dfpath = 'C:/Users/alaka/Work_GDrive/ECM manuscript/raw data/Laminin and nidogen puncta distribution/peakfinder_output/'
imgfiles = fnmatch.filter(os.listdir(fpath), '*.tif')

#PARAMETERS
f = os.path.basename(__file__)                  #stores the filename of the current script in f
toa = str(datetime.datetime.today()).split()
today = toa[0]
now = toa[1]
mu_per_px = 0.126     #pixels to microns conversion factor
prom = 3                #set prominence parameter for find_peaks
df_Parameters = pandas.DataFrame(data={'1. Date of analysis':today, '2. Time of analysis':now, '3. Microns per pixel':mu_per_px, '4. Script name':f, '5. find_peaks prominence':prom}, index=[0])

#specify columns of the pandas dataframe and excel sheets
cols_Data =     ['Date', 'Strain', 'Allele', 'Label', 'Neuron', 'ImageID', 'Distance', 'Normalized distance','Raw intensity', 'Background intensity', 'Neurite intensity']
cols_Peaks =    ['Date', 'Strain', 'Allele', 'Label', 'Neuron', 'ImageID', 'Distance', 'Normalized distance', 'Punctum max intensity', 'Punctum width']
cols_IPDs =     ['Date', 'Strain', 'Allele', 'Label', 'Neuron', 'ImageID', 'Distance', 'Normalized distance', 'Inter-punctum interval']
cols_Analysis = ['Date', 'Strain', 'Allele', 'Label', 'Neuron', 'ImageID', 'Image size', 'Max neurite length', 'Average neurite intensity','Total peaks', 'Average peak intensity', 'Average peak width', 'Average ipd', 'Median ipd']

#initialize Pandas DataFrames
df_Data = pandas.DataFrame()
df_Peaks = pandas.DataFrame()
df_IPDs = pandas.DataFrame()
df_Analysis = pandas.DataFrame()

#%%
for x in imgfiles:                          #create loop for number of images in folder
    img = imageio.imread(fpath+x)    #import image and store it in a list of lists
    
    #extract info from filename
    date=x.split('_')[0]
    strain = x.split('_')[1].split('-')[0]
    allele = strainkey(strain)[0][strainkey(strain)[0].index[0]]
    label =  strainkey(strain)[1][strainkey(strain)[0].index[0]]
    neuron = 'ALM'

    imsize = np.shape(img)                  #calculate image size
    d=np.arange(imsize[1])                  #create list of integers from 0 to length of image for x-axis
    dist = d*mu_per_px                      #pixel to microns conversion
    normdist=dist/dist[-1]
    n = img[7:14, 0:]                       #extract rows to use for neurite
    bg = np.concatenate((img[0:6, 0:], img[15: , 0:]))   #extract rows to use for background
    rawf = np.mean(n, axis=0)               #calculate average raw neurite fluorescence
    bgf = np.mean(bg, axis=0)               #calculate average background fluorescence
    nf = rawf - bgf                         #calculate background subtracted neurite fluorescence

    for i in range(len(nf)):
        if nf[i]<0: nf[i]=0
    
    plt.figure(1, figsize=(0.010*imsize[1],15))
    plt.subplot(311)
    plt.imshow(img, vmin=0,vmax=180)
    
   
    #add image data to pandas dataframe
    all_data1 = pandas.DataFrame({'Date':[date]*imsize[1], 'Strain':[strain]*imsize[1], 'Allele':[allele]*imsize[1], 'Label':[label]*imsize[1], 'Neuron':[neuron]*imsize[1], 'ImageID':[x]*imsize[1], 'Distance':dist, 'Normalized distance':normdist, 'Raw intensity':rawf, 'Background intensity':bgf, 'Neurite intensity':nf}, columns=cols_Data)
    df_Data=df_Data.append(all_data1)


    #find peaks
    peaks = find_peaks(nf[:-397], height=0, prominence=prom, width=0, rel_height=0.5) #Do not use the last 50 um (~397 pixels) for finding peaks
    pd = peaks[0]*mu_per_px
    pnd = pd/dist[-1]
    pmi = [nf[i] for i in peaks[0]]
    ipd = np.diff(pd)
    ipdd = [pd[i]+ipd[i]/2 for i in range(0,len(ipd))]
    ipdnd = ipdd/dist[-1]
#    slope, intercept, r_value, p_value, std_err = stats.linregress(ipdnd,ipd)
    
    
    #create subplot 2
    plt.subplot(312)
    plt.title(x+' peaks')
    plt.xlabel('Distance (um)')
    plt.ylabel('Intensity (AU)')
    plt.axis([0, max(dist), min(nf), 30])
    plt.plot(dist, nf, 'y-')
    plt.plot(pd, pmi, 'go')
    #create subplot3
    plt.subplot(313)
    plt.title(x+' peak widths')
    plt.xlabel('Distance (um)')
    plt.ylabel('Baseline subtracted intensity (AU)')
    plt.axis([0, max(dist), 0, 30])
    plt.plot(dist, nf, 'y-')
    plt.hlines(peaks[1]['width_heights'], peaks[1]['left_ips']*mu_per_px, peaks[1]['right_ips']*mu_per_px)
    plt.savefig(dfpath+x[:-4]+'_pf.png')
    plt.close()
        
    #add data to pandas dataframe
    all_data2 = pandas.DataFrame({'Date':[date]*len(pd), 'Strain':[strain]*len(pd), 'Allele':[allele]*len(pd), 'Label':[label]*len(pd), 'Neuron':[neuron]*len(pd), 'ImageID':[x]*len(pd), 'Distance':pd, 'Normalized distance':pnd, 'Punctum max intensity':pmi, 'Punctum width':peaks[1]['widths']*mu_per_px}, columns=cols_Peaks)
    df_Peaks=df_Peaks.append(all_data2)
    all_data3 = pandas.DataFrame({'Date':[date]*len(ipd), 'Strain':[strain]*len(ipd), 'Allele':[allele]*len(ipd), 'Label':[label]*len(ipd), 'Neuron':[neuron]*len(ipd), 'ImageID':[x]*len(ipd), 'Distance':ipdd, 'Normalized distance':ipdnd, 'Inter-punctum interval':ipd}, columns=cols_IPDs)
    df_IPDs=df_IPDs.append(all_data3)    
    frame = pandas.DataFrame([[date, strain, allele, label, neuron, x, imsize[1], dist[-1], np.mean(nf), len(pd), np.mean(pmi), np.mean(peaks[1]['widths']*mu_per_px), np.mean(ipd), np.median(ipd)]], columns=cols_Analysis)
    df_Analysis = df_Analysis.append(frame)
#%%
plt.figure(2, figsize=(8,8))
sns.boxplot(x='Label', y='Average ipd', hue = 'Allele', data=df_Analysis, color='.8')
sns.stripplot(x='Label', y='Average ipd', hue= 'Allele', data=df_Analysis, alpha=0.5)

#%%
#save data to excel file
timestamp = today.replace('-','')+'-'+now.replace(':','')[:6]
wb = pandas.ExcelWriter(dfpath+timestamp+'_'+strain+'_Data.xlsx', engine='xlsxwriter')
df_Data.to_excel(wb, sheet_name='Data')
df_Peaks.to_excel(wb, sheet_name='Peaks')
df_IPDs.to_excel(wb, sheet_name='IPDs')
df_Analysis.to_excel(wb, sheet_name='Analysis')
df_Parameters.to_excel(wb, sheet_name='Parameters')
wb.save()