#import basic libraries for plotting, data structures and signal processing
import matplotlib.pyplot as plt
import numpy as np
import pandas
from scipy.signal import find_peaks, argrelmin
from scipy.interpolate import CubicSpline
import imageio
import os
import fnmatch
import datetime


folder='Test'
#fpath = '/Users/alkadas/Desktop/Data/IPD/'+folder+'/'        #filepath for Mac where the data is
fpath = 'C:/Users/alaka/Google Drive/IPD/'+folder+'/'
imgfiles = fnmatch.filter(os.listdir(fpath), '*.tif')

#PARAMETERS
f = os.path.basename(__file__)                  #stores the filename of the current script in f
toa = str(datetime.datetime.today()).split()
today = toa[0]
now = toa[1]
mu_per_px = 0.18825     #pixels to microns conversion factor
order=10                #set order parameter for argrelmin for calculating baseline
prom = 5                #set prominence parameter for find_peaks
df_Parameters = pandas.DataFrame(data={'1. Date of analysis':today, '2. Time of analysis':now, '3. Microns per pixel':mu_per_px, '4. Script name':f, '5. argrelmin order':order, '6. find_peaks prominence':prom}, index=[0])

#specify columns of the pandas dataframe and excel sheets
cols_Data =     ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance','Raw intensity', 'Background intensity', 'Neurite intensity', 'Baseline']
cols_Peaks =    ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance', 'Punctum max intensity', 'Punctum width']
cols_IPDs =     ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance', 'Inter-punctum interval']
cols_Analysis = ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Image size', 'Max neurite length', 'Average neurite intensity','Total peaks', 'Average peak intensity', 'Average peak width', 'Average ipd', 'Median ipd']

#initialize Pandas DataFrames
df_Data = pandas.DataFrame()
df_Peaks = pandas.DataFrame()
df_IPDs = pandas.DataFrame()
df_Analysis = pandas.DataFrame()

for x in imgfiles:                          #create loop for number of images in folder
    img = imageio.imread(fpath+x)[:,:,1]    #import image and store it in a list of lists
    
    #extract info from filename
    date=x.split('_')[0]
    strain = folder
    neuron = x.split('_')[3][:3]
    lr = x.split('_')[3][3]

    imsize = np.shape(img)                  #calculate image size
    d=np.arange(imsize[1])                  #create list of integers from 0 to length of image for x-axis
    dist = d*mu_per_px                      #pixel to microns conversion
    normdist=dist/dist[-1]
    n = img[8:13, 0:]                       #extract rows to use for neurite
    bg = np.concatenate((img[0:6, 0:], img[15: , 0:]))   #extract rows to use for background
    rawf = np.mean(n, axis=0)               #calculate average raw neurite fluorescence
    bgf = np.mean(bg, axis=0)               #calculate average background fluorescence
    nf = rawf - bgf                         #calculate background subtracted neurite fluorescence
    
    
    #finding baseline
    minsp = np.concatenate((0, argrelmin(nf, order=order)[0], imsize[1]), axis=None)
    minsh = np.concatenate((nf[minsp[1]], [nf[i] for i in minsp[1:-1]], nf[minsp[-2]]), axis=None)
    #set negative minsh values to zero
    for i in range(len(minsh)):
        if minsh[i]<0: minsh[i]=0
    cs = CubicSpline(minsp, minsh, extrapolate=False)(d)
    
    plt.figure(1, figsize=(0.015*imsize[1],15))
    plt.subplot(311)
    plt.title(x+' Finding baseline')
    plt.xlabel('Distance (um)')
    plt.ylabel('Intensity (AU)')
    plt.axis([0, max(dist), min(nf), 50])
    plt.plot(dist, nf, 'y-')
    plt.plot(minsp*mu_per_px, minsh, 'ro', dist, cs, 'k-')
    
    
    #add image data to pandas dataframe
    all_data1 = pandas.DataFrame({'Date':[date]*imsize[1], 'Strain':[strain]*imsize[1], 'Neuron':[neuron]*imsize[1], 'L/R':[lr]*imsize[1], 'ImageID':[x]*imsize[1], 'Distance':dist, 'Normalized distance':normdist, 'Raw intensity':rawf, 'Background intensity':bgf, 'Neurite intensity':nf, 'Baseline':cs}, columns=cols_Data)
    df_Data=df_Data.append(all_data1)

    #baseline subtracted trace with no negative values
    bsub = nf-cs
    for i in range(len(bsub)):
        if bsub[i]<0: bsub[i]=0

    #find peaks
    peaks = find_peaks(bsub, height=0, prominence=prom, width=0, rel_height=0.5)
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
    plt.axis([0, max(dist), min(nf), 50])
    plt.plot(dist, nf, 'y-')
    plt.plot(pd, pmi, 'go')
    #create subplot3
    plt.subplot(313)
    plt.title(x+' peak widths')
    plt.xlabel('Distance (um)')
    plt.ylabel('Baseline subtracted intensity (AU)')
    plt.axis([0, max(dist), 0, 50])
    plt.plot(dist, bsub, 'y-')
    plt.hlines(peaks[1]['width_heights'], peaks[1]['left_ips']*mu_per_px, peaks[1]['right_ips']*mu_per_px)
    plt.show()
    plt.close()
        
    #add data to pandas dataframe
    all_data2 = pandas.DataFrame({'Date':[date]*len(pd), 'Strain':[strain]*len(pd), 'Neuron':[neuron]*len(pd), 'L/R':[lr]*len(pd), 'ImageID':[x]*len(pd), 'Distance':pd, 'Normalized distance':pnd, 'Punctum max intensity':pmi, 'Punctum width':peaks[1]['widths']*mu_per_px}, columns=cols_Peaks)
    df_Peaks=df_Peaks.append(all_data2)
    all_data3 = pandas.DataFrame({'Date':[date]*len(ipd), 'Strain':[strain]*len(ipd), 'Neuron':[neuron]*len(ipd), 'L/R':[lr]*len(ipd), 'ImageID':[x]*len(ipd), 'Distance':ipdd, 'Normalized distance':ipdnd, 'Inter-punctum interval':ipd}, columns=cols_IPDs)
    df_IPDs=df_IPDs.append(all_data3)    
    frame = pandas.DataFrame([[date, strain, neuron, lr, x, imsize[1], dist[-1], np.mean(nf), len(pd), np.mean(pmi), np.mean(peaks[1]['widths']*mu_per_px), np.mean(ipd), np.median(ipd)]], columns=cols_Analysis)
    df_Analysis = df_Analysis.append(frame)

#save data to excel file
#dfpath = '/Users/alkadas/Desktop/Data/IPD/'        #filepath where the data is
dfpath = 'C:/Users/alaka/Google Drive/IPD/Analysis_20200311/'
timestamp = today.replace('-','')+'-'+now.replace(':','')[:6]
wb = pandas.ExcelWriter(dfpath+timestamp+'_'+folder+'_Data.xlsx', engine='xlsxwriter')
df_Data.to_excel(wb, sheet_name='Data')
df_Peaks.to_excel(wb, sheet_name='Peaks')
df_IPDs.to_excel(wb, sheet_name='IPDs')
df_Analysis.to_excel(wb, sheet_name='Analysis')
df_Parameters.to_excel(wb, sheet_name='Parameters')
wb.save()