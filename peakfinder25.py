#import basic libraries for plotting, data structures and signal processing
import matplotlib.pyplot as plt
import numpy as np
import pandas
from scipy.signal import find_peaks, argrelmin,argrelmax
from scipy.ndimage import maximum_filter
from scipy import stats
import imageio
import os
import fnmatch
import seaborn as sns
import datetime

toa = str(datetime.datetime.today()).split()
today = toa[0]
now = toa[1]

fpath = '/Users/alkadas/Desktop/Data/IPD/GN753/'        #filepath where the data is
#fpath = 'C:/Users/alaka/Google Drive/IPD/GN753/'     #filepath where the data is
imgfiles = fnmatch.filter(os.listdir(fpath), '*.tif')

#PARAMETERS
mu_per_px = 0.18825     #pixels to microns conversion factor
smooth = 2
height = 0.2
prom = 0.1
df_Parameters = pandas.DataFrame(data={'Date of analysis':today, 'Time of analysis':now, 'Microns per pixel':mu_per_px, 'Smoothing type':'Maximum filter', 'Smoothing footprint':smooth, 'Height threshold':height, 'Prominence threshold':prom}, index=[0])

#specify columns of the pandas dataframe and excel sheets
cols_Data =     ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance','Raw intensity', 'Background intensity', 'Neurite intensity', 'Neurite signal-to-noise ratio', 'Filtered snr']
cols_Peaks =    ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance', 'Punctum signal-to-noise ratio', 'Punctum width']
cols_IPDs =     ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance', 'Inter-punctum interval']
cols_Analysis = ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Image size', 'Max neurite length', 'Average neurite intensity', 'Average snr', 'Signal strength','Status','Total peaks', 'Average peak snr', 'Average peak width', 'Average ipd', 'Median ipd', 'Average peaks per distance', 'slope', 'intercept', 'r_value', 'p_value', 'std_err']

#initialize Pandas DataFrames
df_Data = pandas.DataFrame()
df_Peaks = pandas.DataFrame()
df_IPDs = pandas.DataFrame()
df_Analysis = pandas.DataFrame()

for x in imgfiles:                          #create loop for number of images in folder
    img = imageio.imread(fpath+x)[:,:,1]    #import image and store it in a list of lists

    date=x.split('_')[0]
    strain = x.split('_')[1]
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
    for i in range(imsize[1]):
        if nf[i]<0: nf[i]=0                 #set any negative nf value to zero
    snr = nf/bgf                            #calculate normalized fluorescence

    strength = np.mean([snr[i] for i in argrelmax(snr,order=10)[0]])-np.mean([snr[i] for i in argrelmin(snr,order=10)[0]])          #just a quantity to judge quality of signal. Not used for any further calculations
    
    fsnr=maximum_filter(snr, smooth)        #reduce noise by applying a maximum filter

    #add image data to pandas dataframe
    all_data1 = pandas.DataFrame({'Date':[date]*imsize[1], 'Strain':[strain]*imsize[1], 'Neuron':[neuron]*imsize[1], 'L/R':[lr]*imsize[1], 'ImageID':[x]*imsize[1], 'Distance':dist, 'Normalized distance':normdist, 'Raw intensity':rawf, 'Background intensity':bgf, 'Neurite intensity':nf, 'Neurite signal-to-noise ratio':snr, 'Filtered snr':fsnr}, columns=cols_Data)
    df_Data=df_Data.append(all_data1)

    #find peaks
    peaks = find_peaks(fsnr, height=height, prominence=prom)
    pd = peaks[0]*mu_per_px
    pnd = pd/dist[-1]
    psnr = [fsnr[i] for i in peaks[0]]
    ipd = np.diff(pd)
    ipdd = [pd[i]+ipd[i]/2 for i in range(0,len(ipd))]
    ipdnd = ipdd/dist[-1]
    slope, intercept, r_value, p_value, std_err = stats.linregress(ipdnd,ipd)
    
    #peak width measurements
    #Assumption is that the peak is a triangle.
    #Find the closest minima coordinates on each side of the peak.
    #Find the x-intercepts of the lines joining the peak to the minima on each side.
    #Peak width at half maxima is half the distance between these two points
    lex=[]
    rex=[]
    mins = argrelmin(fsnr, order=3)
    for i in peaks[0]:
        j = mins[0][np.where(mins[0]<i)]
        if len(j)==0:j=0
        else: j=j[-1]
        k = np.where((d<i) & (fsnr==0))[0]
        if len(k)==0: k=0
        else: k=k[-1]
        lex = np.append(lex, max(j,k))
        l = mins[0][np.where(mins[0]>i)]
        if len(l)==0: l=imsize[1]-1
        else: l=l[0]
        m = np.where((d>i) & (fsnr==0))[0]
        if len(m)==0: m=imsize[1]-1
        else: m=m[0]
        rex = np.append(rex, min(l,m))
    ley=[fsnr[int(i)]*0.5 for i in lex]
    rey=[fsnr[int(i)]*0.5 for i in rex]
    lex=[i*mu_per_px for i in lex]
    rex=[i*mu_per_px for i in rex]    
    le=[(lex[i]*fsnr[peaks[0][i]]-pd[i]*ley[i])/(fsnr[peaks[0][i]]-ley[i]) for i in range(len(peaks[0]))]
    re=[(rex[i]*fsnr[peaks[0][i]]-pd[i]*rey[i])/(fsnr[peaks[0][i]]-rey[i]) for i in range(len(peaks[0]))]
    pwhm = [(re[i]-le[i])*0.5 for i in range(len(le))]
    le2 = [pd[i]-0.5*(pd[i]-le[i]) for i in range(len(peaks[0]))]
    re2 = [pd[i]+0.5*(re[i]-pd[i]) for i in range(len(peaks[0]))]
    
    print x
    #create plots
    plt.figure(1, figsize=(0.015*imsize[1],5))
    plt.title(x+' peaks')
    plt.xlabel('Distance (um)')
    plt.ylabel('Normalized intensity')
    plt.axis([0, max(dist), 0, max(snr)])
    plt.plot(dist, snr, 'y-')
    plt.plot(dist, fsnr, 'r-')
    plt.plot(pd, psnr, 'go')
    plt.show()
    plt.close()
    
    #add data to pandas dataframe
    all_data2 = pandas.DataFrame({'Date':[date]*len(pd), 'Strain':[strain]*len(pd), 'Neuron':[neuron]*len(pd), 'L/R':[lr]*len(pd), 'ImageID':[x]*len(pd), 'Distance':pd, 'Normalized distance':pnd, 'Punctum signal-to-noise ratio':psnr, 'Punctum width':pwhm}, columns=cols_Peaks)
    df_Peaks=df_Peaks.append(all_data2)
    all_data3 = pandas.DataFrame({'Date':[date]*len(ipd), 'Strain':[strain]*len(ipd), 'Neuron':[neuron]*len(ipd), 'L/R':[lr]*len(ipd), 'ImageID':[x]*len(ipd), 'Distance':ipdd, 'Normalized distance':ipdnd, 'Inter-punctum interval':ipd}, columns=cols_IPDs)
    df_IPDs=df_IPDs.append(all_data3)    
    frame = pandas.DataFrame([[date, strain, neuron, lr, x, imsize[1], dist[-1], np.mean(nf), np.mean(snr), strength, 'Puncta analyzed', len(pd), np.mean(psnr), np.mean(pwhm), np.mean(ipd), np.median(ipd), len(pd)/dist[-1], slope, intercept, r_value, p_value, std_err]], columns=cols_Analysis)
    df_Analysis = df_Analysis.append(frame)

    plt.figure(2, figsize=(10,5))
    sns.scatterplot(x='Normalized distance', y='Inter-punctum interval', data=all_data3)
    plt.show()
    plt.close()
    
#save data to excel file
dfpath = '/Users/alkadas/Desktop/Data/IPD/'        #filepath where the excel file will be stored
#dfpath = 'C:/Users/alaka/Google Drive/IPD/'       #filepath where the excel file will be stored
timestamp = today.replace('-','')+'-'+now.replace(':','')[:6]
wb = pandas.ExcelWriter(dfpath+timestamp+'_'+strain+'_Data.xlsx', engine='xlsxwriter')
df_Data.to_excel(wb, sheet_name='Data')
df_Peaks.to_excel(wb, sheet_name='Peaks')
df_IPDs.to_excel(wb, sheet_name='IPDs')
df_Analysis.to_excel(wb, sheet_name='Analysis')
df_Parameters.to_excel(wb, sheet_name='Parameters')
wb.save()