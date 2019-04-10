#import basic libraries for plotting, data structures and signal processing
import matplotlib.pyplot as plt
import numpy as np
import pandas
from scipy.signal import find_peaks, argrelmin,argrelmax
from scipy.ndimage import gaussian_filter
from scipy import stats
import imageio
import os
import fnmatch
import seaborn as sns

fpath = '/Users/alkadas/Desktop/Data/Microscopy/IPD/GN878/'        #filepath where the data is
#fpath = 'C:/Users/alaka/Google Drive/Microscopy/IPD/GN877/'
imgfiles = fnmatch.filter(os.listdir(fpath), '*.tif')

mu_per_px = 0.18825     #pixels to microns conversion factor


#specify columns of the pandas dataframe and excel sheets
cols_Data = ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance','Raw intensity', 'Background intensity', 'Neurite intensity', 'Neurite signal-to-noise ratio', 'Filtered snr']
cols_Peaks = ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance', 'Punctum signal-to-noise ratio', 'Punctum width']
cols_IPDs = ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Distance', 'Normalized distance', 'Inter-punctum interval']
cols_Metadata = ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Image size', 'Max neurite length', 'Average neurite intensity', 'Average snr', 'Signal strength','Status','Total peaks', 'Average peak snr', 'Average peak width', 'Average ipd', 'Median ipd', 'Average peaks per distance', 'slope', 'intercept', 'r_value', 'p_value', 'std_err']
cols_RejectsMetadata = ['Date', 'Strain', 'Neuron', 'L/R', 'ImageID', 'Image size', 'Max neurite length', 'Average neurite intensity', 'Average snr', 'Signal strength', 'Status']

#initialize Pandas DataFrames
df_Data = pandas.DataFrame()
df_Peaks = pandas.DataFrame()
df_IPDs = pandas.DataFrame()
df_Metadata = pandas.DataFrame()
df_RejectsMetadata = pandas.DataFrame()

for x in imgfiles:                          #create loop for number of images in folder
    img = imageio.imread(fpath+x)[:,:,1]           #import image and store it in a list of lists

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
        if nf[i]<0: nf[i]=0
    snr = nf/bgf                            #calculate signal to noise ratio

    strength = np.mean([snr[i] for i in argrelmax(snr,order=10)[0]])-np.mean([snr[i] for i in argrelmin(snr,order=10)[0]])
    gfsnr=gaussian_filter(snr, 1/strength)

    #add image data to pandas dataframe
    all_data1 = pandas.DataFrame({'Date':[date]*imsize[1], 'Strain':[strain]*imsize[1], 'Neuron':[neuron]*imsize[1], 'L/R':[lr]*imsize[1], 'ImageID':[x]*imsize[1], 'Distance':dist, 'Normalized distance':normdist, 'Raw intensity':rawf, 'Background intensity':bgf, 'Neurite intensity':nf, 'Neurite signal-to-noise ratio':snr, 'Filtered snr':gfsnr}, columns=cols_Data)
    df_Data=df_Data.append(all_data1)

    #rejecting bad quality images
    cutoff=0.6
    if strength<cutoff:
        frame = pandas.DataFrame([[date, strain, neuron, lr, x, imsize[1], dist[-1], np.mean(nf), np.mean(snr), strength, 'Puncta not analyzed']], columns=cols_RejectsMetadata)
        df_RejectsMetadata = df_RejectsMetadata.append(frame)
        continue


    #find peaks
    peaks = find_peaks(gfsnr, height=0.2*strength, prominence=0.1*strength)
    pd = peaks[0]*mu_per_px
    pnd = pd/dist[-1]
    psnr = [gfsnr[i] for i in peaks[0]]
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
    mins = argrelmin(gfsnr, order=3)
    for i in peaks[0]:
        j = mins[0][np.where(mins[0]<i)]
        if len(j)==0:j=0
        else: j=j[-1]
        k = np.where((d<i) & (gfsnr==0))[0]
        if len(k)==0: k=0
        else: k=k[-1]
        lex = np.append(lex, max(j,k))
        l = mins[0][np.where(mins[0]>i)]
        if len(l)==0: l=imsize[1]-1
        else: l=l[0]
        m = np.where((d>i) & (gfsnr==0))[0]
        if len(m)==0: m=imsize[1]-1
        else: m=m[0]
        rex = np.append(rex, min(l,m))
    ley=[gfsnr[int(i)]*0.5 for i in lex]
    rey=[gfsnr[int(i)]*0.5 for i in rex]
    lex=[i*mu_per_px for i in lex]
    rex=[i*mu_per_px for i in rex]    
    le=[(lex[i]*gfsnr[peaks[0][i]]-pd[i]*ley[i])/(gfsnr[peaks[0][i]]-ley[i]) for i in range(len(peaks[0]))]
    re=[(rex[i]*gfsnr[peaks[0][i]]-pd[i]*rey[i])/(gfsnr[peaks[0][i]]-rey[i]) for i in range(len(peaks[0]))]
    pwhm = [(re[i]-le[i])*0.5 for i in range(len(le))]
    le2 = [pd[i]-0.5*(pd[i]-le[i]) for i in range(len(peaks[0]))]
    re2 = [pd[i]+0.5*(re[i]-pd[i]) for i in range(len(peaks[0]))]
    
    print x
    #plot raw and filtered signal-to-noise ratio vs. distance from cell body and identified peaks
    plt.figure(1, figsize=(0.015*imsize[1],5))
    plt.title(x+' peaks')
    plt.xlabel('Distance (um)')
    plt.ylabel('Normalized intensity')
    plt.axis([0, max(dist), 0, max(snr)])
    plt.plot(dist, snr, 'y-')
    plt.plot(dist, gfsnr, 'r-')
    plt.plot(pd, psnr, 'go')
    plt.show()
    plt.close()
    
    
    #add data to pandas dataframe
    all_data2 = pandas.DataFrame({'Date':[date]*len(pd), 'Strain':[strain]*len(pd), 'Neuron':[neuron]*len(pd), 'L/R':[lr]*len(pd), 'ImageID':[x]*len(pd), 'Distance':pd, 'Normalized distance':pnd, 'Punctum signal-to-noise ratio':psnr, 'Punctum width':pwhm}, columns=cols_Peaks)
    df_Peaks=df_Peaks.append(all_data2)
    all_data3 = pandas.DataFrame({'Date':[date]*len(ipd), 'Strain':[strain]*len(ipd), 'Neuron':[neuron]*len(ipd), 'L/R':[lr]*len(ipd), 'ImageID':[x]*len(ipd), 'Distance':ipdd, 'Normalized distance':ipdnd, 'Inter-punctum interval':ipd}, columns=cols_IPDs)
    df_IPDs=df_IPDs.append(all_data3)    
    frame = pandas.DataFrame([[date, strain, neuron, lr, x, imsize[1], dist[-1], np.mean(nf), np.mean(snr), strength, 'Puncta analyzed', len(pd), np.mean(psnr), np.mean(pwhm), np.mean(ipd), np.median(ipd), len(pd)/dist[-1], slope, intercept, r_value, p_value, std_err]], columns=cols_Metadata)
    df_Metadata = df_Metadata.append(frame)

    #plot inter-punctum interval vs. distance from cell body
    plt.figure(2, figsize=(10,5))
    sns.scatterplot(x='Normalized distance', y='Inter-punctum interval', data=all_data3)
    plt.show()
    plt.close()


#save data to excel file
dfpath = '/Users/alkadas/Desktop/Data/Microscopy/IPD/'        #filepath where the excel file will be saved
#dfpath = 'C:/Users/alaka/Google Drive/Microscopy/IPD/'
wb = pandas.ExcelWriter(dfpath+strain+'_Data.xlsx', engine='xlsxwriter')
df_Data.to_excel(wb, sheet_name='Data')
df_Peaks.to_excel(wb, sheet_name='Peaks')
df_IPDs.to_excel(wb, sheet_name='IPDs')
df_Metadata.to_excel(wb, sheet_name='Metadata')
df_RejectsMetadata.to_excel(wb, sheet_name='RejectsMetadata')
wb.save()
