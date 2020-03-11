import random
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
for y in range(10):
    l = random.randint(750, 1000)
    x = np.arange(0,l)
    #Baseline noise and curvature
    bn = [random.randint(0,5) for m in x]
    #bcx=random.randint(25, 50)
    #tau=random.randint(100, 200)
    bc1 = random.randint(25, 50)*signal.exponential(l, center=0, tau=random.randint(100, 200), sym=False)
    bc2 = random.randint(50,100)+random.randint(0,50)*np.sin(0.01*x+random.randint(0,50))
    #Define number of low intensity background peaks
    bspn = l/3
    #Randomize peak height, peak position
    bsp = [0]*l
    for m in range(bspn):
        bsptemp = random.randint(0, 10)*signal.exponential(l, center=random.randint(0, l), tau=random.randint(1,2), sym=False)
        bsp=bsp+bsptemp
        #Define number of prominent puncta
    pn = random.randint(l/25, l/20)
    #Randomize peak height and peak position
    gpeaks = [0]*l
    for m in range(pn):
        gptemp=random.randint(40, 80)*signal.gaussian(10, random.randint(1, 2), sym=True)
        pp=random.randint(0, l-10)
        for n in range(pp, pp+10):
            gpeaks[n]=gpeaks[n]+gptemp[n-pp]
    total=bn+bsp+bc1+bc2+gpeaks
    bs=bn+bc2
    p=bc1+bsp+gpeaks
    plt.figure(1, figsize=(20,5))
    plt.axis([0, l, 0, 255])
    plt.plot(x,bs,'b-')
    plt.plot(x,total,'y-')
    plt.plot(x,p,'r-')
    plt.show()
    plt.close()