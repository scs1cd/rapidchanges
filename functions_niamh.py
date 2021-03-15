import sys
import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from scipy.stats import lognorm
from scipy.stats import norm
mpl.rcParams['agg.path.chunksize'] = 100000
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size': 20})
rc('text', usetex=True)

def define_parameters(variable_no):
    tl = -77.5
    tu = 82.5
    ti = 20
    var= variable_no
    ntl= np.float(tl)
    ntu= np.float(tu)
    return tl, tu, ti, var, ntl, ntu

def define_parameters_all(variable_no):
    tl = -87.5
    tu = 87.5
    ti = 5
    var= variable_no
    ntl= np.float(tl)
    ntu= np.float(tu)
    return tl, tu, ti, var, ntl, ntu

def separate_latitudes_no_time(infile, tl, tu, ti, var, ntl, ntu):
    #This function seperates the model into different latitudes with 20 degree increments
    file   = infile + "_lat=" + str(tl)
    lat0, lon0, dtdt0, dvgpdt0, dvdmdt0, dfdt0  = np.loadtxt(file, usecols=(0,1,2,3,4,5), skiprows=1, unpack='true')
    c = 0 
    l  = len(lat0)

    n  = int(np.ceil(ntu - ntl) / np.abs(ti)+1)
    adtdt  = np.zeros([int(n),l])
    advgpdt= np.zeros([int(n),l])
    advdmdt= np.zeros([int(n),l]) 
    adfdt  = np.zeros([int(n),l])
    name   = np.zeros([int(n)])

    for i in np.arange(ntl,(ntu+1),ti):
        file   = infile + "_lat=" + str((i))
        # Don't need the int() for GGF
        print(i, file)
        try:
            adtdt[c,:],advgpdt[c,:],advdmdt[c,:],adfdt[c,:] = np.loadtxt(file,usecols=(2,3,4,5),skiprows=1,unpack='true')
        except StopIteration:
            print("\033[31m" + file + " is empty!\033[0m")
        name[c] = i
        c = c + 1
    return name, adtdt,advgpdt,advdmdt,adfdt,n

def separate_latitudes(infile, tl, tu, ti, var, ntl, ntu):
    #This function seperates the model into different latitudes with 20 degree increments
    file   = infile + "_lat=" + str(tl)
    lat0, lon0, dtdt0, dvgpdt0, dvdmdt0, dfdt0, time  = np.loadtxt(file, usecols=(0,1,2,3,4,5,6), skiprows=1, unpack='true')
    c = 0 
    l  = len(lat0)

    n  = int(np.ceil(ntu - ntl) / np.abs(ti)+1)
    adtdt  = np.zeros([int(n),l])
    advgpdt= np.zeros([int(n),l])
    advdmdt= np.zeros([int(n),l]) 
    adfdt  = np.zeros([int(n),l])
    atime  = np.zeros([int(n),l])
    name   = np.zeros([int(n)])

    for i in np.arange(0,(ntu+1),ti):
        file   = infile + "_lat=" + str((i))
        
        print(i, file)
        try:
            adtdt[c,:],advgpdt[c,:],advdmdt[c,:],adfdt[c,:],atime[c,:] = np.loadtxt(file,usecols=(2,3,4,5,6),skiprows=1,unpack='true')
        except StopIteration:
            print("\033[31m" + file + " is empty!\033[0m")
        name[c] = i
        c = c + 1
    return name, adtdt,advgpdt,advdmdt,adfdt,atime,n

def separate_latitudes_all(infile, tl, tu, ti, var, ntl, ntu, line_number):
    #This function seperates the model into different latitudes with 20 degree increments
    file   = infile + "_lat=2.5"
    lat0, lon0, dtdt0, dvgpdt0, dvdmdt0, dfdt0, time  = np.loadtxt(file, usecols=(0,1,2,3,4,5,6), skiprows=1, unpack='true')
    c = 0 
    n  = int(np.ceil(ntu - ntl) / np.abs(ti))+1
    print(n)
    l  = len(lat0)
    # Reread the correct file
    lat0, lon0, dtdt0, dvgpdt0, dvdmdt0, dfdt0, time  = np.loadtxt(infile, usecols=(0,1,2,3,4,5,6), skiprows=1, unpack='true')
    
    alat   = np.zeros([int(n),l])
    alon   = np.zeros([int(n),l])
    adtdt  = np.zeros([int(n),l])
    advgpdt= np.zeros([int(n),l])
    advdmdt= np.zeros([int(n),l]) 
    adfdt  = np.zeros([int(n),l])
    atime  = np.zeros([int(n),l])
    name   = np.zeros([int(n)])
    
    start = 0

    for i in np.arange(1,n+1,1):
        
        end = int(((len(lat0)/n) -1)* i)
        #print(end)
        alat[c,:]    = lat0[start:end]
        alon[c,:]    = lon0[start:end]
        adtdt[c,:]   = dtdt0[start:end]
        advgpdt[c,:] = dvgpdt0[start:end]
        advdmdt[c,:] = dvdmdt0[start:end]
        adfdt[c,:]   = dfdt0[start:end]
        atime[c,:]   = time[start:end]
        
        start = end
        
        name[c] = tl + (5 *(i-1))
        c = c + 1
    return name,alat,alon,adtdt,advgpdt,advdmdt,adfdt,atime,n

def P99(adtdt,advgpdt,advdmdt,adfdt,atime,var,n):
    alltop1 = []
    if var==1:
        Q = adtdt
    elif var==2:
        Q = advgpdt
    elif var==3:
        Q = abs(advdmdt/1000.0)
    elif var==4: 
        Q = abs(adfdt/1000.0)
    array = Q[Q>0]
    pts = np.sort(array)
        
    muX       = np.mean(np.log(array))
    sigmaX    = np.std(np.log(array)) 
    logNdist1 = lognorm(sigmaX, scale=np.exp(muX), loc=0)
        
    cdf       = logNdist1.cdf(pts)
    fig, (ax1) = plt.subplots(1, figsize=(16,4))
    ncdf, bins, patches = ax1.hist(pts, bins=500, log=False, density=True, histtype='step', cumulative=True)
    
    cdfindex  = np.argmax(ncdf >= 0.99)
    ax1.axhline(y=ncdf[cdfindex],color='r', linestyle='-', linewidth=0.5)
  
    P99index = cdfindex
    P99value = bins[P99index]
    ax1.axvline(x=bins[P99index],color='r', linestyle='-',linewidth=0.5)
    print(P99value)
    for i in range(0,n,1):
        if var==1:
            Q = adtdt[i,:]
            label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
      
        #Can chang (Q>x) to any value for example 0.2 -> double the modern field for directional changes
        indices = np.where(Q>P99value)
        #print(indices)
       #print(np.max(Q))
        keep      = np.zeros(len(indices[0]))

        for j in range(0,len(indices[0]),1):
            keep[j]  = Q[indices[0][j]]
        alltop1.append(keep)
    return alltop1, P99value

def P99time(alon,adtdt,advgpdt,advdmdt,adfdt,atime,var,n):
   
    alltop1 = []
    alltime1 = []
    alllongitudes1 = []
    if var==1:
        Q = adtdt
    elif var==2:
        Q = advgpdt
    elif var==3:
        Q = abs(advdmdt/1000.0)
    elif var==4: 
        Q = abs(adfdt/1000.0)
    array = Q[Q>0]
    pts = np.sort(array)
        
    muX       = np.mean(np.log(array))
    sigmaX    = np.std(np.log(array)) 
    logNdist1 = lognorm(sigmaX, scale=np.exp(muX), loc=0)
        
    cdf       = logNdist1.cdf(pts)
    fig, (ax1) = plt.subplots(1, figsize=(16,4))
    ncdf, bins, patches = ax1.hist(pts, bins=500, log=False, density=True, histtype='step', cumulative=True)
    cdfindex  = np.argmax(ncdf >= 0.99)
        
    P99index = cdfindex
    P99value = bins[P99index]
    #print(P99value)
    for i in range(0,n,1):
        t = atime[i,:]
        lon = alon[i,:]
        if var==1:
            Q = adtdt[i,:]
            #label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            #label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            #label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            #label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
      
        #Can chang (Q>x) to any value for example 0.2 -> double the modern field for directional changes
        indices = np.where(Q>P99value) # THIS IS ORIGINALLY indices = np.where(Q>P99value)
        #print(indices)
       # print(np.max(Q))
        keep      = np.zeros(len(indices[0]))
        times     = np.zeros(len(indices[0]))
        longitudes= np.zeros(len(indices[0]))
        for j in range(0,len(indices[0]),1):
            keep[j]  = Q[indices[0][j]]
            times[j] = t[indices[0][j]]
            longitudes[j] = lon[indices[0][j]]
        alltop1.append(keep)
        alltime1.append(times)
        alllongitudes1.append(longitudes)
       
    return alltop1, alltime1, alllongitudes1

def events_x4(alon,adtdt,advgpdt,advdmdt,adfdt,atime,var,n):
   #Finds all events that a faster than 4 times the rate of change of the modern field (0.1)
    alltop4 = []
    alltime4 = []
    alllongitudes4 = []
    
    for i in range(0,n,1):
        t = atime[i,:]
        lon = alon[i,:]
        if var==1:
            Q = adtdt[i,:]
            #label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            #label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            #label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            #label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
      
        #Can chang (Q>x) to any value for example 0.2 -> double the modern field for directional changes
        indices = np.where(Q>0.4) 
        #print(indices)
       # print(np.max(Q))
        keep      = np.zeros(len(indices[0]))
        times     = np.zeros(len(indices[0]))
        longitudes= np.zeros(len(indices[0]))
        for j in range(0,len(indices[0]),1):
            keep[j]  = Q[indices[0][j]]
            times[j] = t[indices[0][j]]
            longitudes[j] = lon[indices[0][j]]
        alltop4.append(keep)
        alltime4.append(times)
        alllongitudes4.append(longitudes)
       
    return alltop4, alltime4, alllongitudes4

def calctime(starttime,line_number,year_intervals):
    # Calculate the time of events within the models (for top 1%)
    calctop1  = []
    calctime1 = []
    for i in range(0,n,1):
        #t = atime[i,:]
        if var==1:
            Q     = adtdt[i,:]
            #label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            #label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            #label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            #label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
      
        array = Q[Q>0]
        pts   = np.sort(array)
        
        muX       = np.mean(np.log(array))
        sigmaX    = np.std(np.log(array)) 
        logNdist1 = lognorm(sigmaX, scale=np.exp(muX), loc=0)
        
        #cdf       = logNdist1.cdf(pts)
        fig, (ax1) = plt.subplots(1, figsize=(16,4))
        ncdf, bins, patches = ax1.hist(pts, bins=500, log=False, density=True, histtype='step', cumulative=True)
        cdfindex = np.argmax(ncdf >= 0.99)
        #print(cdfindex)
        P99index = cdfindex
        P99value = pts[P99index]
        print(P99value)
        indices  = np.where(Q>P99value)
        #print(indices)
        #print(np.max(Q))
        keep   = np.zeros(len(indices[0]))
        times  = np.zeros(len(indices[0]))
        for j in range(0,len(indices[0]),1):
            for i in range(0, int(line_number), 1):
                if i<=indices[0][j]/int(line_number)<(i+1):
                    keep[j]  = Q[indices[0][j]]
                    times[j] = starttime +(indices[0][j] - (i*int(line_number))+1)*year_intervals
                    #print(times[j])
        calctop1.append(keep)
        calctime1.append(times)
    return calctop1,calctime1



def plot_figures(adtdt,advgpdt,advdmdt,adfdt,var,n, name):
    # Function plots CDF and PDF histograms
    plt.clf()
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    muXlat = np.zeros([int(n)])
    maxlat = np.zeros([int(n)])
    sigma = np.zeros([int(n)])
    for i in range(2,n,1):
        if var==1:
            Q = adtdt[i,:]
            #label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            #label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            #label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            #label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
      
        if i == 2: c = "yellow"
        if i == 3: c = "orange"
        if i == 4: c = "red"
        if i == 5: c = "green"
        if i == 6: c = "blue"
        if i == 7: c = "purple"
        if i == 8: c = "black"
     
        array     = Q[Q>0]
        muX       = np.mean(np.log(array))
        sigmaX    = np.std(np.log(array)) 
        sigmaraw  = np.std(array)
        logNdist1 = lognorm(sigmaX, scale=np.exp(muX), loc=0)
        print('theta = ', name[i], ' mean =', np.exp(muX), ' sigma = ', sigmaX)
        muXlat[i] = np.exp(muX)
        maxlat[i] = np.max(array)
        sigma[i] = sigmaraw * np.exp(muX)

        pts = np.sort(array)
        print(np.max(pts))
        logNdist1 = lognorm(sigmaX  , scale=np.exp(muX), loc=0)
       #logNdist2 = lognorm(stats[0], scale=stats[2], loc=stats[1])
       #logdist   = norm(scale=sigmaX, loc=muX)

        tmp = logNdist1.cdf(pts)
        print('cdf = ', tmp[-1])
        
        # This is for latitudes 2.5 22.5 42.5
        ax1.plot(pts, logNdist1.cdf(pts), color=c, linestyle='-', label=name[i])
        if  i ==4 or i==5 or i==6 : 
            ax2.plot(pts, logNdist1.pdf(pts), color=c, linestyle='-', label=name[i])
            ax2.hist(pts, bins=1000, log=True, color=c, alpha=0.6, density=True)
    for i in range(0,2,1):
        if var==1:
            Q = adtdt[i,:]
            #label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            #label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            #label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            #label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
        array     = Q[Q>0]
        muX       = np.mean(np.log(array))
        sigmaX    = np.std(np.log(array))
        sigmaraw  = np.std(array)
        muXlat[i] = np.exp(muX)
        maxlat[i] = np.max(array)
        sigma[i]  = sigmaraw * np.exp(muX)
    print(maxlat)
    #plt.plot(pts, 0.95*np.ones(len(pts)))
    ax1.legend(loc='best', borderaxespad=0)
    ax1.set_xlim([1e-3,2e0])
    ax1.set_xscale('log')
    #ax1.set_xlabel(label)
    ax1.set_ylabel("cdf")
    fig1.savefig("cdflat_var="+str(var)+".pdf", format='pdf', bbox_inches='tight',pad_inches = 0.1)

    ax2.legend(loc='best', borderaxespad=0)
    ax2.set_ylim([1e-5,1e2])
    #ax2.set_xlim([-0.01,1e1])
    ax2.set_yscale('log')
    #ax2.set_xlabel(label)
    ax2.set_ylabel("pdf")
    fig2.savefig("pdfall_var="+str(var)+".pdf", format='pdf', bbox_inches='tight',pad_inches = 0.1)
    return fig1, fig2, muXlat, maxlat, sigma

def mean_max(adtdt,advgpdt,advdmdt,adfdt,var,n, name):
    muXlat = np.zeros([int(n)])
    maxlat = np.zeros([int(n)])
    sigma = np.zeros([int(n)])
    for i in range(0,n,1):
        if var==1:
            Q = adtdt[i,:]
            #label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            #label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            #label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            #label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
     
        array     = Q[Q>0]
        muX       = np.mean(np.log(array))
        sigmaX    = np.std(np.log(array)) 
        sigmaraw  = np.std(array)
        muXlat[i] = np.exp(muX)
        maxlat[i] = np.max(array)
        sigma[i] = sigmaraw * np.exp(muX)

    return muXlat, maxlat, sigma

    
def no_excursionGGF(atime,adtdt,advgpdt,advdmdt,adfdt,var,n, name):
    #Find maximum point for GGF when removing the excusrions
    max_no_excursion = np.zeros([int(n)])
    
    for i in range(0,n,1):
        if var==1:
            Q = adtdt[i,:]
            time = atime[i,:]
            #label = "$\partial \\hat{\\mathbf{B}}/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==2:
            Q = advgpdt[i,:]
            time = atime[i,:]
            #label = "$\partial \\hat{\\mathbf{P}}_V/ \partial t$ ($^{\\circ}$ yr$^{-1}$)"
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
            time = atime[i,:]
            #label = "$\partial P_V/ \partial t$ ($ZAm^2$ yr$^{-1}$)"
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            time = atime[i,:]
            #label = "$\partial B/ \partial t$ (${\\mu}$T yr$^{-1}$)"
        array     = Q[Q>0]
        pts = np.sort(array)
        
        indices = [np.where((time > 42450) & (42450<time<37000) & (37000<time<35000) & (time < 30000))]
        print(indices)
        values = np.zeros(len(indices))
        values = pts[indices]
        max_no_excursion[i] = np.max(no_excursion)           
    return max_no_excursion

                           
def longitudes(infile, tl, tu, ti, var, ntl, ntu, name):
    #This function seperates the model into different latitudes with 20 degree increments
    file   = infile + "_lat=2.5" 
    lat0, lon0, dtdt0, dvgpdt0, dvdmdt0, dfdt0, time  = np.loadtxt(file, usecols=(0,1,2,3,4,5,6), skiprows=1, unpack='true')
    c = 0 
    l  = len(lat0)

    n  = int(np.ceil(ntu - ntl) / np.abs(ti)+1)
    alon   = np.zeros([int(n),l])
    
    lat0, lon0, dtdt0, dvgpdt0, dvdmdt0, dfdt0, time  = np.loadtxt(infile, usecols=(0,1,2,3,4,5,6), skiprows=1, unpack='true')

    for i in np.arange(ntl,(ntu+1),ti):
        
        # Don't need the int() for GGF
        print(i, file)
        try:
            alon[c,:] = np.loadtxt(file,usecols=(1),skiprows=1,unpack='true')
        except StopIteration:
            print("\033[31m" + file + " is empty!\033[0m")
        name[c] = i
        c = c + 1
    return alon

def empirical_cdf(adtdt,advgpdt,advdmdt,adfdt,var,n):
    if var==1:
        Q = adtdt
    elif var==2:
        Q = advgpdt
    elif var==3:
        Q = abs(advdmdt/1000.0)
    elif var==4: 
        Q = abs(adfdt/1000.0)
    array = Q[Q>0]
    pts = np.sort(array)
        
    muX       = np.mean(np.log(array))
    sigmaX    = np.std(np.log(array)) 
    logNdist1 = lognorm(sigmaX, scale=np.exp(muX), loc=0)
        
    cdf       = logNdist1.cdf(pts)
    fig, (ax1) = plt.subplots(1, figsize=(16,4))
    ncdf, bins, patches = ax1.hist(pts, bins=500, log=False, density=True, histtype='step', cumulative=True)
    cdfindex  = np.argmax(ncdf >= 0.99)
        
    P99index = cdfindex
    P99value = bins[P99index]
    print(P99value)
    return P99value

def IQR(adtdt,advgpdt,advdmdt,adfdt,var,n,name):
    
    muXlat = np.zeros([int(n)])
    maxlat = np.zeros([int(n)])
    IQR    = np.zeros([int(n)])
    Q1     = np.zeros([int(n)])
    Q3     = np.zeros([int(n)])
    for i in range(2,n,1):
        if var==1:
            Q = adtdt[i,:]
        elif var==2:
            Q = advgpdt[i,:]
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
            
        array     = Q[Q>0]
        muX       = np.mean(np.log(array))
        sigmaX    = np.std(np.log(array))   
        logNdist1 = lognorm(sigmaX, scale=np.exp(muX), loc=0)
        print('theta = ', name[i], ' mean =', np.exp(muX), ' sigma = ', sigmaX)
        muXlat[i] = np.exp(muX)
        maxlat[i] = np.max(array)

        pts = np.sort(array)
        # First quartile
        Q1[i] = np.percentile(pts, 25)
        # Third quartile
        Q3[i] = np.percentile(pts, 75)
        # Interquartile range
        IQR[i] = Q3[i] - Q1[i]

        
    for i in range(0,2,1):
        if var==1:
            Q = adtdt[i,:]
        elif var==2:
            Q = advgpdt[i,:]
        elif var==3:
            Q = abs(advdmdt[i,:]/1000.0)
        elif var==4: 
            Q = abs(adfdt[i,:]/1000.0) # microTesla
        array     = Q[Q>0]
        muX       = np.mean(np.log(array))
        muXlat[i] = np.exp(muX)
        maxlat[i] = np.max(array)
        pts = np.sort(array)
        # First quartile
        Q1[i] = np.percentile(pts, 25)
        # Third quartile
        Q3[i] = np.percentile(pts, 75)
        # Interquartile range
        IQR[i] = Q3[i] - Q1[i]    
    return muXlat, maxlat, IQR, Q1, Q3