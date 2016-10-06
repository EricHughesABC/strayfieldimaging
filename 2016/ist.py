# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:32:18 2016

@author: flph74
"""
import numpy
from scipy import fftpack
from   scipy.fftpack import fft, ifft, hilbert, fftshift



def numPeaks( sss ):
    www = numpy.where( sss > 0 )
    numPeaks = 0
    ppp = www[0]
    for i in range(len(ppp[:-1])):
        v = www[0][i]
        v1 = www[0][i+1]
        if v1-v > 1:
            numPeaks +=1

    return( numPeaks ) 



def compressed_sensing_1VD( fid, mask, num_iter=500, factor=0.95, tol = 0.01, maxPeaks=2 ):
    
    sss = numpy.zeros( len(fid))
    sss0 = numpy.zeros( len(fid))
    
    final_iter_value  = 0
    final_tol_value   = 0
    final_numpeaks_value = 0

    fid1 = fid.copy()
    tol_diff = (numpy.sqrt(((abs(fid1)).sum())/32.))
    
    k=0
    kkk = []
    rrr = fftpack.fft(fid1)
    rss0 = []
    rss0.append(abs(rrr).sum())
    
    tol0 = abs(rrr).sum()
    tol_diff = ( tol0 - abs(rrr).sum() )*100.0 / tol0
    while (tol_diff < tol) and (k < num_iter) and numPeaks( sss ) <maxPeaks:
        
        sss0 = 1.0*sss
        
        rrr = fftpack.fft(fid1)
        m1 = max(rrr.real)
        
        sss_max_old = sss.max()
               
        for i,r in enumerate(rrr):
            if r.real > m1* factor:
                sss[i] = sss[i]+rrr[i].real
                rrr[i] = complex(m1*factor)
        sss_max = sss.max()
        
        rrr_iii = fftpack.hilbert( rrr.real )
        
        rrr = rrr.real + 1j * rrr_iii
        
        fid1 = fftpack.ifft(rrr)
        
        fid1 *= mask
        tol_diff = ( tol0 - abs(rrr).sum() )*100.0 / tol0
        k +=1

    final_iter_value = k
    final_numpeaks_value = numPeaks(sss)
    final_tol_value = tol_diff
    
    return( sss0, [final_iter_value, final_numpeaks_value, final_tol_value ] )
    
    
    
def compressed_sensing_1VD1b( fid, mask, num_iter=500, factor=0.95, maxPeaks=20 ):
    
#    print maxPeaks,num_iter
    sss = numpy.zeros( len(fid))
    sss_prev = numpy.zeros( len(fid))
    
    final_iter_value  = 0
    final_tol_value   = 0
    final_numpeaks_value = 0
    
#    numberOfpeaks = numpy.zeros(num_iter)
    numberOfpeaks = []
    fid1 = fid.copy()
    
    k=0
    rrr = fft(fid1)
    
#    print "Num Peaks", numPeaks( sss )
    
    noise_not_reached = True
#    while  (k < num_iter) and numPeaks( sss[k] ) < maxPeaks and noise_not_reached:

        
    while  noise_not_reached and numPeaks( sss ) <maxPeaks:
        
        sss_prev =  numpy.copy(sss)
        
        rrr = fft(fid1)
        m1 = max(rrr.real)*factor
               
        for i,r in enumerate(rrr):
            if r.real > m1:
                sss[i] = sss[i]+rrr[i].real-m1
                rrr[i] = complex(m1)
        rrr = rrr.real + 1j * hilbert( rrr.real )
        
        fid1 = ifft(rrr)*mask

        
        numberOfpeaks.append( numPeaks( sss ))

#        k +=1
#        print "k, num peaks", k,numberOfpeaks[-1]
        
        npeaks = numpy.array(numberOfpeaks)
    
#        ddd = abs(npeaks-numpy.roll(npeaks,-1))

#        ddd1 = ddd[:]-numpy.roll(ddd,-1)
#        ppp = (numpy.where(ddd[:-2]>0))[0]

#        ppp1 = ppp[1:]-ppp[:-1]
        
        noise_not_reached = True
#        print "ppp,ppp1",ppp,ppp1
        if k>=3:
            ddd = abs(npeaks[1:]-npeaks[:-1])
            ppp = (numpy.where(ddd>0))[0]
            ppp1 = ppp[1:]-ppp[:-1]
            for j,v in enumerate(ppp1):


                if v < 1000:
                    noise_not_reached = False 
                    print "reached noise level",v,npeaks[-1],k
                    break

        k += 1
        if k >  num_iter:
            print "reached max iterations"
            noise_not_reached = False
        
    final_iter_value = k

#    print k,
                    
    return( sss_prev, [final_iter_value, numberOfpeaks[-1] ] )





def compressed_sensing_1VD1c( fid, mask, num_iter=500, factor=0.95, maxPeaks=20, peak_separation = 100 ):
    
    sss = numpy.zeros( len(fid))
    sss_prev = numpy.zeros( len(fid))
    
    final_iter_value  = 0
    final_tol_value   = 0
    final_numpeaks_value = 0
    
    numberOfpeaks = []
    fid1 = fid.copy()
    
    k=0
    rrr = fft(fid1)
    
    
    noise_not_reached = True

        
    while  noise_not_reached and numPeaks( sss ) <maxPeaks:
        
        sss_prev =  numpy.copy(sss)
        
        rrr = fft(fid1)
        m1 = max(rrr.real)*factor
               
        for i,r in enumerate(rrr):
            if r.real > m1:
                sss[i] = sss[i]+rrr[i].real-m1
                rrr[i] = complex(m1)
        rrr = rrr.real + 1j * hilbert( rrr.real )
        
        fid1 = ifft(rrr)*mask

        
        numberOfpeaks.append( numPeaks( sss ))

        
        npeaks = numpy.array(numberOfpeaks)
    
#        noise_not_reached = True
        if k>=5:
            ddd = abs(npeaks[1:]-npeaks[:-1])
            ppp = (numpy.where(ddd>0))[0]
            ppp1 = ppp[1:]-ppp[:-1]
            for j,v in enumerate(ppp1):


                if v < peak_separation and npeaks[-1]>5:
                    noise_not_reached = False 
                    print "reached noise level",v,npeaks[-1],k
                    break

        k += 1
        if k >  num_iter:
            print "reached max iterations"
            noise_not_reached = False
        
    final_iter_value = k

#    print k,
                    
    return( sss_prev, [final_iter_value, numberOfpeaks[-1] ] )




def compressed_sensing_1VD1d( fid, mask, num_iter=500, factor=0.95, maxPeaks=20, peak_separation = 100, tolerance=1e-3 ):
    
    sss = numpy.zeros( len(fid))
    sss_prev = numpy.zeros( len(fid))
    
    final_iter_value  = 0
    final_tol_value   = 0
    final_numpeaks_value = 0
    
    numberOfpeaks = []
    sig_difference = [0,]
    fid1 = fid.copy()
    
    k=0
    rrr = fft(fid1)
    
    
    noise_not_reached = True

        
    while  noise_not_reached and numPeaks( sss ) <maxPeaks:
        
        sss_prev =  numpy.copy(sss)
        
        rrr = fft(fid1)
        m1 = max(rrr.real)*factor
               
        for i,r in enumerate(rrr):
            if r.real > m1:
                sss[i] = sss[i]+rrr[i].real-m1
                rrr[i] = complex(m1)
        sig_difference.append( float(sss.sum()))
        rrr = rrr.real + 1j * hilbert( rrr.real )
        
        fid1 = ifft(rrr)*mask

        
        numberOfpeaks.append( numPeaks( sss ))

        
        npeaks = numpy.array(numberOfpeaks)
    
#        noise_not_reached = True
        if k>=5:
            ddd = abs(npeaks[1:]-npeaks[:-1])
            ppp = (numpy.where(ddd>0))[0]
            ppp1 = ppp[1:]-ppp[:-1]
            for j,v in enumerate(ppp1):


                if v < peak_separation and npeaks[-1]>maxPeaks:
                    noise_not_reached = False 
                    print "reached noise level",v,npeaks[-1],k
                    break

        k += 1
        if k >  num_iter:
            print "reached max iterations"
            noise_not_reached = False
            
        if k>2:
            
#            tol = (float(sig_difference[-1])-float(sig_difference[-2]))/float(sig_difference[-1])
            
            tol = (float(sig_difference[-1])-float(sig_difference[-2]))/float(sig_difference[-1])
            
            if tol < tolerance:
                
                noise_not_reached = False
            
        
    final_iter_value = k

#    print k,
                    
    return( sss_prev, [final_iter_value, numPeaks( sss ), tol ] )