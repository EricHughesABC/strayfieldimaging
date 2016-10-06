# -*- coding: utf-8 -*-
"""
Created on Mon May 09 09:37:35 2016

@author: flph74
"""
import pylab
import scipy
import os
import nmrglue


import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter



class TwoD_NMR_MAT_plot:
    
    def __init__(self, exp, pinfo, info, axis_domains = ["w","w"], dimensions_ppm=[], contourlevels_range=[0.01,1.1] ):
        
#        print exp.shape
        
        self.exp = exp
        self.pinfo = pinfo
        self.info  = info
        self.dimensions = dimensions_ppm
        self.contourlevels_range = contourlevels_range
        self.axis_domains = axis_domains # possible values [w,w] or [t, t] or [w,t] or [t,w]
        
        self.plot_size = [6,6]
        
        self.rr,self.cc = exp.shape
        self.X = np.zeros(exp.shape)
        self.Y = np.zeros(exp.shape)
        
#        r1=0
        r2=self.rr

#        c1=0
        c2=self.cc
        
#        print r2,c2
        
        self.nucleus_labels = {"1H":"$^1$H ", "2H":"$^2$H ", "13C":"$^{13}$C ", "87Rb":"$^{87}$Rb "}
        
#        self.indirect_nucleus = ""
#        self.observe_nucleus = ""

        
    def display_plot(self):
        
        self.create_axes(  self.pinfo, self.info, self.rr, self.cc, self.dimensions, self.axis_domains )
          
        self.create_plot_layout(self.dimensions_index, self.axis_domains)
        
        self.plot_plots()
        

        
    def create_axes( self, pinfo, info, rr,cc, dimensions_ppm, axis_domains):
                
        self.f1_offset_p = pinfo['procs' ]['OFFSET']
        self.f1_sw_hz     = pinfo['procs' ]['SW_p']
        self.f1_omega    = pinfo['procs' ]['SF']
        self.f1_sw_ppm   = self.f1_sw_hz/self.f1_omega
        
        self.f2_offset_p = pinfo['proc2s' ]['OFFSET']
        self.f2_sw_hz     = pinfo['proc2s' ]['SW_p']
        self.f2_omega    = pinfo['proc2s' ]['SF']
        self.f2_sw_ppm   = self.f2_sw_hz/self.f2_omega
        print self.f1_sw_ppm
        
        self.observe_nucleus = info['acqus']['NUC1']
        if info['acqus']['NUC2'] == "off":
            self.indirect_nucleus = info['acqus']['NUC1']
        else:
            self.indirect_nucleus = info['acqus']['NUC2']
            
#        print self.observe_nucleus, self.indirect_nucleus
        
        if axis_domains[1] == "w":
            self.f1 = np.linspace(self.f1_offset_p, self.f1_offset_p-self.f1_sw_ppm,  self.rr)
        else:
            self.f1 = np.linspace(0,self.rr, self.rr)
            
        if axis_domains[0] == "w":
            self.f2 = np.linspace(self.f2_offset_p, self.f2_offset_p-self.f2_sw_ppm,  self.cc)
        else:
            self.f2 = np.linspace(0,self.cc, self.cc)
        
        self.dw_f1_ppm = self.f1[1]-self.f1[0]
        self.dw_f2_ppm = self.f2[1]-self.f2[0]
        
        for r in range(self.rr):
    
            for c in range( self.cc):
        
                self.Y[r,c] = self.f1[r]
                self.X[r,c] = self.f2[c] 
                
        print dimensions_ppm       
        if dimensions_ppm == []:
            self.dimensions_index = scipy.array([0,self.rr-1,0,self.cc-1])
        else:
            if axis_domains[1] == "w":
                r1 = int( (dimensions_ppm[1]-self.f1_offset_p)/self.dw_f1_ppm)
                r2 = int( (dimensions_ppm[0]-self.f1_offset_p)/self.dw_f1_ppm)
            else:
                r1 = int( (dimensions_ppm[1]))
                r2 = int( (dimensions_ppm[0]))

            if  axis_domains[0] == "w":
                c1 = int( (dimensions_ppm[2]-self.f2_offset_p)/self.dw_f2_ppm)
                c2 = int( (dimensions_ppm[3]-self.f2_offset_p)/self.dw_f2_ppm)
            else:
                c1 = int( (dimensions_ppm[2]))
                c2 = int( (dimensions_ppm[3]))
                
            
            self.dimensions_index = scipy.array([r1,r2,c1,c2 ])
            
#        print "self.dimensions_index", self.dimensions_index
            
        self.Z1 = self.exp[self.dimensions_index[0]:self.dimensions_index[1],self.dimensions_index[2]:self.dimensions_index[3]]
        self.X1 =   self.X[self.dimensions_index[0]:self.dimensions_index[1],self.dimensions_index[2]:self.dimensions_index[3]]
        self.Y1 =   self.Y[self.dimensions_index[0]:self.dimensions_index[1],self.dimensions_index[2]:self.dimensions_index[3]]


    def savefig( self, filename ):

        self.ppplot.savefig( filename)
        
        
        
    def create_plot_layout( self, dimensions_index, axis_domains):
        
#        print "dimensions_index",dimensions_index
        
        nullfmt   = NullFormatter()         # no labels

        # definitions for the axes
        left, width = 0.20, 0.55
        bottom, height = 0.2, 0.55
        bottom_h = left_h = left+width+0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.18]
        rect_histy = [left_h, bottom, 0.18, height]

        # start with a rectangular Figure
        self.ppplot = plt.figure(2, figsize=self.plot_size)

        self.axScatter = plt.axes(rect_scatter)
        self.axHistx = plt.axes(rect_histx)
        self.axHisty = plt.axes(rect_histy)

# no labels
        self.axHistx.xaxis.set_major_formatter(nullfmt)
        self.axHisty.yaxis.set_major_formatter(nullfmt)
        
        self.axScatter.tick_params(axis='x', labelsize=11)
        self.axScatter.tick_params(axis='y', labelsize=11)
        if axis_domains[0] == "w":
            self.axScatter.set_xlabel(self.nucleus_labels[self.observe_nucleus]+'[ppm]',fontsize=12)
        else:
            self.axScatter.set_xlabel(self.nucleus_labels[self.observe_nucleus]+'[pts]',fontsize=12)
        #ax.set_xlim(-60, 60)
        if axis_domains[1] == "w":
            self.axScatter.set_ylabel(self.nucleus_labels[self.indirect_nucleus]+'[ppm]', fontsize=12)
        else:
            self.axScatter.set_ylabel(self.nucleus_labels[self.indirect_nucleus]+'[expt #]', fontsize=12)
            
        
        self.axHistx.axis('off')
        self.axHisty.axis('off')

        f1_start = self.f1[dimensions_index[0]]
        f1_end   = self.f1[dimensions_index[1]]

        f2_start = self.f2[dimensions_index[2]]
        f2_end   = self.f2[dimensions_index[3]]

        self.axScatter.set_ylim( (f1_start, f1_end) )
        self.axScatter.set_xlim( (f2_start, f2_end) )
        
        

    def plot_plots(self):
        
        # the scatter plot:
        cl = np.linspace(self.Z1.max()*self.contourlevels_range[0], self.Z1.max()*self.contourlevels_range[1],10)
#        print "Z1.shape",self.Z1.shape
        

        sum_f1 = self.Z1.sum(axis=0)
        max_f1 = self.Z1.max(axis=0)
#        print "len(sum_f1)",len(sum_f1)
        sum_f2 = self.Z1.sum(axis=1)
        max_f2 = self.Z1.max(axis=1)
#        print "len(sum_f2)",len(sum_f2)

        cset = self.axScatter.contour(self.X1, self.Y1, self.Z1,  cl,  colors='red')
        #
        
#        self.axHistx.plot(sum_f1, 'r-')
#        self.axHisty.plot(sum_f2,range(len(sum_f2)),'r')


#       self.axHistx.set_xlim( (0,len(sum_f1)-1) )
#       self.axHisty.set_ylim( (0,len(sum_f2)-1) )        

        self.axHistx.plot(max_f1, 'r-')
        self.axHisty.plot(max_f2,range(len(max_f2)),'r')


        self.axHistx.set_xlim( (0,len(max_f1)-1) )
        self.axHisty.set_ylim( (0,len(max_f2)-1) )        

        
if __name__ == "__main__":

    expt_no = "27"  
    
    #/home/Service/data/bdat06/eeh_20151007_biaxial_2H
    
    
    
    directory = os.path.join("V:", os.sep,"Backup", "Avance",  "bdat06","eeh_20151007_biaxial_2H" )
    expt_no = "10"  

    directory = os.path.join("V:", os.sep,"Backup", "Avance",  "bdat08","eeh_20160508" )
    
    file_rrr = os.path.join( directory, expt_no, "pdata", "1")
    
    file_ser = os.path.join(directory,expt_no)
                             
    cs_file  = os.path.join( directory,expt_no,"vdlist")
    file1    = os.path.join( directory,expt_no,"pdata","1","2rr" )
    
    pinfo,expt = nmrglue.bruker.read_pdata( file_rrr )
    info,expt_fid = nmrglue.bruker.read( file_ser )
    

    
    pplot = TwoD_NMR_MAT_plot( expt_fid, pinfo, info, axis_domains = ["t","t"])
#    print "dir(pplot)",dir(pplot)
    pplot.plot_size=[5,5]

    pplot.observe_nucleus = info['acqus']['NUC1']
    
#    if info['acqus']['NUC2'] == "off":
#        pplot.indirect_nucleus = info['acqus']['NUC1']
#    else:
#        pplot.indirect_nucleus = info['acqus']['NUC2']
        
#    print pplot.observe_nucleus, pplot.indirect_nucleus
    pplot.display_plot()
    pylab.show()

        