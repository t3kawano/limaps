# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 13:14:31 2017

"""
"""
###########################################################################

Copyright (c) 2018, Taizo kawano

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

###########################################################################
lethargus and inter mold duration analysis python script

20180223
start develop based on laps.


This script process one experiment at a time and output the several measurements;
lethargus duration, total quiescent, number of bout etc.
For further analysis and figure preparation, 
you should use your own program, excel or whatever you like.

How to use
Before run this script, Prepare .csv file using imagesubtandmeaaure.py
Name the .csv file as "date_experimentname_experimentnumber_interval.xls" format.
e.g. 170719_n2rem5comparison_1_2.xls

Put limaps.py and limaps_classes.py (and dotplot.py, if you want) into a folder.
Run this script. 
It will show file choose dialog.
Choose the .csv file you made.
The script read the file and shows foq graphs of each sample.


The script outputs following files.

1. date_groupname_experimentnumber_samplenumber_foq.png
    Each sample's graph of fractoin of q

2. date_groupname_experimentnumber_samplenumber_foq.csv
    Each sample's fractoin of q.
    The second and third rows indicate beginning and end of lethargus

3. date_groupname_experimentnumber_summary.csv
    Table of measurements.
    fqlt_start: start point of lethargus
    fqlt_end: end point of lethargus
    fqlt_duration: duration (hrs) of lethargus
    totalq: total quiescence during lethargus(min)
    meanquiescent: quiescent time/lethargus duration (totalq/fqlt_duration/60)
    meanquiescentout: quiescent time/non-lethargus time
    numberofbout: quiescent number during lethargus
    qmean: mean duration of quiescent bout (sec)
    amean: mean duration of active bout(sec)
    qmedian: median duration of quiescent bout (frame)
    amedian: median duration of active bout (frame)
    qfreq: numberofbout/fqlt_duration (number/hr)

4. date_groupname_experimentnumber_fqalll_lethargus.png
    A graph contains foq of each sample.
    
5. date_groupname_experimentnumber_fqaligned_lethargus.png
    foq aligned at beginning of lethargus
    dotted lines: each sample
    black line: mean
    dashed line: mean +- sd
    
6. date_groupname_experimentnumber_fqalignedtail_lethargus.png
    foq aligned at end of lethargus

7. date_groupname_experimentnumber__fqlt_df.csv
    foq table during lethargus +- 1hr
    
The script may fail to detect lethargus correctly,
You have to be careful to use such data for analysis.    

20240517
    fix some comments

20240304
    separate class files

20200603
    qapixthreash to eliminate larvae contaminate
20200206 
    added code for saving samplegroups using pickle
20180123 ver.1.5.1
    lethargus detection criteria; longer than 2hrs.

20171211 ver.1.5
    1. foqthreshold is changeable at the top of the script.
    2. fqlt_df.csv and summary.csv contains interval data at the filename

20171122 ver.1.4
    1. foq plot of all sample 
    2. area rate of all sample 
    are saved

20170823 ver.1.3 
    1. save foq csv all data
    2. plot and save area figure 
    3. samplegroup functions for manual analysis;

    showcandidates(self, _samplenum): This method create a plot of 
    foq and onset exit candidate. also print the number of frame.
    
    manualanalysis(self, _samplenum, _start, _end): input the frame number.
    
    eg.
    : sg.showcandidates(1)
    start 764,7497
    end 235,6026,8018 
    :sg.manualanalysis(1,764,6026)
    totalq 95 min
    qend + aduration 7372.0
    ltend+inoutshift 6326
    trim the lastrow data
    numberofbout 306
    qmean 18.61sec
    amean 15.27sec
    qmedian 3.0frame
    amedian 2.0frame
    qfreq 104.68/hour
    meanquiescent 0.542189281642
    meanquiescentout 0.0067892790191
    ------------------
    man_start 764
    man_end 6026
    man_duration 2.9233333333333333
    totalq 95.1
    meanquiescent 0.542189281642
    meanquiescentout 0.0067892790191
    numberofbout 306
    qmean 18.607843137254903
    amean 15.265573770491804
    qmedian 3.0
    amedian 2.0
    qfreq 104.67502850627137
    
20170810 ver.1.2 clicking outside of lethargus period ignore the sample
20170807 ver.1.1 handle .csv as well
20170720 ver.1

#need to comment out following lines to use other env before upload.
#dont forget!!
#sys.path.append(os.path.join(os.path.dirname(__file__), 
#                                "..","tkmodules"))


"""

import os
import sys
import time
import datetime
import re
import csv
import pandas
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import tkinter
import tkinter.filedialog

import limaps_classes

#region parameters 
#180309 following usr defined parameters has to be set

#171211 implemented foqthreshold. 
#previously it was 0.05. 
#0.2 seems gives consistent result? when use low mag 0.8x
foqthreshold = 0.2
#foqthreshold = 0.3#ire-1;pYH330 rescue? -> still significant.seems ok
#foqthreshold = 0.5#200214 0.5 to detect only high foq part of rem5
#l3l4 24x hole use 0.05? rpi also 0.05 might better
#foqthreshold = 0.05
#foqthreshold = 0.02#221206 for cmc 1x1 olympus x2

qapixthreash = 1

#210618 for sfiwt 4x obj high mag, increase size detect as moving 800 is ok?
#qapixthreash = 800

#for aged HIS, to ignore larvae from old adult put some value here?
#smaller animals tend affect much than bigger one... may not good for use?
#qapixthreash = 50

#240516 to change minimal duration of lethargus, 
#change the minduration = 1.5 of class Samplegroup in limaps_classes.py

#region experimental sample settings
#180309 to handle an experiments contains multiple genotype..
#in the future, set these info at imageprocessing part might be better.
#6x8 config
colnum = 8
rownum =6
#6x6
#colnum = 6
#rownum =6
#8x6 vertical
#colnum = 6
#rownum =8

#group order horizontal or vertical?
grouporder = "h"
#grouporder = "v"

uniquegroupnames = ["n2","rem5","pek1","rem5pek1"]
uniquegroupnames = ["rem5","s206_3_3","s324_7","s324_11",
                    "s326_8","empty"]
uniquegroupnames = ["s337","s338","s339","s340"]



#in this case from 1st to 2nd columns are initial group, 
#3rd and 4th are another group and so on.
gindex = [(1,2),(3,4),(5,6),(7,8)]
gindex = [(1),(2),(3),(4),(5),(6),(7),(8)]

#180309 above usr defined parameters has to be set

#region prep index
#################################################################
indexgrid = np.array(np.arange(colnum*rownum)+1).reshape(rownum,colnum)
if grouporder == "h":
    indexgrid = np.fliplr(indexgrid.T)
"""
gsamplenumlist = []
for i in range(len(gindex)):
    ilist = list(range(gindex[i][0]-1,gindex[i][1]))
    templist = []
    for j in ilist:
        templist.extend(indexgrid[j,])
    print(uniquegroupnames[i]+" "+str(templist))
    gsamplenumlist.append(templist)
"""
#181226 mod for single column group
gsamplenumlist = []
for i in range(len(gindex)):
    if type(gindex[i]) is int:
        ilist = list(range(gindex[i]-1, gindex[i]))
    else:
        ilist = list(range(gindex[i][0]-1,gindex[i][1]))
    templist = []
    for j in ilist:
        templist.extend(indexgrid[j,])
    print(uniquegroupnames[i]+" "+str(templist))
    gsamplenumlist.append(templist)



scriptname = "limaps"
print(scriptname)
#print(str(time.time()))
print(str(datetime.datetime.today()))
print(os.getcwd())


# file choose dialog
tk = tkinter.Tk()
tk.withdraw()

# choose .csv file made by imagesubtandmeasure.py
targetfile = tkinter.filedialog.askopenfilename()

print("targetfile "+targetfile)
#if the file is not .csv/xls stop the process

targetdir = "/".join(targetfile.split("/")[0:-1])
targetextention = ""
separater = ""

if (not ".xls" in targetfile.split("/")[-1]) and (not ".csv" in targetfile.split("/")[-1]):
    print("this is not .xls file")
    sys.exit()
else:
    if ".xls" in targetfile.split("/")[-1]:
        targetextention=".xls"
        separater = "\t"
    if ".csv" in targetfile.split("/")[-1]:
        targetextention=".csv"
        separater = ","
        


#the file name should have date_geynotype_interval format, use the info
if targetfile.split("Result")[0]!="":
    #this format must concider bit more. "Result" may not included.
    #".xls" extension automaticcaly added at imagej.
    #date_sample goupe name_number(incase multiple sample at same day)_interval
    #filenameinfo = targetfile.split("/")[-1].split("Result")[0].split("_")
    filenameinfo = targetfile.split("/")[-1].split(targetextention)[0].split("_")
    thedate = filenameinfo[0]
    #240229 gropuname here is confusing. actually experiment name?
    #240516 changed to experiment_name
    experiment_name = filenameinfo[1]
    expnum = int(filenameinfo[2])
    #20180201 need in case that msec. also interval as float
    # msec interval suffix m. eg. 180201_n2_1_500m
    if "m" in filenameinfo[3]:
        #interval = int(filenameinfo[3])
        interval = float(filenameinfo[3].strip("m"))/1000
    else:
        interval = float(filenameinfo[3])
        

print("thedate: " + thedate + "\ngexperiment_name: "+experiment_name+
      "\nexpnum: "+ str(expnum)+"\ninterval: " + str(interval))
    
df = pandas.read_csv(targetfile, sep= separater)
runmedwindow = int(60/interval)+1

#extruct area data only
areacolnames = []
for a in df.columns:
    if "Area" in a:
        print(a)
        areacolnames.append(a)
areacolnames

areadf = df[areacolnames]

# 240222 considering separate configs in a limaps_config.py file.
# If the config file exist in the same dir, use it. -> not impolemented yet





#sys.exit()




#20180309 seting group name here is pending for now
"""
#20170925 to handle data contains multiple samplegroup.
# the roi name should contains groupname and separate the number by _
# e.g. n2_1
uniquegroupnames = []
groupnamearray =[]

# if the colname contains area()xx_, like Area(n2_1)...
if len(re.split("\(|\)|_",areacolnames[0])) > 1:
    groupnamearray = [re.split("\(|\)|_",a)[1] for a in areacolnames]
    for i in groupnamearray:
        if i not in uniquegroupnames:
            uniquegroupnames.append(i)
else:
    uniquegroupnames.append(groupname)
    groupnamearray = [groupname for a in range(len(areadf.columns))]
            
print(uniquegroupnames)
print(groupnamearray)
"""

#areadf.iloc[:,["test" == a for a in groupnamearray]]

"""
make class storing a data of an individual
Individual
date
sample groupe name(genotype andor treatments)
interval

rawdata
qaboolean
foq 10min
normalized area
active rate 1min

foq depend lethargus
onset,exit candidates
final lethargus onset, exit

active rate depend lethargus
onset,exit candidates
final lethargus onset, exit

output each figure.
foq values
lethargus on, off points

about durations, qbout durations,
mean foq at lethargus, out of lethargus
totalq, let duration, transition #


[list of the class]
-> make summary figure
foq graph
    arrign at onset, and exit
-> make csv files
qaboolean
foq values
about durations, qbout durations,
mean foq at lethargus, out of lethargus
totalq, let duration, transition #

    


"""

#Samplegroup.makeamatrix = makeamatrix		
#Samplegroup.deletealethargus = deletealethargus
#Samplegroup.appendanindividual = appendanindividual
#Samplegroup.calcmaxduration =calcmaxduration
#Samplegroup.makeallfigurefoq = makeallfigurefoq
#Samplegroup.makeallfigureofrawarea = makeallfigureofrawarea
#Samplegroup.prepallfigureformat = prepallfigureformat
#Samplegroup.calcinitialqlatency = calcinitialqlatency
#Samplegroup.makeallfigurewithinitialq = makeallfigurewithinitialq
#Samplegroup.manualanalysis = manualanalysis
#Samplegroup.saveafig = saveafig      
#Samplegroup.deleteanindividual = deleteanindividual      
#Samplegroup.makedf = makedf      
#Samplegroup.savesummarydf = savesummarydf      
#Samplegroup.makeallfigure = makeallfigure      
#Samplegroup.saveafig = saveafig      
#Samplegroup.makeltmatrix = makeltmatrix      
#Samplegroup.makealignfigure = makealignfigure      
#Samplegroup.saveltmatrix = saveltmatrix      
#Samplegroup.makeallfigureofarea = makeallfigureofarea      
#Samplegroup.makeallfigureofareacropped = makeallfigureofareacropped      
#Samplegroup.showcandidates = showcandidates      
#Samplegroup.processagroup = processagroup      

#sg.calcinitialqlatency()
#sg.makeallfigurewithinitialq()
#sg.showcandidates(2)
#sg.makeallfigureofarea(30)
#sg.makeallfigureofareacropped()
        
#df = sg.makedf()
#sg.savesummarydf()
#fig = sg.makeallfigure()
#sg.saveafig(fig)
#temp = sg.makeltmatrix(True)
#temp = sg.makeltmatrix(False)
#temp = sg.ltfoqmatrixtail

#areadf.columns[0]

#sys.exit()


#region sample group process

samplegroups = []

#for gn in uniquegroupnames:
for (gn, indice) in zip(uniquegroupnames, gsamplenumlist):
    print(gn, indice)
    #this groupname format is date_groupname_expnum
    #this may not good? just use groupname?
    sg = limaps_classes.Samplegroup("_".join([thedate, gn, str(expnum)]),targetdir)
    #print(uniquegroupnames[i]+" "+str(templist))
    #gsamplenumlist.append(templist)

    #tempdf = areadf.loc[:,[gn == a for a in groupnamearray]]
    tempdf = areadf.iloc[:,np.array(indice)-1]
    
    #Individual.plotlethargus = plotlethargus
    #def processagroup(self, thedate, groupname, expnum, interval, dataframe):    
    #sg.processagroup(thedate, gn, expnum, interval,
    #                 tempdf, threshold= foqthreshold)
    sg.processagroup(thedate, gn, expnum, interval,
                     tempdf, threshold= foqthreshold, pixthreash = qapixthreash)
    

    #summary output    
    df = sg.makedf(ltlist = sg.fullindlist)
    if not df.empty:
        sg.savesummarydf()
    
    #fig = sg.makeallfigure()
    #sg.saveafig(fig, "fqalll")
    
    #temp = sg.makeltmatrix(True)
    temp = sg.makeltmatrix(ordinalnum = 0)
    tempfig = sg.makealignfigure(sg.ltfoqmatrix, True, "mean")
    sg.saveafig(tempfig, "fqaligned")
    
    temp = sg.makeltmatrix(align = "tail", ordinalnum = 0)
    temp = sg.ltfoqmatrixtail
    tempfig = sg.makealignfigure(sg.ltfoqmatrixtail, False, "mean")
    sg.saveafig(tempfig, "fqalignedtail")
    
    if not df.empty:
        sg.saveltmatrix(sg.ltfoqmatrix, "fqlt",ordinalnum = 0)
    
        areacropfig = sg.makeallfigureofareacropped(sg.indlist, ordinalnum = 0)
        sg.saveafig(areacropfig, "areacropped")
    #sg.saveafig(areacropfig, "areacropped_del21")
    """
    tempsg = samplegroups[0]
    areacropfig = tempsg.makeallfigureofareacropped(tempsg.indlist, ordinalnum = 0)
    tempsg.saveafig(areacropfig, "areacropped")
    """

    fullfoqfig = sg.makeallfigurefoq(sg.fullindlist)#
    fullfoqfig.subplots_adjust(hspace = .3)
    fullfoqfig.subplots_adjust(wspace = .01)
    sg.saveafig(fullfoqfig, "foqfull")

    areafig = sg.makeallfigureofrawarea(60, sg.fullindlist,"median")
    areafig.subplots_adjust(hspace = .1)
    areafig.subplots_adjust(wspace = .01)                                    
    sg.saveafig(areafig, "rawareafull")

    """
    fullfoqwolfig = sg.makeallfigurefoq(sg.fullindlist,
                             overlay = [1,2])#
    sg.saveafig(fullfoqwolfig, "foqfullwol")

    fullfoqwolfig.subplots_adjust(hspace = .1)
    fullfoqwolfig.subplots_adjust(wspace = .01)

    fullfoqfig = sg.makeallfigurefoq(sg.fullindlist,
                             overlay = [[0,0.25],[1,1.25],[2,2.25],[3,3.25],[4,4.25]])#
    fullfoqgridfig = sg.makeallfigurefoq(sg.fullindlist,
                                     figsize = (16,4),
                                     col = 6, row = 4,
                                     labeloff = True)
                                     #overlay = [[0,0.25],[1,1.25],[2,2.25],[3,3.25],[4,4.25]])
                                     #overlay = [1,3])#8x6
                                     #overlay = [1,2])#8x6
    fullfoqgridfig.subplots_adjust(hspace = .1)
    fullfoqgridfig.subplots_adjust(wspace = .01)
    sg.saveafig(fullfoqgridfig, "foqfullgrid")
    """
    #sg.saveafig(fullfoqfig, "foqfull")
    """
    fullfoqfig = sg.makeallfigurefoq(sg.fullindlist[0:16], 
                                     lethargus ="off")#
    sg.saveafig(fullfoqfig, "foqfullsub")
    """    
    #areafig = sg.makeallfigureofrawarea(60, sg.fullindlist, "median")#60 sec
    #areafig = sg.makeallfigureofrawarea(3, sg.fullindlist, "median")#6 sec
    """
    areafig = sg.makeallfigureofrawarea(60, sg.fullindlist[0:16],
    areafig = sg.makeallfigureofrawarea(60, sg.fullindlist,
                                        "median",
                                        overlay = [[0,0.25],[1,1.25],[2,2.25],[3,3.25],[4,4.25]])
                                        #overlay = [4,5])#60 sec
    areafig.subplots_adjust(hspace = .1)
    areafig.subplots_adjust(wspace = .01)
    sg.saveafig(areafig, "rawareafullsub")
    """
    """
    #initial 2 hours. for heat induced sleep
    areagridfig = sg.makeallfigureofrawarea(60, sg.fullindlist, "median",
                                        figsize = (16,4),
                                        col = 6, row = 4,
                                        xlim = [0, 1*60*60/interval],
                                        labeloff = True,
                                        overlay = [1,3])#120sec, 8x6
    areagridfig.subplots_adjust(hspace = .1)
    areagridfig.subplots_adjust(wspace = .01)
    sg.saveafig(areagridfig, "rawarea2hrsgrid")
    
    fullfoqgridfig = sg.makeallfigurefoq(sg.fullindlist,
                                     figsize = (16,4),
                                     col = 8, row = 6,
                                     xlim = [0, 1*60*60/interval],
                                     labeloff = True)
                                     #overlay = [1,3])#8x6
    fullfoqgridfig.subplots_adjust(hspace = .1)
    fullfoqgridfig.subplots_adjust(wspace = .01)
    sg.saveafig(fullfoqgridfig, "foq2hrsgrid")
    """    
    """
    areagridfig = sg.makeallfigureofrawarea(60, sg.fullindlist, "median",
                                        figsize = (16,4),
                                        col = 6, row = 4,
                                        labeloff = True)
                                        #overlay = [[0,0.25],[1,1.25],[2,2.25],[3,3.25],[4,4.25]])
                                        #overlay = [1,3])#120sec, 8x6
    areagridfig.subplots_adjust(hspace = .1)
    areagridfig.subplots_adjust(wspace = .01)
    #sg.saveafig(areagridfig, "rawarea120fullallgrid")
    sg.saveafig(areagridfig, "rawarea60fullallgrid")
    """
    """
    sg.saveafig(areafig, "rawareafullallsample")
    
    np.nanmean(sg.summarydf.fqlt_duration)
    np.nanstd(sg.summarydf.fqlt_duration)
    """
    samplegroups.append(sg)


# 200206 using pandas and pickle to save samplegroups object made by limaps
objfilename = "_".join(filenameinfo)+"_samplegroups.pkl"
objfilepath = os.path.join(targetdir,objfilename )
pandas.to_pickle(samplegroups,objfilepath)
############################################################################





#gridindex; start at top left and goes right (row)
#same order of pyplot subplot number
def gridtosgindex(gridindex):
    gnum = None
    theindex = None
    for i, alist in enumerate(gsamplenumlist):
        if gridindex in alist:
            gnum = i
            theindex = alist.index(gridindex)
            #returns group index and sample index within the group
            return gnum, theindex

def plotgridfig(**kwargs):
    gridfig = plt.figure(figsize = (16,6))
    plotdata = "foq"
    if "plotdata" in kwargs:
        plotdata = kwargs["plotdata"]
    window = 60
    if "window" in kwargs:
        window = kwargs["window"]
    overlayparam = None
    if "overlayparam" in kwargs:
        overlayparam = kwargs["overlayparam"]
    meanduration = 2.64
    if "meanduration" in kwargs:
        meanduration = kwargs["meanduration"]
    sdduration = 0.11
    if "sdduration" in kwargs:
        sdduration = kwargs["sdduration"]
    if "xlim" in kwargs:
        xlim = kwargs["xlim"]
    else:
        xlim = None
    for i in range(colnum * rownum):
        #figmultifoq.clear()
        ax = gridfig.add_subplot(rownum,colnum,i+1)
        sgn, ti = gridtosgindex(i+1)
        _ind = samplegroups[sgn].fullindlist[ti]
        if plotdata == "foq":
            ax.plot(_ind.foq,
                    linestyle = "-", linewidth = 1,
                    color = "black")
            if (_ind.letharguslist is not None) and (len(_ind.letharguslist)>0):
                for lt in _ind.letharguslist:
                    ax.axvline(x= lt.start, color= "black")
                    ax.axvline(x= lt.end, color= "black", linestyle = "-")
        elif plotdata == "area":
            tempraw = _ind.calcrawarea(60)
            ax.plot(tempraw.rolling(window=int(window/_ind.interval),
                                     center= True).median().values,
                    linestyle = "-", linewidth = 0.5,
                    color = "black")
        if overlayparam is not None:
            if overlayparam == "fq_duration":
                if (_ind.letharguslist is not None) and (len(_ind.letharguslist)>0):
                    lt = _ind.letharguslist[0]
                    fontcolor = "black"
                    if lt.fq_duration > meanduration + 3*sdduration:
                        fontcolor = "red"
                    elif lt.fq_duration < meanduration - 3*sdduration:
                        fontcolor = "blue"
                    ax.annotate("{0:4.3}".format(lt.fq_duration),
                                xy = (lt.start, 0.1),
                                fontsize=10,
                                color = fontcolor)
                

        ax.set_ylim(0, 1.1)
        #ax.set_xlim(0, 5*60*60/interval)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        hrtick = np.arange(len(_ind.foq)*_ind.interval/60/60).astype(int)
        ax.set_xticks(hrtick*60*60/_ind.interval)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if xlim is not None:
            ax.set_xlim(0, xlim)
            #print("xlim set")
        #ax.set_xlabel("time (h)")
        #ax.legend(groupnames)
        label = "_".join([_ind.date, _ind.groupname, str(_ind.expnum),
                                  str(_ind.samplenum)])
        
        ax.annotate(label,\
                        xy=(0.01, 0.9),\
                         xycoords='axes fraction', fontsize=8,\
                         horizontalalignment='left', verticalalignment='bottom')
        
                
    ax.set_xticklabels(hrtick)
    gridfig.tight_layout()
    gridfig.subplots_adjust(hspace = .1)
    gridfig.subplots_adjust(wspace = .01)
    return gridfig

gridfig = plotgridfig(plotdata = "foq")

picformat = ".png"
filename = "_".join([thedate, experiment_name, "gridfoq", picformat])
figfilepath = os.path.join(targetdir,filename)
gridfig.savefig(figfilepath,dpi=100)

#240201 for svg format
picformat = ".svg"
filename = "_".join([thedate, experiment_name, "gridfoq", picformat])
figfilepath = os.path.join(targetdir,filename)
gridfig.savefig(figfilepath,dpi=100)


gridfig = plotgridfig(plotdata = "area", window = 120)
gridfig.savefig(targetdir +"/"+"gridarea.png",dpi=100)

# for heat induced sleep and hsp over expression
stimflag = False
#stimflag = True
crophrs = 6
#crophrs = 3
cropframe = int(crophrs*60*60/2)
if stimflag:
    gridfigxh = plotgridfig(plotdata = "foq", xlim = 1800*3)
    gridfigxh.savefig(targetdir +"/"+"gridfoq{0}h.png".format(crophrs),dpi=100)

    
    gridfigxh = plotgridfig(plotdata = "area", window = 120, xlim = 1800*3)
    gridfigxh.savefig(targetdir +"/"+"gridarea{0}h.png".format(crophrs),dpi=100)
    
    foqmatlist = []
    areamatlist = []
    for sg in samplegroups:
        #sg.indlist[0].
        foqfigcrop = sg.makeallfigurefoq(sg.fullindlist, xlim=(0,cropframe))
        foqfigcrop.subplots_adjust(hspace = .3)
        foqfigcrop.subplots_adjust(wspace = .01)
        sg.saveafig(foqfigcrop, "foq{0}hrs".format(crophrs))
        
        #temmat = sg.makeamatrix(datatype = "foq", xlim =(0, cropframe))
        #180904 
        #to eliminate some failed saples.... need to impliment later
        temmat = sg.makeamatrix(datatype = "foq", xlim =(0, cropframe),autoelim = True)

        tempmatfig = sg.makealignfigure(temmat, True, "mean")
        sg.saveafig(tempmatfig, "fqalignedat0_{0}hr".format(crophrs))
        foqmatlist.append(temmat)
        
        temmat = sg.makeamatrix(datatype = "area", xlim =(0, cropframe))
        tempmatfig = sg.makealignfigure(temmat, True, "mean")
        sg.saveafig(tempmatfig, "areaalignedat0_{0}hr".format(crophrs))
        areamatlist.append(temmat)

    #at time point 5min, 10 min...
    timepoints = [5,10,15,20,25,30,35,40,45,50,55,60]
    for t in timepoints:
        data0 = foqmatlist[0][:,30*t]
        data1 = foqmatlist[1][:,30*t]
        stat, pval = stats.ttest_ind(data0, data1, equal_var = False)
        print("foq {0}min pval ".format(t), str(pval))
    for t in timepoints:
        data0 = areamatlist[0][:,30*t]
        data0 = data0[~np.isnan(data0)]
        data1 = areamatlist[1][:,30*t]
        data1 = data1[~np.isnan(data1)]
        stat, pval = stats.ttest_ind(data0, data1, equal_var = False)
        print("area {0}min pval ".format(t), str(pval))
        
#gridfig = plotgridfig(overlayparam = "fq_duration")
#gridfig.savefig(targetdir +"/"+"gridfoqwithflag.png",dpi=100)



#180501 for mt screen. 

#Put NaN where the col value are outlier.
#upper = quantile[0.75]+1.5*(q0.75-q0.25)
#lower = quantile[0.25]-1.5*(q0.75-q0.25)
def deloutlier(df, col):
    thequantiles = df.loc[:,[col]].quantile([.25, .5, .75])
    upperlimit = thequantiles.loc[0.75]+1.5*(thequantiles.loc[0.75]-thequantiles.loc[0.25])
    lowerlimit = thequantiles.loc[0.25]-1.5*(thequantiles.loc[0.75]-thequantiles.loc[0.25])
    #return upperlimit, lowerlimit
    dfcopy = df.copy()
    amask = dfcopy.loc[:,["fqlt_duration"]]>upperlimit
    dfcopy.loc[:,["fqlt_duration"]]=dfcopy.loc[:,["fqlt_duration"]].mask(amask)
    amask = dfcopy.loc[:,["fqlt_duration"]]<lowerlimit
    dfcopy.loc[:,["fqlt_duration"]]=dfcopy.loc[:,["fqlt_duration"]].mask(amask)
    
    return dfcopy, upperlimit, lowerlimit

#dfwool, upperlimit, lowerlimit = deloutlier(durationdf, "fqlt_duration")
#meanwoolduration = dfwool.loc[:,["fqlt_duration"]].mean()
#sdwoolduration = dfwool.loc[:,["fqlt_duration"]].std()
 



param = "fqlt_duration"

drdf = pandas.DataFrame(columns = ["groupname",param])
for sg in samplegroups: 
    if sg.summarydf is not None:
        #drdf = drdf.append(sg.summarydf.loc[:,["groupname",param]])
        drdf = pandas.concat([drdf,sg.summarydf.loc[:,["groupname",param]]])

#grouplist = anadf.groupname.unique()

dfwool, upperlimit, lowerlimit = deloutlier(drdf, "fqlt_duration")

#240220 nan cause error? no groupname is unnecessary
#meantheex = dfwool.loc[:,["groupname",param]].mean().values[0]
meantheex = dfwool.loc[:,[param]].mean().values[0]
#stdtheex = dfwool.loc[:,["groupname",param]].std().values[0]
stdtheex = dfwool.loc[:,[param]].std().values[0]
gridfig = plotgridfig(overlayparam = "fq_duration", 
                      meanduration = meantheex)


picformat = ".png"
filename = "_".join([thedate, experiment_name, "gridfoqwithflag_meancalc", picformat])
figfilepath = os.path.join(targetdir,filename)
gridfig.savefig(figfilepath,dpi=100)


#240201 for svg format
picformat = ".svg"
#filename = "_".join(["foqmarged_mean","marging",str(marginhr),
#                     ugname, timesuff, picformat])
filename = "_".join([thedate, experiment_name, "gridfoqwithflag_meancalc", picformat])
figfilepath = os.path.join(targetdir,filename)
gridfig.savefig(figfilepath,dpi=100)


#sys.path.append is not recommended way to set path for module... 
# but for now use easy way
# In the readme.md on git, suggesting that 
# put dotplot.py in the same directory.
# so this line is only for my devlopment env.
#sys.path.append(os.path.join(os.path.dirname(__file__), 
#                                "..","tkmodules"))
#import importlib
#importlib.reload(dotplot)
import dotplot
dotplotfig = dotplot.dotplots(drdf,outlier = "oc",
                              col="gray",
                              ylim = (0,max(drdf[param])),
                              binnum = 25, size = 5, thickness = 1.5,
                              figsize = (2,4))
                              #figsize = (4,8))


dotplotfig.axes[0].set_ylabel(param)
#dotplotfig.axes[0].annotate("p={0:<8.3g}".format(pval),
#uplvalue = upperlimit.values[0]
#dotplotfig.axes[0].annotate("upper={0:<3.2g}".format(uplvalue),
#                 xy=(0.03, 0.92),\
#                 xycoords='axes fraction', fontsize=10,\
#                 horizontalalignment='left', verticalalignment='bottom')
dotplotfig.tight_layout()

picformat = ".png"
filename = "_".join([thedate, experiment_name, "dotplot", param, "all", picformat])
figfilepath = os.path.join(targetdir,filename)
dotplotfig.savefig(figfilepath,dpi=100)


#240201 for svg format
picformat = ".svg"
#filename = "_".join(["foqmarged_mean","marging",str(marginhr),
#                     ugname, timesuff, picformat])
filename = "_".join([thedate, experiment_name, "dotplot", param, "all", picformat])
figfilepath = os.path.join(targetdir,filename)
dotplotfig.savefig(figfilepath,dpi=100)





#sys.exit()
#############################################################################


