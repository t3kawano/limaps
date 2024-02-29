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
Before run this script, Prepare .xls file using imageJ's Multi measure function.
Name the .xls file as "date_groupname_experimentnumber_interval.xls" format.
e.g. 170719_oy59_1_2.xls

Open Spyder and laps.py (this script).
Run this script (click triangle or push F5). 
It will show file choose dialog.
Choose the .xls file you made.
The script read the file and shows foq graphs of each sample.

In some case, it will detect multiple possible lethargus periods.
You have to choose the one that most likely to be correct by clicking the graph.
Also sometimes it mistakenly detect quiescent period as lethargus.
You should always inspect the result by your eyes,
 and eliminate such data from further analysis.

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
    
As mentioned above, the script sometimes fails lethargus detection,
so these files contains incorrect data.
You have to eliminate the incorrect data for publication quality analysis.    

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

#180309 following usr defined parameters has to be set
#from scipy import stats
#171211 implemented foqthreshold. 
#previously it was 0.05. 
#0.2 seems gives consistent result? when use low mag 0.8x
foqthreshold = 0.2
#foqthreshold = 0.3#ire-1;pYH330 rescue?
#foqthreshold = 0.5#200214 0.5 to detect only high foq part of rem5
#l3l4 24x hole use 0.05? rpi also 0.05 might better
#foqthreshold = 0.05
#foqthreshold = 0.02#221206 for cmc 1x1 olympus x2

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


uniquegroupnames = ["num1","num2","num3","num4","num5","num6","num7","num8"]
uniquegroupnames = ["n2","rem5","pek1","rem5pek1"]

qapixthreash = 1

#210618 for sfiwt 4x obj high mag, increase size detect as moving 800 is ok?
#qapixthreash = 800

#for aged HIS, to ignore larvae from old adult put some value here?
#smaller animals tend affect much than bigger one... may not good for use?
#qapixthreash = 50

#in this case from 1st to 2nd columns are initial group, 
#3rd and 4th are another group and so on.
gindex = [(1,2),(3,4),(5,6),(7,8)]
#sgindex = [(1),(2),(3),(4),(5),(6),(7),(8)]

#180309 above usr defined parameters has to be set

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

# choose .xls (actually csv) file made by imagejs multi measuer
targetfile = tkinter.filedialog.askopenfilename()

print("targetfile "+targetfile)
#if the file is not .xls, stop the process

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
        


#the file named should have date_geynotype_interval format, use the info
if targetfile.split("Result")[0]!="":
    #this format must concider bit more. "Result" may not included.
    #".xls" extension automaticcaly added at imagej.
    #date_sample goupe name_number(incase multiple sample at same day)_interval
    #filenameinfo = targetfile.split("/")[-1].split("Result")[0].split("_")
    filenameinfo = targetfile.split("/")[-1].split(targetextention)[0].split("_")
    thedate = filenameinfo[0]
    #240229 gropuname here is confusing. actually experiment name?
    groupname = filenameinfo[1]
    expnum = int(filenameinfo[2])
    #20180201 need in case that msec. also interval as float
    # msec interval suffix m. eg. 180201_n2_1_500m
    if "m" in filenameinfo[3]:
        #interval = int(filenameinfo[3])
        interval = float(filenameinfo[3].strip("m"))/1000
    else:
        interval = float(filenameinfo[3])
        

print("thedate: " + thedate + "\ngroupname: "+groupname+
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

#Lethargus class. contains area and foq data of individual,
# and calculate the lethargus measurement.
#180226 meanquiescentout is not correct way if multiple lethargus are deteced
class Lethargus:
    def __init__(self, _ind, lethargusperiod):
        self.ind = _ind
        self.start = lethargusperiod[0]
        self.end = lethargusperiod[1]
        #measurements
        self.meanquiescent = None
        #self.meanquiescentout = None
        self.totalq = None
        self.numberofbout = None
        self.qmean = None 
        self.amean = None
        self.qmedian = None 
        self.amedian = None 
        self.fq_duration = None
        self.qfreq = None
        self.qaboutdataframe = None
        self.lethargusperiodqadf = None
        
        #set true if distance is closer tha mindistance
        self.distanceproblem = False

    #if use fraction of q, inoutshift must be 10min. because foq is boxcar average
    #if areabased detection, inoutshift is 0. becase it running median
    #def calcfoqlethargusmeasures(self, _ltstart, _ltend):
    def calcmeasurements(self):
        inoutshift = int(600/self.ind.interval)
        #measures = self.calcmeasures(self.start+inoutshift,
        #                             self.end+inoutshift)
        measures = self.calcmeasures()
        #outputdata = [meanquiescent, meanquiescentout, totalq, numberofbout,
        #              qmean, amean, qmedian, amedian, fq_duration, qfreq,
        #              qaboutdataframe,lethargusperiodqadf]
        
        self.meanquiescent = measures[0]
        self.meanquiescentout = measures[1]
        self.totalq = measures[2]
        self.numberofbout = measures[3]
        self.qmean = measures[4] 
        self.amean = measures[5]
        self.qmedian = measures[6] 
        self.amedian = measures[7] 
        self.fq_duration = measures[8]
        self.qfreq = measures[9]
        self.qaboutdataframe = measures[10]
        self.lethargusperiodqadf = measures[11]
        #return outputdata

    #def calcmeasures(self, _ltstart, _ltend):
    def calcmeasures(self):
        #qa = calcqa(_rawdata, interval, bminmax)
        #shift 10min  when shift onset, also shift exit.
        #inoutshift = int(600/self.interval)
        #ltstart = self.fq_finallethargusperiod[0]
        #ltend = self.fq_finallethargusperiod[1]
        ltstart = self.start
        ltend = self.end
        #total quiescence during the lethargus
        #totalq = np.sum(self.qaboolean[(ltstart+inoutshift):
        #                            (ltend+inoutshift)])*self.interval/60
        totalq = np.sum(self.ind.qaboolean[ltstart:ltend])*self.ind.interval/60
        #totalq all imaging duration
        totalqall = np.sum(self.ind.qaboolean)*self.ind.interval/60
        #totalq out of lethargus
        totalqout = totalqall - totalq
        #here cause some trouble, is arawdata was seriise totalq is not list.s
        #totalq = totalq[0]
        print("totalq "+str(int(totalq))+" min")
    
        #calc q and a bout duration
        qabooleandiff = self.ind.qaboolean.astype(int).diff()
        qstart = np.where(qabooleandiff == 1)[0]
        qend = np.where(qabooleandiff == -1)[0]
    
        #fix always qstart < qend
        if qstart[0] > qend[0]:
            qend = qend[1:]
        if qstart[-1] > qend[-1]:
            qstart = qstart[:-1]
    
        qaboutdataframe = pandas.DataFrame(columns =["qstart",
                                                     "qend", "qduration","aduration"])
        qaboutdataframe["qstart"]=qstart
        qaboutdataframe["qend"]=qend
        qaboutdataframe["qduration"]=qend-qstart
        #.astype(int) convert nan to -2147483648
        #qaboutdataframe["aduration"]=np.append(qstart[1:]-qend[:-1], np.nan).astype(int)
        qaboutdataframe["aduration"]=np.append(qstart[1:]-qend[:-1], np.nan)
    
    
        #lethargusperiodboolean = (qaboutdataframe["qstart"] > ltstart+inoutshift)\
        #                            &\
        #                            (qaboutdataframe["qend"] < ltend+inoutshift)
        lethargusperiodboolean = (qaboutdataframe["qstart"] > ltstart)\
                                    &\
                                    (qaboutdataframe["qend"] < ltend)
        
        lethargusperiodqadf = qaboutdataframe[lethargusperiodboolean].copy()
        lastrow = lethargusperiodqadf.iloc[-1]
        #if active bout beyond lethargus, eliminate it
        print("qend + aduration "+str(lastrow["qend"] + lastrow["aduration"]))
        #print(" ltend+inoutshift " + str( ltend+inoutshift))
        print(" ltend " + str( ltend))
        #if (lastrow["qend"] + lastrow["aduration"]) >  ltend+inoutshift:
        if (lastrow["qend"] + lastrow["aduration"]) >  ltend:
            print("trim the lastrow data")
            dfname = lethargusperiodqadf.iloc[-1].name
            lethargusperiodqadf.loc[dfname,"aduration"] = np.nan
    
        numberofbout = len(lethargusperiodqadf)
        print("numberofbout " +str(numberofbout) )
    
        qmean = np.mean(lethargusperiodqadf["qduration"])*self.ind.interval
        amean = np.mean(lethargusperiodqadf["aduration"])*self.ind.interval
    
        print("qmean "+str(np.round(qmean,decimals=2))+"sec")
        print("amean "+str(np.round(amean,decimals=2))+"sec")
    
        qmedian = np.median(lethargusperiodqadf["qduration"].dropna())
        amedian = np.median(lethargusperiodqadf["aduration"].dropna())
        print("qmedian "+str(qmedian)+"frame")
        print("amedian "+str(amedian)+"frame")
      
    
        #frequency /hour
        fq_duration = (ltend-ltstart)*self.ind.interval/60/60
        qfreq = numberofbout/fq_duration
        print("qfreq " + str(np.round(qfreq,2)) +"/hour")
        #quiescent / time neary= mean foq
        meanquiescent = totalq/((ltend-ltstart)*self.ind.interval/60)
        print("meanquiescent " +str(meanquiescent))
        #180226 meanquiescentout is not correct way if multiple lethargus are deteced
        meanquiescentout = totalqout/((len(self.ind.qaboolean)-(ltend-ltstart))*self.ind.interval/60)
        print("meanquiescentout " +str(meanquiescentout))
        
        outputdata = [meanquiescent, meanquiescentout, totalq, numberofbout,
                      qmean, amean, qmedian, amedian, fq_duration, qfreq,
                      qaboutdataframe,lethargusperiodqadf]
        """
        self.meanquiescent = meanquiescent
        self.meanquiescentout = meanquiescentout
        self.totalq = totalq
        self.numberofbout = numberofbout
        self.qmean = qmean 
        self.amean = amean
        self.qmedian = qmedian 
        self.amedian = amedian 
        self.fq_duration = fq_duration
        self.qfreq = qfreq
        self.qaboutdataframe = qaboutdataframe
        self.lethargusperiodqadf = lethargusperiodqadf
        """

        return outputdata

#a kind of processor. when fed a Ind (contains area data),
# it will spit where is the lethargus
# self.ind.letharguslist = letharguslist
class Lethargusdetector:
    #foq threshod, 0.05 old way, 0.2 for 48x holes may be 0.1 is for 24x?
    #minimum duration. 1 or 2hrs?
    #minimul interval between lethargus. 6hrs? 10 hrs
    def __init__(self, foqthreshold, minduration, mininterval):
        self.foqthreshod = foqthreshold
        self.minduration = minduration
        self.mininterval = mininterval

    def setadata(self, ind):
        self.ind = ind
        self.interval = self.ind.interval
        
    def processdata(self, **kwargs):
        #here the method to detect lethargus. may be multiple.
        fig, ax, ax2 = self.prepfig()
        #ax2.plot(self.ind.foq, linewidth =0.5, 
        #        color = "black", linestyle ="-")
        ax.plot(self.ind.foq, linewidth =0.5, 
                color = "black", linestyle ="-")
        rawdata = self.ind.calcrawarea(60)
        ax.plot(rawdata.rolling(window = 60).median(), linewidth =0.5, 
                color = "gray", linestyle ="-")
        fig.tight_layout()
        ax2.axhline(y = self.foqthreshod, linewidth =0.5, 
                  color = "black", linestyle =":")
        
        self.detectlethargus1stscreen_fq()

        if len(self.ind.fq_onsetcandidates) >0:
            #screen periods that have 1h pre/post lethargus and over 1h duration
            
            self.ind.fq_oescreendmatrix = self.detectlethargus2ndscreen(self.ind.fq_qbooleanvec, 
                                                                        self.ind.fq_onsetcandidates,
                                                                        self.ind.fq_exitcandidates,
                                                                        prepostmargin = 0.5)
            if len(self.ind.fq_oescreendmatrix) > 0:
                self.plotcandidates(ax)
                #even after above, still may exist multiple. use human eyes.
                #20180223 to detect l3l4 duration at least two requied.
                """
                lths = self.detectlethargus(self.ind.fq_oescreendmatrix)
                if len(lths)>0:
                    for lt in lths:
                        rect = plt.Rectangle((lt[0],0), lt[1]-lt[0],0.9,alpha = 0.2, color = "blue")
                        ax.add_patch(rect)
                """
                
                letharguslist = []
                for lt in self.ind.fq_oescreendmatrix:
                    alt = Lethargus(self.ind, lt)
                    alt.calcmeasurements()
                    letharguslist.append(alt)
                
                self.checkdistance(letharguslist)
                while self.problematicdistance:
                    self.filterbytotalq(letharguslist)
                    
                if len(letharguslist)>0:
                    for lt in letharguslist:
                        rect = plt.Rectangle((lt.start,0), lt.end-lt.start,0.9,alpha = 0.2, color = "blue")
                        ax.add_patch(rect)
                self.ind.letharguslist = letharguslist
                
                """
                self.ind.fq_finallethargusperiod = self.detectlethargus(ind.fq_oescreendmatrix)
                if len(ind.fq_finallethargusperiod) > 0:
                    ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
                    ltstart = ind.fq_finallethargusperiod[0]
                    ltend = ind.fq_finallethargusperiod[1]        
                    ind.calcfoqlethargusmeasures(ltstart, ltend)
                """
        return fig, ax

    def plotcandidates(self, _ax):
        for oei in self.ind.fq_oescreendmatrix:
            #rect = plt.Rectangle((oei[0],0), oei[1]-oei[0],1,alpha = 0.2,hatch="/", color="blue")
            #rect = plt.Rectangle((oei[0],0), oei[1]-oei[0],1,alpha = 0.2,hatch="/",fill=False)
            rect = plt.Rectangle((oei[0],0), oei[1]-oei[0],0.5,alpha = 0.2, color = "gray")
            _ax.add_patch(rect)
        #return _ax

    def plotlethargus(self, _ax, _finallethargusperiod):
        #ax = _fig.get_axes()[0]
        #rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],1,alpha = 0.2,hatch="/", color="blue")
        #rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],0.05, alpha = 0.2, color="blue")
        rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],1, alpha = 0.2, color="gray")
        _ax.add_patch(rect)
        lthresh = 0.05
        _ax.axhline(y = lthresh, color = "gray", linewidth = 0.5)
        
        return _ax





        
    #################################################################
    #instead of fraction of q, use area rate(activity)
    #to determine onset/exit of lethargus
    def detectlethargus1stscreen_ar(self):
        #self.normalizedarea is return of def normalizearea(self):
        #20min run med
        runmedactivity = self.ind.normalizedarea.rolling(window=int(1200/self.interval),
                                    center= True).median()
        acthreash = 0.75
        self.ind.ar_qbooleanvec = runmedactivity < acthreash
        arbooleandiff = self.ind.ar_qbooleanvec.astype(int).diff()
        self.ind.ar_onsetcandidates = np.where(arbooleandiff == 1)[0]
        self.ind.ar_exitcandidates = np.where(arbooleandiff == -1)[0]
        


    #################################################################
    # Hayashi lab critera
    # quiescence onset time was defined as a time point 
    # after which the fractional quiescence remains > 0.05 for at least 60 minutes. 
    # >0.2, 1.5 hours
    # The quiescence exit time was defined as the time point 
    # after which the fractional quiescence had reached < 0.1. (actually 0.05)
    def detectlethargus1stscreen_fq(self, **kwargs):
    #def detectlethargus(_qbooleanvec, onsetcandidates, exitcandidates):
        #foq threashold to determine lethargus
        self.ind.fq_qbooleanvec = self.ind.foq > self.foqthreshod
        #convert series...
        #qbooleanvec = qbooleanvec.ix[:,0]
        #print(sum(qbooleanvec))
        qbooleandiff = self.ind.fq_qbooleanvec.astype(int).diff()
        self.ind.fq_onsetcandidates = np.where(qbooleandiff == 1)[0]
        self.ind.fq_exitcandidates = np.where(qbooleandiff == -1)[0]

    #screen by duration length        
    #foq based or area rate based are defined by input arg _onsetcandidates etc.
    def detectlethargus2ndscreen(self, _qbooleanvec, _onsetcandidates,
                                 _exitcandidates, **kwargs):
        #1h (3600 sec) need  lethargus duration, pre/post
        #20180123 changed to longer than 2hrs
        #continuouslength = int(2*60*60/self.interval)
        #20180223 may be 1.5?
        continuouslength = int(self.minduration*60*60/self.interval)
        #pre and post requirments could be adjustable.
        prepostmargin = int(1*60*60/self.interval)
        if "prepostmargin" in kwargs:
            prepostmarginhr = kwargs["prepostmargin"]
            prepostmargin  = int(prepostmarginhr*60*60/self.interval)

        temponsetcandidates=[]
        for oi in _onsetcandidates:
            #print("oi in onsetcandidates "+str(oi))
            if oi < prepostmargin :
                amessage = "pre-period imaging duraion is short"
                print(amessage)
                self.ind.foqdetectionlog.append(amessage)
                continue
            elif len(_qbooleanvec)- oi < prepostmargin:
                amessage = "post-period imaging duraion is short"
                print(amessage)
                self.ind.foqdetectionlog.append(amessage)
                continue
            sumofq = sum(_qbooleanvec[oi:(oi+continuouslength)])
            print("sumofq " + str(sumofq))
            if sumofq == continuouslength:
                amessage =  "suit the criteria"
                print(amessage)
                self.ind.foqdetectionlog.append(amessage)
                temponsetcandidates.append(oi)                
        print("temponsetcandidates " + str(temponsetcandidates))
        
        #exit when goes under threshold
        tempexitcandidates = _exitcandidates
        
        print("tempexitcandidates " +str(tempexitcandidates))
        
        onsetexitmatrix = []
        for i in range(len(temponsetcandidates)):
            #open is no foud exit yet
            openperiodfrag = True
            tempstart = temponsetcandidates[i]
            print("slice "+str(tempstart))
            if len(onsetexitmatrix) !=0:
                for j in range(len(onsetexitmatrix)):
                    if onsetexitmatrix[j][0] < tempstart < onsetexitmatrix[j][1]:
                        openperiodfrag = False
                        break
                
            if openperiodfrag:
                for j in range(len(tempexitcandidates)):
                    if tempexitcandidates[j] > tempstart:
                        onsetexitmatrix.append([tempstart, tempexitcandidates[j]])
                        openperiodfrag = False
                        break
        
        print("onsetexitmatrix "+str(onsetexitmatrix))
        return onsetexitmatrix

    #foq based or area rate based are defined by input arg _oescreendmatrix.    
    #self.ind.fq_oescreendmatrix
    def detectlethargus(self, _oescreendmatrix):
        # if there is one, [[aaa, bbb]], two [[aaa,bbb],[ccc,ddd]]
        finallethargusperiods=[]
        if len(self.ind.fq_oescreendmatrix) ==0:
            print("This sample desnt much lethargus criteria")
            #dataname = thedate + "_" + genotype+ "_" + str(samplenum +1) +  "_" +str(interval)
            #fig.savefig(dataname+"foqandlethargus.png",figsize=(8,2),dpi=100)
            #continue
            return finallethargusperiods
        else:
            #if ther are multiple lethargus state, choose by human eye?
            #if there are more than 2, need to check interval between them
            if len(self.ind.fq_oescreendmatrix)>1:
                #.ginput return [()]
                """
                print("please click lethargus period")
                ax.annotate("please click lethargus period",\
                            xy=(0.05, 0.8),\
                            xycoords='axes fraction', fontsize=8,\
                            horizontalalignment='left', verticalalignment='bottom')
                clickpos = plt.ginput(n=1)[0]
                print(clickpos)
                for oei in self.ind.fq_oescreendmatrix:
                    if oei[0] < clickpos[0] < oei[1]:
                        finallethargusperiod=oei.copy()
                        rect = plt.Rectangle((finallethargusperiod[0],0), 
                                    finallethargusperiod[1]-finallethargusperiod[0],
                                    1,alpha = 0.2,hatch="/", color="blue")
                        ax.add_patch(rect)
                if finallethargusperiod == []:
                    print("you didnt choose any candidate")
                    ax.annotate("you didnt choose any candidate",\
                            xy=(0.05, 0.7),\
                            xycoords='axes fraction', fontsize=8,\
                            horizontalalignment='left', 
                            verticalalignment='bottom', color = "red")
                """    
                
                # here is some filter to check if they are actual lethargus.
                # interval, mean foq? etc...
                #1st, filter by rawdata.rolling > threshold 0.01 rate
                # -> seems not work well l3 lethargus tend have high rate?
                """
                rawdata = self.ind.calcrawarea(60)
                rmrawdata = rawdata.rolling(window = 60).median()
                threshold = 0.01
                for lt in self.ind.fq_oescreendmatrix:
                    print(lt)
                    #1min running median over the thredhold
                    sumoverth = sum((rmrawdata[lt[0]: lt[1]]) > threshold)
                    #if overthe threshold time is more than 80 %
                    if (sumoverth/(lt[1]-lt[0])) > 0.8:
                        print(sumoverth/(lt[1]-lt[0]))
                        print("This has high motion rate")
                    else:
                        print("This one seems lethargus")
                        finallethargusperiods.append(lt)
                """
                        
                
                
            else:
                finallethargusperiods=self.ind.fq_oescreendmatrix.copy()
            
            
            print("finallethargusperiods "+str(finallethargusperiods))
            #rect = plt.Rectangle((finallethargusperiod[0],0), finallethargusperiod[1]-finallethargusperiod[0],1,alpha = 0.2,hatch="/", color="blue")
            #ax.add_patch(rect)
        return finallethargusperiods
    
        #dataname = thedate + "_" + genotype+ "_" + str(samplenum+1) +  "_" +str(interval)
        #fig.savefig(dataname+"foqandlethargus.png",figsize=(8,2),dpi=100)
        

    #if several closer each other, discard less totalq one?
    def checkdistance(self, letharguslist):
        ll = letharguslist
        self.problematicdistance = False
        for i in range(len(ll)-1):
            currentend = ll[i].end
            postend = ll[i+1].end
            distance = postend- currentend
            if distance < self.mininterval*60*60/self.ind.interval:
                ll[i].distanceproblem = True
                #ll[i+1].distanceproblem = True
                self.problematicdistance = True
                
    #if several closer each other, discard less totalq one?
    def filterbytotalq(self, letharguslist):
        ll = letharguslist
        for i in range(len(ll)-1):
            print("len(ll) " + str(len(ll)) +" i "+str(i))
            if len(ll) != i:
                if ll[i].distanceproblem:
                    print("ll start at "+str(ll[i].start/(60*60/interval))+" has problem")
                    #ll[i].calcmeasurements()
                    #ll[i+1].calcmeasurements()
                    if ll[i].totalq > ll[i+1].totalq:
                        ll.remove(ll[i+1])
                    else:
                        ll.remove(ll[i])
                    ll[i].distanceproblem = False
                    #ll[i+1].distanceproblem = False
        self.checkdistance(ll)
        #181225 debug. screen2 #244_9 cause error
        """
        ld = Lethargusdetector(0.2,1.5,10)
        ind = Individual(thedate, groupname, expnum, interval, 9, 
                            tempdf.loc[:,tempdf.columns[8]])
        ld.setadata(ind)        
        datalabel = "_".join([ind.date,ind.groupname, str(ind.expnum), 
                                  str(ind.samplenum)])
            print("---------------------")
            print(datalabel)
            
            
            #low level process
            ind.qaboolean = ind.calcqa(ind.rawdata)
            ind.foq = ind.calcfoq(ind.rawdata)
            ind.normalizedarea = ind.normalizearea()
            ind.arearate = ind.calcarearate()
            
            ld.setadata(ind)
            fig, ax = ld.processdata()         
            
        Lethargusdetector.filterbytotalq = filterbytotalq
        """

        
    def prepfig(self, **kwargs):
        figsize  = (8,2)
        if "figsize" in kwargs:
            figsize = kwargs["figsize"]
        
        fig = plt.figure(figsize = figsize)
        return self.prepplot(fig, **kwargs)

    def prepplot(self, fig, **kwargs):        
        labeloff = False
        overlay = None
        xlim = (0, len(self.ind.rawdata))
        ylim = (0, 1.1)
        figsize  = (8,2)
        if "labeloff" in kwargs:
            labeloff = kwargs["labeloff"]
        if "overlay" in kwargs:
            overlay = kwargs["overlay"]
        if "xlim" in kwargs:
            xlim = kwargs["xlim"]
        if "ylim" in kwargs:
            ylim = kwargs["ylim"]
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.spines['right'].set_visible(False)
        #ax.spines['top'].set_visible(False)
        ax2 = ax.twiny()
        ax2.xaxis.tick_top()
        #originalxtick = ax.get_xticklabels()
        hrtick = np.arange(xlim[1]*self.ind.interval/60/60).astype(int)
        ax.set_xticks(hrtick*60*60/self.ind.interval)
        ax.set_xticklabels([]) 
        if labeloff:
            ax.set_yticklabels([]) 

        label = "_".join([self.ind.date, self.ind.groupname,
                          str(self.ind.expnum),
                          str(self.ind.samplenum)])
        if overlay is not None:
            if type(overlay[0]) != list:
                ollist = [overlay]
            else:
                ollist = overlay
            for aol in ollist:
                ostart = aol[0]*60*60/ind.interval
                oend = aol[1]*60*60/ind.interval
                rect = plt.Rectangle((ostart,0), oend-ostart,1, 
                                 alpha = 0.2, color="blue")
                ax.add_patch(rect)
                    
        ax.annotate(label,\
        xy=(0.01, 0.9),\
         xycoords='axes fraction', fontsize=8,\
         horizontalalignment='left', verticalalignment='bottom')

        #ax.set_xticklabels([int(x) for x in ax.get_xticks()])
        fig.tight_layout()
        #ax.xaxis.tick_top()
        ax.set_xticklabels(hrtick)
        return fig, ax, ax2
    
"""
foqthreshold = 0.2
foqthreshold = 0.05
minduration = 1.5
mininterval = 10

ld = Lethargusdetector(foqthreshold,minduration,mininterval)
ld.setadata(sg.fullindlist[23])
ld.processdata()

fig, ax = ld.prepfig()
ax.plot(sg.fullindlist[0].calcrawarea(30))
ax.clear()
ld.setadata(sg.fullindlist[21])
fig, ax = ld.prepplot(fig)
ax.plot(sg.fullindlist[21].foq)
    
"""

#Individual.calcrawarea =calcrawarea

class Individual:
    
    def __init__(self, date, groupname, expnum, interval, samplenum, rawdata):
        #data obtained from .xls file name
        self.date = date
        self.groupname = groupname
        self.expnum = expnum
        self.interval = interval
        
        #one .xls file contains several individuals
        self.samplenum = samplenum
        
        #data        
        self.rawdata = rawdata
        
        #calcurated from rawdata
        self.qaboolean = None
        self.foq = None #10min running average
        self.curmaxarea = None
        self.normalizedarea = None
        self.arearate = None #20min running median
               
        #foq depend lethargus
        self.foqdetectionlog = []
        #if there are multiple candidate, human correction required
        self.fq_onsetcandidates = None
        self.fq_exitcandidates = None
        self.fq_qbooleanvec = None
        self.fq_oescreendmatrix = None
        self.fq_finallethargusperiod = None
        
        #20180226 lethargusdetector implemented
        self.letharguslist = None
        #end to end
        self.interlethargus = []
        
        
        #self.fq_finallethargusexit = None
        #measurements
        self.meanquiescent = None
        self.meanquiescentout = None
        self.totalq = None
        self.numberofbout = None
        self.qmean = None 
        self.amean = None
        self.qmedian = None 
        self.amedian = None 
        self.fq_duration = None
        self.qfreq = None
        self.qaboutdataframe = None
        self.lethargusperiodqadf = None

        #area rate depend lethargus
        self.ar_onsetcandidates = None
        self.ar_exitcandidates = None
        self.ar_qbooleanvec = None
        self.ar_oescreendmatrix = None
        self.ar_finallethargusperiod = None
        #self.ar_finallethargusexit = None
       
    
    #################################################################
    # calc (determine ) quiescent or active
    # _dataseries is the rawdata which consist of pixel number above threshold
    def calcqa(self, _dataseries, **kwargs):
        pixthreash = 1
        if "pixthreash" in kwargs:
            pixthreash = kwargs["pixthreash"]
        #false true array < 1 are q.
        _qaboolean = _dataseries < pixthreash    
        
        return _qaboolean
    
    #################################################################
    # calculate fraction of q
    # _dataseries is the rawdata which consist of pixel number above threshold
    def calcfoq(self, _dataseries, **kwargs):
        _qaboolean = self.calcqa(_dataseries, **kwargs)
    
        #10min 600sec window
        #use reset_index(drop=True) to eliminate nan
        _foq=_qaboolean.rolling(window = int(600/self.interval), center=False).mean()
        #_foq=self.qaboolean.rolling(window = int(600/self.interval), center=False).mean()
        _foq=_foq.dropna().reset_index(drop=True)
        return _foq

    #################################################################
    #area normalized with max area. not curmax
    def calcrawarea(self, window):
        medianraw = self.rawdata.rolling(window=int(window/self.interval),
                            center= True).median()

        maxraw = np.nanmax(medianraw)
        rawdataproportion = self.rawdata/maxraw
        return rawdataproportion

    #normalize the area by accumulation max of runmed area
    def calccurmaxarea(self):
        windowsize = int(60/self.interval)
        rawrunmed = self.rawdata.rolling(window=windowsize, 
                                         center= False).median()            
        #data = rawrunmed.values.flatten()       
        data = rawrunmed.dropna().reset_index(drop=True)       
        tempmax = 0
        returnval = np.zeros(len(self.rawdata))
        for i in range(len(data)):
            if data[i] > tempmax:
                tempmax = data[i]
            returnval[windowsize -1 + i] = tempmax
        for i in range(windowsize-1):
            returnval[i] = returnval[windowsize-1]
            
        return returnval

    def normalizearea(self):
        self.curmaxarea = self.calccurmaxarea()
        correctedrunmedarea = self.rawdata/self.curmaxarea
        return correctedrunmedarea
    
    def calcarearate(self):
        return self.normalizedarea.rolling(window=int(1200/self.interval), center= True).median()
        
        
    #################################################################
    def preparefig(self):
        fig = plt.figure(figsize = (8,2))
        #ax.cla()
        #delieat top/right axis
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
        plt.ylim(0,1.1)
        plt.xlim(0, len(self.rawdata))
        ax = fig.add_subplot(1,1,1)
        #fig.tight_layout()

        return fig
        
    def plotnormalizedarea(self, _ax, _window):
        _ax.plot(self.normalizedarea.rolling(window=int(_window/self.interval), center= True).median(), color = "gray", linewidth =0.5)
        return _ax
        
    def plotfoq(self, _ax):
        _ax.plot(self.foq, color = "black")
        #foq threashold to determine lethargus
        #lthresh = 0.05
        #ax.axhline(y = lthresh, color = "gray")
        return _ax

    def plotlowlevel(self, _ax):
        #_ax.plot(self.normalizedarea.rolling(window=int(600/interval), center= True).median(), color = "gray", linewidth =0.5)
        _ax.plot(self.arearate, color = "gray", linewidth =0.5)
        
        _ax.plot(self.foq, color = "black")
        #foq threashold to determine lethargus
        #lthresh = 0.05
        #ax.axhline(y = lthresh, color = "gray")
        return _ax
    
    def plotcandidates(self, _ax, _oescreendmatrix):
        for oei in _oescreendmatrix:
            #rect = plt.Rectangle((oei[0],0), oei[1]-oei[0],1,alpha = 0.2,hatch="/", color="blue")
            #rect = plt.Rectangle((oei[0],0), oei[1]-oei[0],1,alpha = 0.2,hatch="/",fill=False)
            rect = plt.Rectangle((oei[0],0), oei[1]-oei[0],0.5,alpha = 0.2, color = "gray")
            _ax.add_patch(rect)
        return _ax

    def plotlethargus(self, _ax, _finallethargusperiod):
        #ax = _fig.get_axes()[0]
        #rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],1,alpha = 0.2,hatch="/", color="blue")
        #rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],0.05, alpha = 0.2, color="blue")
        rect = plt.Rectangle((_finallethargusperiod[0],0), _finallethargusperiod[1]-_finallethargusperiod[0],1, alpha = 0.2, color="gray")
        _ax.add_patch(rect)
        lthresh = 0.05
        _ax.axhline(y = lthresh, color = "gray", linewidth = 0.5)
        
        return _ax
        

    def saveafig(self, _fig, opt = None):
        figexp = ""
        if opt is None:
            figexp = "foq"
        else:
            figexp = opt
        dataname = "_".join([self.date, self.groupname, str(self.expnum), str(self.samplenum),"_"])
        #dataname = thedate + "_" + genotype+ "_" + str(samplenum+1) +  "_" +str(interval)
        #_fig.savefig(dataname+"rawfoqar.png",figsize=(8,2),dpi=100)
        #need to add directory
        #_fig.savefig(targetdir +"/"+dataname+"rawfoqar.png",figsize=(8,2),dpi=100)
        _fig.set_size_inches(8, 2)#240220 py311 need this?
        #_fig.savefig(targetdir +"/"+dataname+figexp +".png",figsize=(8,2),dpi=100)
        _fig.savefig(targetdir +"/"+dataname+figexp +".png",dpi=100)
        
    #180226 not yet handling multiple lethargus
    def savefoqincsvfile(self):        
        #save the lethargus period and foq with the date_groupname_expnum_# name
        dataname = "_".join([self.date, self.groupname, str(self.expnum), 
                             str(self.samplenum),"_"])
        with open(targetdir +"/"+dataname+"foq"+".csv","w", newline='') as file:
            #lineterminator=os.linesep need?
            writer = csv.writer(file)
            writer.writerow(["_".join([self.date, self.groupname, str(self.expnum),
                                      str(self.samplenum), str(self.interval)])])
            if self.fq_finallethargusperiod is not None and \
                self.fq_finallethargusperiod != []:
                writer.writerow([self.fq_finallethargusperiod[0]])
                writer.writerow([self.fq_finallethargusperiod[1]])
            else:
                writer.writerow([None])
                writer.writerow([None])                    
            #writer.writerows(zip(self.rawdata))
            writer.writerows(zip(self.foq))
            


class Samplegroup:   
    def __init__(self, _groupname):
        #date_groupname_expnum format. may be need change to only groupname later?
        self.groupname = _groupname
        #include all sample. some are not have lethargus period
        self.fullindlist = []
        #only samples that have one lethargus period
        self.indlist = []
        self.summarydf = None
        self.fq_maxduration = None
        self.ltfoqmatrix = None
        self.ltfoqmatrixtail = None
        
        
    def processagroup(self, thedate, groupname, expnum, interval, dataframe, **kwargs):
        i=0
        if "threshold" in kwargs:
            threshold = kwargs["threshold"]
        else:
            threshold = 0.05
            
        #foqthreshold = 0.2
        #foqthreshold = 0.05
        minduration = 1.5
        #190305 test to detect rem30
        #minduration = 0.5
        mininterval = 10
        ld = Lethargusdetector(threshold,minduration,mininterval)
        #for i in range (len(areadf.columns)):
        for i in range (len(dataframe.columns)):
            #i=6
            #samplenum start from 1
            ind = Individual(thedate, groupname, expnum, interval, i+1, 
                             dataframe.loc[:,dataframe.columns[i]])
            datalabel = "_".join([ind.date,ind.groupname, str(ind.expnum), 
                                  str(ind.samplenum)])
            print("---------------------")
            print(datalabel)
            
            
            #low level process
            if "pixthreash" in kwargs:
                qapixthreash = kwargs["pixthreash"]
            
            ind.qaboolean = ind.calcqa(ind.rawdata,pixthreash = qapixthreash)
            ind.foq = ind.calcfoq(ind.rawdata,pixthreash = qapixthreash)
            ind.normalizedarea = ind.normalizearea()
            ind.arearate = ind.calcarearate()
            
            ld.setadata(ind)
            fig, ax = ld.processdata()

            """
            #lethargus detectoin by foq
            #find slice where foq over/go under the threshold
            ind.detectlethargus1stscreen_fq(threshold = threshold)
            fig = ind.preparefig()
            ax = ind.plotlowlevel(fig.get_axes()[0])
            if len(ind.fq_onsetcandidates) >0:
                #screen periods that have 1h pre/post lethargus and over 1h duration
                ind.fq_oescreendmatrix = ind.detectlethargus2ndscreen(ind.fq_qbooleanvec, ind.fq_onsetcandidates, ind.fq_exitcandidates)
                if len(ind.fq_oescreendmatrix) > 0:
                    ax = ind.plotcandidates(ax, ind.fq_oescreendmatrix)
                    #even after above, still may exist multiple. use human eyes.
                    ind.fq_finallethargusperiod = ind.detectlethargus(ind.fq_oescreendmatrix)
                    if len(ind.fq_finallethargusperiod) > 0:
                        ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
                        ltstart = ind.fq_finallethargusperiod[0]
                        ltend = ind.fq_finallethargusperiod[1]        
                        ind.calcfoqlethargusmeasures(ltstart, ltend)
                        #ind.savefoqincsvfile()
                        #sg.appendanindividual(ind)
            #fig.tight_layout()
            """
            ind.savefoqincsvfile()

            
            self.appendanindividual(ind)
                
            ax.annotate(datalabel,\
                        xy=(0.01, 0.9),\
                         xycoords='axes fraction', fontsize=8,\
                         horizontalalignment='left', verticalalignment='bottom')
                    
            ind.saveafig(fig)
            plt.close(fig)
            


    def appendanindividual(self, _ind):
        self.fullindlist.append(_ind)
        #ind.letharguslist contains multiple lethargus
        #if(_ind.fq_finallethargusperiod is not None and\
        #   _ind.fq_finallethargusperiod != []):
        if _ind.letharguslist is not None:
            self.indlist.append(_ind)
        
    #assign by samplenum
    def deleteanindividual(self, _samplenum):
        #del self.indlist[index-1]
        samplenumlist = [a.samplenum for a in self.indlist]
        del self.indlist[samplenumlist.index(_samplenum)]

    #assign by samplenum
    def deletealethargus(self, _samplenum, lethargusnum):
        #del self.indlist[index-1]
        samplenumlist = [a.samplenum for a in self.fullindlist]
        ind = self.fullindlist[samplenumlist.index(_samplenum)]
        del ind.letharguslist[lethargusnum]
        
    def makedf(self, **kwargs):
        df = pandas.DataFrame()       
        ltlist = self.indlist
        if "ltlist" in kwargs:
            ltlist = kwargs["ltlist"]
            
        for ind in ltlist:
            dataname = "_".join([ind.date, ind.groupname,
                                 str(ind.expnum), str(ind.samplenum)])
            print(dataname, ind.meanquiescent)
            if ind.letharguslist is None:
                print("No lethargus detected")
            else:                
                for lt in ind.letharguslist:
                    print(str(ind.letharguslist.index(lt)))
                    tempdf = pandas.DataFrame(
                            {#"name": dataname,
                             "date": ind.date,
                             "groupname": ind.groupname,
                             "expnum": str(ind.expnum),
                             "samplenum": str(ind.samplenum),
                             "ltnum": str(ind.letharguslist.index(lt)),
                             "fqlt_start": [lt.start],
                             "fqlt_end": [lt.end],
                             "fqlt_duration": [lt.fq_duration],
                             "totalq": [lt.totalq],
                             "meanquiescent": [lt.meanquiescent],
                             "meanquiescentout": [lt.meanquiescentout],
                             "numberofbout": [lt.numberofbout],
                             "qmean": [lt.qmean],
                             "amean": [lt.amean],
                             "qmedian": [lt.qmedian],
                             "amedian":[lt.amedian],
                             "qfreq": [lt.qfreq]
                             })
                    #df = df.append(tempdf)
                    df = pandas.concat([df,tempdf])
            #df = df.append(pandas.Series(ind.meanquiescent),ignore_index=True)
        #df = pandas.DataFrame({"name" : namelist, measurement: datalist})
        #df = df.append(pandas.Series([ind.meanquiescent], index = [ind.groupname]))
        #set the order of colname.
        if not df.empty:
            #df = df[["name",
            df = df[["date","groupname","expnum","samplenum","ltnum",
                     "fqlt_start","fqlt_end","fqlt_duration","totalq",
                     "meanquiescent","meanquiescentout","numberofbout",
                     "qmean","amean","qmedian","amedian","qfreq"]]
            self.summarydf = df
        return df
    
    
    def calcmaxduration(self, **kwargs):
        self.fq_maxduration = 0
        #for ind in self.indlist:
        ordinalnum = 0
        if "ordinalnum" in kwargs:
            ordinalnum = kwargs["ordinalnum"]
        self.numofadequatesample = 0
        for ind in self.fullindlist:
            #ltstart = ind.fq_finallethargusperiod[0]
            #ltend = ind.fq_finallethargusperiod[1]
            if ind.letharguslist is not None:
                if len(ind.letharguslist) > ordinalnum:
                    self.numofadequatesample = self.numofadequatesample+1
                    ltstart = ind.letharguslist[ordinalnum].start
                    ltend = ind.letharguslist[ordinalnum].end
                    if ltend-ltstart+1 > self.fq_maxduration:
                        self.fq_maxduration = ltend-ltstart+1

    #adataseries could be foq, arearate etc. alignhead True or tail False
    #def makeltmatrix(self, adataseries, start, end, alignhead):
    #def makeltmatrix(self, alignhead, **kwargs):
    #**kwargs could have which lethargus uses (L3 or L4..)
    def makeltmatrix(self, **kwargs):
        align = "head"
        if "align" in kwargs:
            align = kwargs["align"]
        ordinalnum = 0
        if "ordinalnum" in kwargs:
            ordinalnum = kwargs["ordinalnum"]
        #namelist = []
        if self.fq_maxduration is None:
            self.calcmaxduration(**kwargs)
                
       
        #emptymatrix = np.zeros(shape = (len(self.indlist), 
        emptymatrix = np.zeros(shape = ( sg.numofadequatesample, 
                            self.fq_maxduration+int(60*60/self.fullindlist[0].interval)*2))
        emptymatrix.fill(np.nan)
        ltfoqmatrix = emptymatrix.copy()
        
        j=0
        for ind in self.fullindlist:
            if (ind.letharguslist is not None) and (len(ind.letharguslist) > 0):
                ltstart = ind.letharguslist[ordinalnum].start
                ltend = ind.letharguslist[ordinalnum].end
                #margin 1hr
                swithmargin = int(ltstart-60*60/ind.interval)
                ewithmargin = int(ltend+60*60/ind.interval)
                #in case that pre/post lethargus not enough margin
                emptymat = np.empty((1,ewithmargin-swithmargin+1))
                emptymat.fill(np.nan)
                ltfoq = emptymat[0]
                #200330 if swithmargin is smallter than 0, it return empty series
                foqextracted = ind.foq[swithmargin:ewithmargin]
                ltfoq[:len(foqextracted)] = foqextracted       
                #if not alignhead:
                if align == "tail":
                    #ltfoq = ltfoq.reverse()
                    #ltfoq = ltfoq.iloc[::-1]
                    ltfoq = ltfoq[::-1]
                ltfoqmatrix[j][0:len(ltfoq)] = np.array(ltfoq).flatten()
                #something to make sub list having only adequate samples
                #self.indlist.
                j = j+1
                
        
        if align == "head":
            self.ltfoqmatrix = ltfoqmatrix
        else:            
            self.ltfoqmatrixtail = np.fliplr(ltfoqmatrix)
        return ltfoqmatrix

    def makeamatrix(self, datatype, **kwargs):
        matrixrow = len(self.fullindlist)
        matrixcol = len(self.fullindlist[0].rawdata)
        xlim = (0, matrixcol)
        indlist = self.fullindlist
        if "xlim" in kwargs:
            xlim = kwargs["xlim"]
            matrixcol = xlim[1] - xlim[0]
        else:
            #matrixcol = matrixcol-1
            pass
        if "samplenumlist" in kwargs:
            samplenumlist = kwargs["samplenumlist"]
            matrixrow = len(samplenumlist)
            indlist = []
            #for 
        autoflag = False
        #automatically eliminate high foq sample since it is no sample
        if "autoelim" in kwargs:
            autoflag = kwargs["autoelim"]
            
            
        
        ematrix = np.zeros(shape = (matrixrow, matrixcol))
        ematrix.fill(np.nan)
        
        for i, ind in enumerate(indlist):
            if datatype == "foq":
                data = ind.foq[xlim[0]:(xlim[1]-1)]
                        
                    
            elif datatype == "area":
                tempraw = ind.calcrawarea(60)
                data = tempraw[xlim[0]:(xlim[1]-1)]
        
            #ematrix[i][xlim[0]:(xlim[1]-1)] = data
            ematrix[i][0:len(data)] = data
            if autoflag:
                if np.mean(data) > 0.9:#arbitrary value....
                    ematrix[i][0:len(data)] = np.nan
            
        rmatrix = ematrix
        
        return rmatrix
        
    
    def savesummarydf(self):
        dataname = "_".join([self.groupname])
        interval = self.fullindlist[0].interval
        self.summarydf.to_csv(targetdir +"/"+dataname+"_"+str(int(interval))+"_summary"+".csv", index=False )
        """
        with open(targetdir +"/"+dataname+"_summary"+".csv","w", newline='') as file:
            #lineterminator=os.linesep need?
            writer = csv.writer(file)
            writer.writerow(["_".join([self.date, self.groupname, str(self.expnum),
                                      str(self.samplenum), str(self.interval)])])
            writer.writerow([self.fq_finallethargusperiod[0]])
            writer.writerow([self.fq_finallethargusperiod[1]])
            #writer.writerows(zip(self.rawdata))
            writer.writerows(zip(self.foq))
        """    
    #_amatrix is self.ltfoqmatrix or ltfoqmatrixtail
    def makealignfigure(self, _amatrix, alignhead, representtype):
        fig = plt.figure(figsize = (8,4))
        ax = fig.add_subplot(1,1,1)
        ax.set_ylim(-0.1,1)
        ax.set_xlim(0, _amatrix.shape[1])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        interval = self.fullindlist[0].interval
        hrtick = np.arange(_amatrix.shape[1]*interval/60/60).astype(int)
        ax.set_xticks(hrtick*60*60/interval)
        ax.set_xticklabels(hrtick-1)
        if not alignhead:
             shift = _amatrix.shape[1]-(hrtick*60*60/interval)[-1]
             ax.set_xticks(hrtick*60*60/interval + shift)
             ax.set_xticklabels(hrtick-max(hrtick)+1)
        for adata in _amatrix:
            ax.plot(adata, linestyle = ":", linewidth =0.5, color = "gray")
            
        if representtype == "mean":
            mean = np.nanmean(_amatrix, axis = 0)
            sd = np.nanstd(_amatrix, axis = 0)
            ax.plot(mean+sd, linestyle = "--", linewidth = 1, color = "gray")
            ax.plot(mean-sd, linestyle = "--", linewidth = 1, color = "gray")
            ax.plot(mean, linestyle = "-", linewidth = 1, color = "black")
        elif representtype == "median":
            median = np.nanmedian(_amatrix, axis = 0)
            sd = np.nanstd(_amatrix, axis = 0)
            ax.plot(median, linestyle = "-", linewidth = 1, color = "black")
        
        fig.tight_layout()   
        ax.annotate(self.groupname+" ("+str(_amatrix.shape[0])+")",\
                xy=(0.01, 0.9),\
                 xycoords='axes fraction', fontsize=8,\
                 horizontalalignment='left', verticalalignment='bottom')
        return fig


    def makeallfigure(self):
        fig = self.prepallfigureformat([0,len(self.indlist[0].rawdata)],
                                        [0,1],self.indlist)
        #ax = ind.plotfoq(ax)
        #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        axlist = fig.get_axes()
        i=0
        for ax in axlist:
            ind = self.indlist[i]
            ax = ind.plotfoq(ax)
            ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
            i = i + 1
        return fig
    
    def makeallfigurefoq(self, _indlist, **kwargs):
        window = 60
        xlim = [0,len(_indlist[0].rawdata)]
        ylim = [0,1.1]
        if "xlim" in kwargs:
            xlim = kwargs["xlim"]
        if "ylim" in kwargs:
            ylim = kwargs["ylim"]
            
        fig = self.prepallfigureformat(xlim,ylim,_indlist,**kwargs)

        letopt = "line"
        #ax = ind.plotfoq(ax)
        #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        axlist = fig.get_axes()
        i=0
        for ax in axlist:
            ind = _indlist[i]
            ax = ind.plotfoq(ax)
            #ax = ind.plotnormalizedarea(ax,window)
            #ax.plot(ind.normalizedarea.rolling(window=int(window/ind.interval), center= True).median(), color = "black", linewidth =0.5)
            ##ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
            if "lethargus" in kwargs:
                letopt = kwargs["lethargus"]
            if letopt != "off":
                if ind.letharguslist is not None:
                #if len(ind.letharguslist) >0:
                    for lt in ind.letharguslist:
                        ax.axvline(x= lt.start, color= "black")
                        ax.axvline(x= lt.end, color= "black", linestyle = "-")

            i = i + 1
        return fig        

    def makeallfigurewithinitialq(self, _indlist):
        window = 60
        fig = self.prepallfigureformat([0,len(_indlist[0].rawdata)],
                                        [0,1.1], _indlist)
        #ax = ind.plotfoq(ax)
        #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        axlist = fig.get_axes()
        i=0
        for ax in axlist:
            ind = _indlist[i]
            ax = ind.plotfoq(ax)
            #ax = ind.plotnormalizedarea(ax,window)
            #ax.plot(ind.normalizedarea.rolling(window=int(window/ind.interval), center= True).median(), color = "black", linewidth =0.5)
            ##ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
            
            #initical q
            ax.axvline(x= ind.fq_onsetcandidates[0], color= "black")
            initialafterlethargus = ind.fq_onsetcandidates[ind.fq_onsetcandidates > ind.fq_finallethargusperiod[1]][0]
            ax.axvline(x= ind.fq_finallethargusperiod[1], color= "black", linestyle = "--")
            ax.axvline(x= initialafterlethargus, color= "black")
            i = i + 1
        return fig        
    
    def calcinitialqlatency(self):
        print("initial reactivate")
        for ind in self.indlist:
            initiallatency = ind.fq_onsetcandidates[0]
            
            initialafterlethargus = ind.fq_onsetcandidates[ind.fq_onsetcandidates > ind.fq_finallethargusperiod[1]][0]
            print(str(initiallatency*ind.interval) +" "+str((initialafterlethargus-ind.fq_finallethargusperiod[1])*ind.interval))
        #return
    
    #_xlim and _ylim two vlues list, _figsize is taple
    #def prepallfigureformat(self,_xlim, _ylim, _figsize):
    #_xlim and _ylim two vlues list, _indlist could be self.indlist or fullindlist
    #kwargs: col row
    def prepallfigureformat(self,_xlim, _ylim, _indlist, width = None,
                            **kwargs):
        #fig = plt.figure(figsize = _figsize)
        if width is None:
            width = 8
        col = 1
        row = len(_indlist)
        labeloff = False
        overlay = None
        if "col" in kwargs:
            col = kwargs["col"]
        if "row" in kwargs:
            row = kwargs["row"]
        if "labeloff" in kwargs:
            labeloff = kwargs["labeloff"]
        if "overlay" in kwargs:
            overlay = kwargs["overlay"]
        if "figsize" in kwargs:
            figsize = kwargs["figsize"]
        else:
            figsize  = (width,row)
        
            
        fig = plt.figure(figsize = figsize)
        #ax.cla()
        #delieat top/right axis
        #plt.gca().spines['right'].set_visible(False)
        #plt.gca().spines['top'].set_visible(False)
        #plt.gca().yaxis.set_ticks_position('left')
        #plt.gca().xaxis.set_ticks_position('bottom')
        
        #for i in range(len(self.indlist)):
        for i in range(len(_indlist)):
            #ind = self.indlist[i]
            ind = _indlist[i]
            ax = fig.add_subplot(row,col,i+1)
            #ax.set_ylim(0,1)
            #ax.set_xlim(0, len(ind.rawdata))
            ax.set_xlim(_xlim[0], _xlim[1])
            ax.set_ylim(_ylim[0],_ylim[1])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            #hrtick = np.arange(len(ind.rawdata)*ind.interval/60/60).astype(int)
            hrtick = np.arange(_xlim[1]*ind.interval/60/60).astype(int)
            ax.set_xticks(hrtick*60*60/ind.interval)
            ax.set_xticklabels([]) 
            if labeloff:
                ax.set_yticklabels([]) 
            if overlay is not None:
                if type(overlay[0]) != list:
                    ollist = [overlay]
                else:
                    ollist = overlay
                for aol in ollist:
                    ostart = aol[0]*60*60/ind.interval
                    oend = aol[1]*60*60/ind.interval
                    rect = plt.Rectangle((ostart,0), oend-ostart,1, 
                                     alpha = 0.2, color="blue")
                    ax.add_patch(rect)
            #ax = ind.plotfoq(ax)
            #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
            
            fig.tight_layout()        
            label = "_".join([ind.date, ind.groupname, str(ind.expnum),
                              str(ind.samplenum)])
            ax.annotate(label,\
                    xy=(0.01, 0.9),\
                     xycoords='axes fraction', fontsize=8,\
                     horizontalalignment='left', verticalalignment='bottom')

        #ax.set_xticklabels([int(x) for x in ax.get_xticks()])
        fig.tight_layout()
        ax.set_xticklabels(hrtick)
        return fig
    
    def makeallfigureofarea(self, window, _indlist,**kwargs):
        fig = self.prepallfigureformat([0,len(_indlist[0].rawdata)],
                                        [0,1.1],_indlist,**kwargs)
        #ax = ind.plotfoq(ax)
        #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        axlist = fig.get_axes()
        i=0
        for ax in axlist:
            ind = _indlist[i]
            #ax = ind.plotnormalizedarea(ax,window)
            ax.plot(ind.normalizedarea.rolling(window=int(window/ind.interval), center= True).median(), color = "black", linewidth =0.5)
            #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
            i = i + 1
        return fig
    
    def makeallfigureofrawarea(self, window, _indlist, rollingtype,**kwargs):
        xlim = [0,len(_indlist[0].rawdata)]
        ylim = [0,1.1]
        if "xlim" in kwargs:
            xlim = kwargs["xlim"]
        if "ylim" in kwargs:
            ylim = kwargs["ylim"]
            
        fig = self.prepallfigureformat(xlim,ylim,_indlist,**kwargs)
        #ax = ind.plotfoq(ax)
        #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
        axlist = fig.get_axes()
        i=0
        for ax in axlist:
            ind = _indlist[i]
            #ax = ind.plotnormalizedarea(ax,window)
            rawdataproportion = ind.rawdata/np.nanmax(ind.rawdata.rolling(window=int(window/ind.interval), center= True).median())
            if rollingtype == "median":
                ax.plot(rawdataproportion.rolling(window=int(window/ind.interval), center= True).median(), color = "black", linewidth =0.5)
            elif rollingtype == "mean":
                ax.plot(rawdataproportion.rolling(window=int(window/ind.interval), center= True).mean(), color = "black", linewidth =0.5)
            #ax = ind.plotlethargus(ax, ind.fq_finallethargusperiod)
            i = i + 1
        return fig

    def makeallfigureofareacropped(self,_indlist,**kwargs):
        if self.fq_maxduration is None:
            self.calcmaxduration()
        maxduration = self.fq_maxduration
        if "maxduration" in kwargs:
             maxduration = kwargs["maxduration"]
        ordinalnum = 0
        if "ordinalnum" in kwargs:
            ordinalnum = kwargs["ordinalnum"]
        #plus minus 1 h
        margin = int(60*60/_indlist[0].interval)
        #foq lethargus correction. it is box car 10 min. so, eliminate 10min
        corcoef = int(60*10/_indlist[0].interval)
        #fig = self.prepallfigureformat([0,self.fq_maxduration + margin*2 - corcoef*2],
        #                                [0,1.2],(16,len(self.indlist)))
        #fig = self.prepallfigureformat([0,self.fq_maxduration + margin*2 - corcoef*2],
        #                                [0,1.2],(self.fq_maxduration/60/30*6,
        #                                        len(self.indlist)))
        xlim = [0,maxduration + margin*2 - corcoef*2]
        ylim = [0,1.2]
        fig = self.prepallfigureformat(xlim,
                                        ylim,_indlist,
                                        maxduration/60/30*6)
        axlist = fig.get_axes()
        i=0
        for ax in axlist:
            ind = _indlist[i]
            #if (ind.letharguslist is not None) and (len(ind.letharguslist) > 0):
            ltstart = ind.letharguslist[ordinalnum].start
            ltend = ind.letharguslist[ordinalnum].end
            #plus minus 1 h
            #margin = int(60*60/ind.interval)
            #foq lethargus correction. it is box car 10 min. so, eliminate 10min
            #corcoef = int(60*10/ind.interval)
            #plotstart = ltstart - margin + corcoef
            #plotend = ltend + margin - corcoef
            #200203 foq already trancated at the initial 10min,
            # so, no need to adjust at start? but >0.2 foq draw line near the 
            # down edge of area. so add corcoef
            #and not shift at the end of period.
            #because this is just for visualization of area accroding to lethargus
            plotstart = ltstart - margin + corcoef
            plotend = ltend + margin
            #ax = ind.plotnormalizedarea(ax,window)
            tempdata = ind.normalizedarea[plotstart:plotend]
            ax.plot(tempdata.values, linestyle = "-",
                    linewidth =0.5, color = "gray")
            
        
            ax.plot(tempdata.rolling(window=int(60/ind.interval),
                                     center= True).median().values,
                    linestyle = "-",
                    linewidth =1, color = "black")
            i = i + 1
            
        #x tick is now 0,1,2... so start it from -1 hr.
        shiftedtick = [int(a.get_text())-1 for a in ax.get_xticklabels()]
        ax.set_xticklabels(shiftedtick)
        xlimit = maxduration + margin*2 - corcoef*2
        (maxduration/60/30*6,
                                                len(_indlist))
        fig.set_size_inches(xlimit/60/60*_indlist[0].interval*3,
                            len(_indlist))
        
        return fig
        
         
    
    def saveafig(self, _fig, op, **kwargs):
        dataname = "_".join([self.groupname])
        #dataname = thedate + "_" + genotype+ "_" + str(samplenum+1) +  "_" +str(interval)
        #_fig.savefig(dataname+"rawfoqar.png",figsize=(8,2),dpi=100)
        #need to add directory
        #_fig.savefig(targetdir +"/"+dataname+"_"+op+"_lethargus.png", **kwargs)
        _fig.savefig(targetdir +"/"+dataname+"_"+op+"_.png", **kwargs)
        
    #_amatix is self.ltfoqmatrix etc.
    def saveltmatrix(self, _amatrix, op, **kwargs):
        ordinalnum = 0
        if "ordinalnum" in kwargs:
            ordinalnum = kwargs["ordinalnum"]
        dataname = "_".join([self.groupname])
        df = pandas.DataFrame(_amatrix.T)
        
        tempdf = self.summarydf[self.summarydf.ltnum == str(ordinalnum)]
        labeldf = tempdf.loc[:,["date", "groupname", "expnum","samplenum"]]
        #colnames =  ["_".join([x.date, x.groupname,
        #                       str(x.expnum),
        #                       str(x.samplenum)]) for x in self.indlist]
        colnames =  ["_".join([x.date, x.groupname,
                               str(x.expnum),
                               str(x.samplenum)]) for i, x in labeldf.iterrows()]
        df.columns = colnames
        interval = self.fullindlist[0].interval
        df.to_csv(targetdir +"/"+dataname+"_"+str(interval)+"_"+op+"_df"+".csv", index=False )

    def showcandidates(self, _samplenum):
        ind = self.fullindlist[_samplenum-1]
        fig = ind.preparefig()
        ax = ind.plotfoq(fig.get_axes()[0])
        #ax = ind.plotcandidates(ax, ind.fq_oescreendmatrix)
        for x in ind.fq_onsetcandidates:
            ax.axvline(x = x, color = "red", linewidth = 0.5, linestyle="--")
        for x in ind.fq_exitcandidates:
            ax.axvline(x = x, color = "blue", linewidth = 0.5, linestyle="--")

        self.printcandidates(_samplenum)
        return fig

    def printcandidates(self, _samplenum):
        print("start "+ ",".join(map(str,self.fullindlist[_samplenum-1].fq_onsetcandidates)))
        print("end "+ ",".join(map(str,self.fullindlist[_samplenum-1].fq_exitcandidates)))
        
    def manualanalysis(self, _samplenum, _start, _end):
        #samplenum start from 1
        tempind = self.fullindlist[_samplenum-1]
        print("tempind.samplenum "+ str(tempind.samplenum))
        measures = tempind.calcmeasures( _start, _end)        
        """
        outputdata = [meanquiescent, meanquiescentout, totalq, numberofbout,
                      qmean, amean, qmedian, amedian, fq_duration, qfreq,
                      qaboutdataframe,lethargusperiodqadf]
        """
        print("------------------")
        print("man_start " + str(_start))
        print("man_end " + str(_end))
        print("man_duration " + str(measures[8]))
        print("totalq " + str(measures[2]))
        print("meanquiescent " + str(measures[0]))
        print("meanquiescentout " + str(measures[1]))
        print("numberofbout " + str(measures[3]))
        print("qmean " + str(measures[4]))
        print("amean " + str(measures[5]))
        print("qmedian " + str(measures[6]))
        print("amedian " + str(measures[7]))
        print("qfreq " + str(measures[9]))

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



samplegroups = []

#for gn in uniquegroupnames:
for (gn, indice) in zip(uniquegroupnames, gsamplenumlist):
    print(gn, indice)
    #this groupname format is date_groupname_expnum
    #this may not good? just use groupname?
    sg = Samplegroup("_".join([thedate, gn, str(expnum)]))
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
filename = "_".join([thedate, groupname, "gridfoq", picformat])
figfilepath = os.path.join(targetdir,filename)
gridfig.savefig(figfilepath,dpi=100)

#240201 for svg format
picformat = ".svg"
filename = "_".join([thedate, groupname, "gridfoq", picformat])
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
filename = "_".join([thedate, groupname, "gridfoqwithflag_meancalc", picformat])
figfilepath = os.path.join(targetdir,filename)
gridfig.savefig(figfilepath,dpi=100)


#240201 for svg format
picformat = ".svg"
#filename = "_".join(["foqmarged_mean","marging",str(marginhr),
#                     ugname, timesuff, picformat])
filename = "_".join([thedate, groupname, "gridfoqwithflag_meancalc", picformat])
figfilepath = os.path.join(targetdir,filename)
gridfig.savefig(figfilepath,dpi=100)


#sys.path.append is not recommended way to set path for module... 
# but for now use easy way
# In the readme.md on git, suggesting that 
# put dotplot.py in the same directory.
# so this line is only for my devlopment env.
sys.path.append(os.path.join(os.path.dirname(__file__), 
                                "..","tkmodules"))
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
filename = "_".join([thedate, groupname, "dotplot", param, "all", picformat])
figfilepath = os.path.join(targetdir,filename)
dotplotfig.savefig(figfilepath,dpi=100)


#240201 for svg format
picformat = ".svg"
#filename = "_".join(["foqmarged_mean","marging",str(marginhr),
#                     ugname, timesuff, picformat])
filename = "_".join([thedate, groupname, "dotplot", param, "all", picformat])
figfilepath = os.path.join(targetdir,filename)
dotplotfig.savefig(figfilepath,dpi=100)





sys.exit()
#############################################################################


