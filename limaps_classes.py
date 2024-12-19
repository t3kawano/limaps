"""
Creation time: 2024/03/04 10:33:18

@author: tk
"""

import numpy as np
import pandas
import csv
import matplotlib.pyplot as plt


#region Lethargus class. contains area and foq data of individual,
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

#region Lethargusdetector
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
    ## Hayashi lab critera (old; around 2018?)
    ## quiescence onset time was defined as a time point 
    ## after which the fractional quiescence remains > 0.05 for at least 60 minutes. 
    ## The quiescence exit time was defined as the time point 
    ## after which the fractional quiescence had reached < 0.1. 
    # 
    # currentry, 2024, use >0.2 foq >1.5 hr as default
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
        #20180123 changed to longer than 2hrs
        #continuouslength = int(2*60*60/self.interval)
        #20180223 may be 1.5?
        continuouslength = int(self.minduration*60*60/self.interval)
        #need 1h (3600 sec) pre/post lethargus duration
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
                    print("ll start at "+str(ll[i].start/(60*60/self.ind.interval))+" has problem")
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

#region Individual.calcrawarea =calcrawarea

class Individual:
    
    def __init__(self, date, groupname, expnum, interval, samplenum, rawdata, targetdir):
        #data obtained from .xls file name
        self.date = date
        self.groupname = groupname
        self.expnum = expnum
        self.interval = interval

        
        #one .xls file contains several individuals
        self.samplenum = samplenum
        
        #data        
        self.rawdata = rawdata
        
        #240304 according to class file separation
        self.targetdir = targetdir

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
        _fig.savefig(self.targetdir +"/"+dataname+figexp +".png",dpi=100)
        
    #180226 not yet handling multiple lethargus
    def savefoqincsvfile(self):        
        #save the lethargus period and foq with the date_groupname_expnum_# name
        dataname = "_".join([self.date, self.groupname, str(self.expnum), 
                             str(self.samplenum),"_"])
        with open(self.targetdir +"/"+dataname+"foq"+".csv","w", newline='') as file:
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
            

#region Samplegroup
class Samplegroup:   
    def __init__(self, _groupname, targetdir):
        #date_groupname_expnum format. may be need change to only groupname later?
        self.groupname = _groupname
        #240304 class file separation realted
        self.targetdir = targetdir

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
                             dataframe.loc[:,dataframe.columns[i]],
                             self.targetdir)
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
        #emptymatrix = np.zeros(shape = ( sg.numofadequatesample, 
        emptymatrix = np.zeros(shape = (self.numofadequatesample, 
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
        self.summarydf.to_csv(self.targetdir +"/"+dataname+"_"+str(int(interval))+"_summary"+".csv", index=False )
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
        _fig.savefig(self.targetdir +"/"+dataname+"_"+op+"_.png", **kwargs)
        
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
        df.to_csv(self.targetdir +"/"+dataname+"_"+str(interval)+"_"+op+"_df"+".csv", index=False )

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
