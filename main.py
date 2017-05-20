'''
Created on Sep 30, 2016
@author: Maximilian Pfeiffer
@argv: filename, -seperator: file to be opened and unicode string seperating data entries, respectively.
@argv: quantileWidth, baselineWidth: size of buckets to be scanned for the baseline / quantiles, respectively.

    general process outline:
    
    1. finish deconvolution
    2. generation of plots / transient property files into seperate functions, as well as INPUT
    3. All normalization / correction methods shouldn't take more input than necessary and have defaults
    4. 2 & 3 should collapse main.py to an overseeable format! -> if not: rinse & repeat
    5. COMPARE PERFORMANCES
        -how fast is it compared to alternatives, OPTIMIZE
        -how close is the result? OPTIMIZE
    6. Put methods in a library-esque format
    7. Documentation! 
        - README for the done script
        - explanation for methods
TIMELINE:
    1: by (including) 24th of April
    2: 26.4
    3: 28.4
    4: 1.5
    5: 7.5
    6. 9.5
    7. 11.5
'''
import time, argparse, cal_neuroIm
from numpy import mean, isnan, squeeze
#'/home/maximilian/unistuff/paris_ens/cal_neuroim/driftData/a5 KO anaesthetized mice data/2014.02.04 a5 cage 18/TSeries-02042014-1019-4562/Analysis/Results.xls'
start_time = time.time() # get sys time to calculate runtime
 
descriptionString = "Takes a list of itensity matrices.\nNormalizes the input (for -q != 0), defines events in the data and converts those to spike trains.\nReturns spike trains and transient/event properties in respective files."
inputParser = argparse.ArgumentParser(description=descriptionString)
inputParser.add_argument('inputFile',nargs="+", help="Input matrices containing mean Intensities",type=str)        
inputParser.add_argument('outputFilePath',help="Full file path to save the output to", type=str)
inputParser.add_argument("-q",'-quantileWidth',"-qWidth", help="Width of the window to be considered for quantile normalization. Set to 0 to skip quantile Norm", nargs='?', type=int,default=150)        
inputParser.add_argument("-b",'-baselineWidth',"-bWidth",help="Width of the window to be considered for baseline normalization", nargs='?',type=int,default=150)                
inputParser.add_argument('-seperator',"-sep",help="Seperator used for entries in the input matrix - if it is human readable (e.g. csv format).", nargs='?',default="False")        
args = inputParser.parse_args()

for filename in args.inputFile:
    
    rawMatrix, numOfVals = cal_neuroIm.importMatrix(filename,args.seperator)
    print("on file:" + str(filename) + " with " + str(numOfVals) + " ROIs")
    time1 = time.time()
    baselineMatrix, baselineCoordinates = cal_neuroIm.pushToBaseline(rawMatrix,args.b)
    print ("time for baseline norm: %f") % (time.time() - time1)
    time1 = time.time()
    
    transientMatrix, quantileMatrix = cal_neuroIm.eventDetect(baselineMatrix, args.q,len(rawMatrix)/150)
    print ("time for event detection: %f") % (time.time() - time1)
    
    time1 = time.time()
    meanKernel = cal_neuroIm.createMeanKernel(transientMatrix)
    print ("time for mean kernel calculation: %f") % (time.time() - time1)
    #apTimings,uselessValue = cal_neuroIm.importMatrix('/home/maximilian/unistuff/paris_ens/cal_neuroim/simdata/apTimings.csv',args.seperator,args.csv)
    #TODO: array containing only the transient data -> check negative transients first and filter!
    decontime = 0
    time2 = time.time()
    for i in range(numOfVals):
        import matplotlib.pyplot as plt
        f,(axarr0,axarr1,axarr2) = plt.subplots(3, sharex=True)
        eventCoordinates = [j.getCoordinates() for j in transientMatrix[i]] # ugly hack that defeats the purpose. REFACTOR
        time1 = time.time()
        spikeSignal = cal_neuroIm.deconvolve(transientMatrix[i], meanKernel)
        spikeTrain = cal_neuroIm.generateSpiketrainFromSignal(spikeSignal)
        
        axarr2.plot(spikeTrain)
        
        axarr2.plot(cal_neuroIm.aucDeconvolve(transientMatrix[i], meanKernel),'^g')
        decontime += (time.time() - time1)
        plt.grid()
        axarr1.plot(spikeSignal)
        '''
        for x in squeeze(apTimings[i]):
            if ~isnan(x):
                axarr1.plot([int(x),int(x)],[0,5],'g--',lw=2)
        '''
        axarr0.plot(rawMatrix[:,i],"r")
        thisMean2 = mean(baselineMatrix[:,i])
        if(transientMatrix[i]): # is there a transient in this time series?
            for j in transientMatrix[i]: # yes, iterate over transient objects
                coords = j.getCoordinates()
                axarr0.plot(coords,[0,max(j.getData())],'k^', lw=2)
        axarr0.plot(baselineCoordinates[i],[thisMean2,thisMean2],'k.',lw=1)
        plt.savefig(str(filename) + str(i)+".png")
        plt.close()
    '''try:    
        with open(str(args.outputFilePath) + str(filename).split('/')[-1] + ".csv", 'w') as outFile:
            outFile.write("Amplitude\tRiseTime\tDecayTime\tmeanIntensity\ttotalLength\tnumOfPeaks\tstartTime\tendTime\n")    
            for num,t in enumerate(transientList):
                outFile.write(str(t.getAmplitude()) + "\t" + str(t.getRiseTime())+ "\t" + str(t.getDecayTime()) + "\t" + str(t.getMeanIntensity()) + "\t" + str(t.getTotalTime()) +"\t" + str(t.getNumOfPeaks())+"\t" + str(t.getStartTime()) + "\t"+ str(t.getEndTime()) + "\n")
                if(t.getIsLast()):
                    outFile.write("__________________________ROI" + str(num) + "__________________________\n")
    except IOError:
        print(str(args.outputFilePath) +"output" +  str(filename).split('\\')[-1] + " is not a valid file path for output!")
        exit()
        '''
print("decon time: %f vs. total plotting loop time: %f" % (decontime,time.time() - time2))
print("--- %s seconds ---" % (time.time() - start_time))

