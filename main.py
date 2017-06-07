'''
Created on Sep 30, 2016
@author: Maximilian Pfeiffer
@argv: filename, -seperator: file to be opened and unicode string seperating data entries, respectively.
@argv: quantileWidth, baselineWidth: size of buckets to be scanned for the baseline / quantiles, respectively.

    general process outline:
    
    7. Documentation! 
        - README for the done script
        - explanation for methods
        - check transient AUC

'''
import time, argparse, cal_neuroIm
from numpy import mean, isnan, squeeze, zeros
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
    
    optiMethodNorm = zeros([5000,numOfVals])
    optiMethodSpikes = zeros([5000,numOfVals])
    optiTime = time.time()
    for i in range(numOfVals):
        optiMethodObject = cal_neuroIm.constrained_foopsi(rawMatrix[:,i])
        optiMethodNorm[i] = optiMethodObject[0]
        optiMethodSpikes[i] = optiMethodObject[-1]
    print ("time for optimization approach: %f" % (time.time() - optiTime))
    
    baselineMatrix, baselineCoordinates = cal_neuroIm.pushToBaseline(rawMatrix,args.b)
    transientMatrix, quantileMatrix, slopeDistributions = cal_neuroIm.eventDetect(baselineMatrix, args.q,len(rawMatrix)/150)
    
    meanKernel = cal_neuroIm.createMeanKernel(transientMatrix)
    #print ("time for mean kernel calculation: %f") % (time.time() - time1)
    apTimings,uselessValue = cal_neuroIm.importMatrix('/home/maximilian/unistuff/paris_ens/cal_neuroim/simdata/apTimings.csv',args.seperator)
    #TODO: array containing only the transient data -> check negative transients first and filter!
    #ALSO: plot slope distribution, pack decon techniques into one picture -> where should cutoff be?

    for i in range(numOfVals):
        import matplotlib.pyplot as plt
        transientMatrix[i] = cal_neuroIm.aucCorrect(transientMatrix[i], meanKernel)
        f,(axarr0,axarr1,axarr2) = plt.subplots(3, sharex=False)
        eventCoordinates = [j.getCoordinates() for j in transientMatrix[i]] # ugly hack that defeats the purpose. REFACTOR
        time1 = time.time()
        spikeSignal = cal_neuroIm.deconvolve(transientMatrix[i], meanKernel)
        
        spikeTrain = cal_neuroIm.generateSpiketrainFromSignal(spikeSignal)
        '''
        axarr2.hist(slopeDistributions[i][2], 100, normed=1, facecolor='green', alpha=0.75)
        axarr2.plot(slopeDistributions[i][0],slopeDistributions[i][1],'r--', linewidth=2)
        axarr2.plot([slopeDistributions[i][3],slopeDistributions[i][3]],[0,50],'k--',lw=3)
        '''
        axarr2.plot(optiMethodNorm[i])
        axarr1.plot(spikeSignal,'b')
        plt.grid(True)
        axarr1.plot(spikeTrain,'r-')
        
        axarr0.plot(optiMethodSpikes[i],"r")
    
        for x in squeeze(apTimings[i]):
            if ~isnan(x):
                axarr1.plot([int(x),int(x)],[0,1],'g--',lw=2)
        
        thisMean2 = mean(baselineMatrix[:,i])
        '''for j in transientMatrix[i]: # yes, iterate over transient objects
            if(transientMatrix[i]): # is there a transient in this time series?
                coords = j.getCoordinates()
                axarr0.plot(coords,[0,1],'g^', lw=2)
        axarr0.plot(baselineCoordinates[i],[thisMean2,thisMean2],'k.',lw=1)
        '''
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
#print("decon time: %f vs. total plotting loop time: %f" % (decontime,time.time() - time2))
print("--- %s seconds ---" % (time.time() - start_time))

