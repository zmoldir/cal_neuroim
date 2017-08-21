'''
Created on Sep 30, 2016
@author: Maximilian Pfeiffer
@argv: filename, -seperator: file to be opened and unicode string seperating data entries, respectively.
@argv: quantileWidth, baselineWidth: size of buckets to be scanned for the baseline / quantiles, respectively.

    pictures to generate:
    - drift correction & transient smoothing for smaller q-widths 
    - baseline correction / scaling of different ROIs
    - multi-peaked vs. single-peak transients
    - graph of alpha function / alpha kernel results
    - image for the deconvolved time series
    - screenshots for command-line usage of the script // need to go over input args again

'''
import time, argparse, cal_neuroIm, os
from numpy import savetxt, array,sum, isnan,delete, zeros
#'/home/maximilian/unistuff/paris_ens/cal_neuroim/driftData/a5 KO anaesthetized mice data/2014.02.04 a5 cage 18/TSeries-02042014-1019-4562/Analysis/Results.xls'
start_time = time.time() # get sys time to calculate runtime
 
descriptionString = "Takes a list of itensity matrices.\nNormalizes the input (for -q != 0), defines events in the data and converts those to spike trains.\nReturns spike trains and transient/event properties in respective files."
inputParser = argparse.ArgumentParser(description=descriptionString)
inputParser.add_argument('inputFile',nargs="+", help="Input matrices containing mean Intensities",type=str)        
inputParser.add_argument('-o','-outputFilePath',nargs='?',help="Full file path to save the output to", type=str,default=os.path.dirname(os.path.realpath(__file__)))
inputParser.add_argument("-q",'-quantileWidth',"-qWidth", help="Width of the window to be considered for quantile normalization. Set to 0 to skip quantile Norm", nargs='?', type=int,default=0)        
inputParser.add_argument("-b",'-baselineWidth',"-bWidth",help="Width of the window to be considered for baseline normalization", nargs='?',type=int,default=300)                
inputParser.add_argument('-seperator',"-sep",help="Seperator used for entries in the input matrix - if it is human readable (e.g. csv format).", nargs='?',default="False")        
args = inputParser.parse_args()

for filename in args.inputFile:
    
    rawMatrix, numOfVals = cal_neuroIm._importMatrix(filename,args.seperator)
    print("on file:" + str(filename) + " with " + str(numOfVals) + " ROIs")
    baselineMatrix, baselineCoordinates = cal_neuroIm.pushToBaseline(rawMatrix,args.b)
    
    transientMatrix,nopeLol = cal_neuroIm.eventDetect(baselineMatrix,args.q,3)
    #TODO: what is a good range for the slope calculations? should I hard-code this too? maybe add parameter sheet option for stuff
    meanKernel = cal_neuroIm.createMeanKernel(transientMatrix)
    
    #apTimings,uselessValue = cal_neuroIm._importMatrix('/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/cai-2/cai2_CAMP_full_aps.csv'," ")
    #apTimings = delete(apTimings,(0),axis=1)
    for i in range(numOfVals):#print("decon time: %f vs. total plotting loop time: %f" % (decontime,time.time() - time2))

        import matplotlib.pyplot as plt
        #transientMatrix[i] = cal_neuroIm.aucCorrect(transientMatrix[i], meanKernel)
        #transientMatrix[i] = cal_neuroIm.negativeTransientsCorrect(transientMatrix[i])
        f,(axarr0,axarr1,axarr2) = plt.subplots(3, sharex=True)
        axarr0.plot(rawMatrix[:,i],lw=0.5)
        spikeSignal = cal_neuroIm.deconvolve(transientMatrix[i], meanKernel)
        
        #axarr2.hist(slopeDistributions[i][2], 100, normed=1, facecolor='green', alpha=0.75)
        #axarr2.plot(slopeDistributions[i][0],slopeDistributions[i][1],'r--', linewidth=2)
        #axarr2.plot([slopeDistributions[i][3],slopeDistributions[i][3]],[0,max(slopeDistributions[i][2])],'k--',lw=3)
        #apList = apTimings[:,i][~isnan(apTimings[:,i])]
        
        
        #theseTimings = zeros(len(baselineMatrix[:,i]))
        #for j in apList:
        #    theseTimings[int(j)]= j
        #titleString = str(sum(spikeSignal)) + " " + str(sum(apList))
        #axarr1.text(1,2,titleString)
        #axarr1.plot(slopeDistribution[i][2],lw=0.1)
        axarr1.plot(spikeSignal,'r',lw=0.1)
        #axarr1.plot(apList,'black',alpha=0.8,lw=0.1)
        
        for j in transientMatrix[i]: # yes, iterate over transient objects
            if(transientMatrix[i]): # is there a transient in this time series?
                coords = j.coordinates
                axarr0.plot(range(coords[0],coords[1]),baselineMatrix[coords[0]:coords[1],i],'r')
        axarr0.plot(baselineCoordinates[i],[0,0],'k.',lw=1)
        
        outputFile = args.o+filename.split("/")[-1] +str(i)+".png"
        
        plt.savefig(outputFile)
        plt.close()
        
    '''deconvolvedMatrix = array(cal_neuroIm.deconvolve(transientMatrix, meanKernel)).transpose()
    savetxt("/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/deconvolvedCampAps.csv", deconvolvedMatrix, fmt="%i", delimiter=" ")'''
    
print("--- %s seconds runtime --- output is in %s" % (time.time() - start_time, args.o))