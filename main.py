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
from numpy import savetxt, array,sum, isnan

start_time = time.time() # get sys time to calculate runtime
 
descriptionString = "Takes a list of itensity matrices.\nNormalizes the input (for -q != 0), defines events in the data and converts those to spike trains.\nReturns spike trains and transient/event properties in respective files."
inputParser = argparse.ArgumentParser(description=descriptionString)
inputParser.add_argument('inputFile',nargs="+", help="Input matrices containing mean Intensities",type=str)        
inputParser.add_argument('-o','-outputFilePath',nargs='?',help="Full file path to save the output to", type=str,default=os.path.dirname(os.path.realpath(__file__)))
inputParser.add_argument("-q",'-quantileWidth',"-qWidth", help="Width of the window to be considered for quantile normalization. Set to 0 to skip quantile Norm", nargs='?', type=int,default=0)        
inputParser.add_argument("-b",'-baselineWidth',"-bWidth",help="Width of the window to be considered for baseline normalization", nargs='?',type=int,default=300)                
inputParser.add_argument("-s",'-slopeWidth',help="Distance of frames considered for slope calculation in eventDetect", nargs='?',type=int,default=3)
inputParser.add_argument('-c','-cutoff',help="Distance from slope distribution mean in standard deviance, used for eventDetect", nargs='?',type=float,default=2.0)
inputParser.add_argument('-n','-startNum', help='Number of times slope threshold (set with -c) has to be passed consecutively for a transient to be declared',nargs='?',type=int,default=5)
inputParser.add_argument('-m','-minEventLength', help="Minimum number of data points above transient start fluorescence - if not met, discard transient", nargs="?",type=int,default=30)
inputParser.add_argument('-seperator',"-sep",help="Seperator used for entries in the input matrix - if it is human readable (e.g. csv format).", nargs='?',default="False")        
args = inputParser.parse_args()

apFileList = ["/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/cai-1/9ChenAPs",
              "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/cai-1/11ChenAPs",
              "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/cai-1/akerboomAPs",
              "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/cai-2/cai2_CAMP_full_aps.csv",
              "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/cai-2/cai2_GECO_full_aps.csv"]
outFileList = ["/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/output/9chenDeconvo",
               "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/output/11chenDeconvo",
               "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/output/akerboomDeconvo",
               "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/output/campDeconvo",
               "/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/output/gecoDeconvo",]

counter = -1
for filename in args.inputFile:
    counter += 1
    rawMatrix, numOfVals = cal_neuroIm._importMatrix(filename,args.seperator)
    print("on file:" + str(filename) + " with " + str(numOfVals) + " ROIs")
    time1 = time.time()
    baselineMatrix, baselineCoordinates = cal_neuroIm.pushToBaseline(rawMatrix,args.b)
    print("time for baseline correction: %f" %(time.time() - time1))
    time2 = time.time()
    #transientMatrix = cal_neuroIm.thresholdEventDetect(baselineMatrix,args.q,5,2)
    transientMatrix, whatevsLOL = cal_neuroIm.eventDetect(baselineMatrix, args.q, 2, 1, 1, 10,6)
    print("time for event detection: %f" %(time.time() - time2))
    time3 = time.time()
    
    #TODO: what is a good range for the slope calculations? should I hard-code this too? maybe add parameter sheet option for stuff
    meanKernel = cal_neuroIm.createMeanKernel(transientMatrix)
    print("time for mean kernel: %f" %(time.time() - time3))
    
    #transientMatrix = cal_neuroIm.aucCorrect(transientMatrix, meanKernel)
    #meanKernel = cal_neuroIm.createMeanKernel(transientMatrix)
    apTimings,uselessValue = cal_neuroIm._importMatrix(apFileList[counter]," ")
    '''
    for i in range(numOfVals):#print("decon time: %f vs. total plotting loop time: %f" % (decontime,time.time() - time2))

        import matplotlib.pyplot as plt
        #transientMatrix[i] = cal_neuroIm.aucCorrect(transientMatrix[i], meanKernel)
        #transientMatrix[i] = cal_neuroIm.negativeTransientsCorrect(transientMatrix[i])
        
        f,(axarr0,axarr1,axarr2) = plt.subplots(3, sharex=False)
        axarr0.plot(rawMatrix[:,i],lw=0.2)
        spikeSignal = cal_neuroIm.deconvolve(transientMatrix[i], meanKernel)
        theseTimings = apTimings[:,i][~isnan(apTimings[:,i])]
        titleString =  str(sum(spikeSignal)) + " " + str(sum(theseTimings)) #
        plt.xlabel(titleString)
        #axarr2.hist(slopeDistributions[i][2], 100, normed=1, facecolor='green', alpha=0.75)
        axarr1.plot(spikeSignal,'r',lw=0.2)
        axarr2.set_ylim([0,1])
        axarr2.plot(theseTimings,lw=0.2)
        #axarr2.plot(meanKernel[0:200])
        #axarr2.plot([slopeDistributions[i][3],slopeDistributions[i][3]],[0,max(slopeDistributions[i][2])],'k--',lw=3)
        
        for j in transientMatrix[i]: # yes, iterate over transient objects
            if(transientMatrix[i]): # is there a transient in this time series?
                coords = j.coordinates
                axarr0.plot(range(coords[0],coords[1]),baselineMatrix[coords[0]:coords[1],i],'r',lw=0.2)
        #axarr2.plot(theseTimings,lw=0.1,color = 'black')
        axarr0.plot(baselineCoordinates[i],[0,0],'k.',lw=1)
        
        outputFile = args.o+filename.split("/")[-1] +str(i)+".png"
        
        plt.savefig(outputFile)
        plt.close()
    '''
    time4 = time.time()
    deconvolvedMatrix = array(cal_neuroIm.deconvolve(transientMatrix, meanKernel)).transpose()
    print("time for  deconvolution: %f" %(time.time() - time4))
    
    savetxt(outFileList[counter], deconvolvedMatrix, fmt="%i", delimiter=" ")
    
print("--- %s seconds runtime --- output is in %s" % (time.time() - start_time, args.o))