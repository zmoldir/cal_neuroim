'''
Created on Sep 30, 2016
@author: Maximilian Pfeiffer
For a detailed explanation of how to run the script, see the PDF that came with this
@argv: filename, -seperator: file to be opened and unicode string seperating data entries, respectively.
@argv: quantileWidth, baselineWidth: size of buckets to be scanned for the baseline / quantiles, respectively.

'''

import time, argparse, cal_neuroIm, os
from numpy import savetxt, array

start_time = time.time() # get sys time to calculate runtime
 
# all of the below: parse command line input into python arguments
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

# this loop iterates over the number of input files
for filename in args.inputFile:
    
    rawMatrix, numOfVals = cal_neuroIm._importMatrix(filename,args.seperator) # cast the file to a python numpy array
    
    print("on file:" + str(filename) + " with " + str(numOfVals) + " ROIs") 
    
    baselineMatrix, baselineCoordinates = cal_neuroIm.pushToBaseline(rawMatrix,args.b) # baseline correction and coordinates for visualization
    
    transientMatrix, visualizationVar = cal_neuroIm.eventDetect(baselineMatrix, args.q, 2, 1, 1, 10,6) # annotation of transients
    
    meanKernel = cal_neuroIm.createMeanKernel(transientMatrix) # create kernel for the deconvolution from the list of transients
    
    deconvolvedMatrix = array(cal_neuroIm.deconvolve(transientMatrix, meanKernel)).transpose() # deconvolve the input, yields a spike train
    
    savetxt(args.o, deconvolvedMatrix, fmt="%i", delimiter=" ") # write the spike train to a file
    
print("--- %s seconds runtime --- output is in %s" % (time.time() - start_time, args.o))