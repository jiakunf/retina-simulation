from numpy import *
from scipy import *


def RGCresponse(stimulus, RFradius, RFloc, RGCtype, sample_rate):
    # INPUTS: stimulus x*y*time/sampleRate
    # INPUTS: RFradius in [0,1] space
    # INPUTS: RFloc, RF center in [0,1]x[0,1] space
    # INPUTS: RGCtype, 0 if ON, 1 if OFF
    
    # parameters
    T = stimulus.shape[2] * sample_rate   # duration of stimulus in s
    filter_duration = int(32 / (1000 / (1 / sample_rate))) # must be integer
    nlPoint = -0.9
    nlSlope = 2.
    max_firing_rate = 100.
    #sampleRate = T/stimulus.shape[2]
    on_weights  = .7  # center-surround
    off_weights = -.2    # center-surround
    frequency  = 1.   # of biphasic filter
    phase      = -0.5    # of biphasic filter
    variance   = 1.    # variance of biphasic filter

    # generate stimulus

    #pdb.set_trace()
    # assuming RFradius and RFloc are in [0,1] x [0,1] space
    center_radius = 0.5*RFradius*mean([stimulus.shape[0],stimulus.shape[1]])
    yLength = stimulus.shape[0]
    xLength = stimulus.shape[1]
    pixel_center = [int(RFloc[0]*xLength), int(RFloc[1]*yLength)]
    pixel_width  = [max(int(pixel_center[0]-center_radius),int(0)), min(int(pixel_center[0]+center_radius),int(stimulus.shape[0]-1))]
    pixel_height = [max(int(pixel_center[1]-center_radius),int(0)), min(int(pixel_center[1]+center_radius),int(stimulus.shape[1]-1))]

    if RGCtype == 0:    # ON type
        spatial_output = zeros(stimulus.shape[2])
        for i in arange(pixel_width[0], pixel_width[1]+1):
            for j in arange(pixel_height[0], pixel_height[1]+1):
                if ((i-pixel_center[0])**2+(j-pixel_center[1])**2) < (center_radius/2)**2:
                    spatial_output += stimulus[i,j,:]*on_weights
                elif ((i-pixel_center[0])**2+(j-pixel_center[1])**2) < (center_radius)**2:
                    spatial_output += stimulus[i,j,:]*off_weights
    elif RGCtype == 1:  # OFF type
        spatial_output = zeros(stimulus.shape[2])
        for i in arange(pixel_width[0], pixel_width[1]+1):
            for j in arange(pixel_height[0], pixel_height[1]+1):
                if ((i-pixel_center[0])**2+(j-pixel_center[1])**2) < (center_radius/2)**2:
                    spatial_output += stimulus[i,j,:]*off_weights
                elif ((i-pixel_center[0])**2+(j-pixel_center[1])**2) < (center_radius)**2:
                    spatial_output += stimulus[i,j,:]*on_weights
   #       spatialFilter = zeros((stimulus.shape[0],stimulus.shape[1]))
 #       for i in     arange(pixel_width[0],  pixel_width[1]-1):
 #           for j in arange(pixel_height[0], pixel_height[1]-1):
 #               if sqrt((i-pixel_center[0])**2+(j-pixel_center[1])**2) < center_radius:
 #                   spatialFilter[i,j]     = on_weights
 #                   if sqrt((i-pixel_center[0])**2+(j-pixel_center[1])**2) < center_radius/2:
 #                       spatialFilter[i,j] = off_weights

#    spatial_output = zeros(stimulus.shape[2])
#    for t in range(stimulus.shape[2]):
#        projection = array(stimulus[:,:,t])*array(spatialFilter)
#        spatial_output[t] = sum(projection)

    #pdb.set_trace()

    kernel = linearKernel(frequency,phase,variance,filter_duration)
    linearOutput = convolve(spatial_output,kernel,'valid')   # note this is max - min + 1
    #nonlinear_output = threshold(linearOutput, nlPoint, nlSlope)

    # normalize nonlinear output to be between 0 and 1
    #nonlinear_output = nonlinear_output - min(nonlinear_output)
    #nonlinear_output = nonlinear_output/max(nonlinear_output)
    nonlinear_output = exp(linearOutput)/(1+exp(linearOutput))

    spikes = []
    spike_times = []
    #timesGen = drange(0,T,sampleRate)
    times = arange(0,len(nonlinear_output))  # was stimulus.shape[2]
    times = times * sample_rate
    #times = ["%g" % x for x in timesGen]
    print("{} time samples".format(len(times)))
    for t in range(len(nonlinear_output)):
        spikes.append(random.poisson(max_firing_rate * sample_rate * nonlinear_output[t]))
        if spikes[-1]:
            now = 1000.*array(times[t], float64) + random.rand(1)
            spike_times.append(now)

    return spike_times

def linearKernel(freq,phase,variance,width):
    #needs to be fixed to be a proper high pass filter
    x = linspace(0,2*pi,width)
    sinusoid = sin(freq*x + phase)
    y = gaussian(0,variance,10e3)
    gauss = histogram(y,width)[0]*1.0/10e3
    kernel = sinusoid*gauss
    kernel = kernel/max(kernel)
    return kernel


def threshold(linearResponse, nlPoint, nlSlope):
    result = zeros(size(linearResponse))

    for i in range(size(linearResponse)):
        if linearResponse[i] >= nlPoint:
            result[i] = nlSlope*(linearResponse[i] - nlPoint)
    return result


def gaussian(mean,var,numSamples):
    samples = sqrt(var)*random.randn(numSamples) + mean
    return samples

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

