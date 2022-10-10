#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 13:54:11 2022
@author: cvaneck

The POSSUM channel flagging script.
This script is broken into two parts:
    -identify bad channels (using Stokes I, Q, and U), and make a list of bad channels
    -modify a set of Stokes cubes (I,Q, U, maybe V?) to NaN out bad channels.


"""
import argparse
import astropy.io.fits as pf
import numpy as np
#from math import ceil,floor
import os


def find_bad_channels(cubename,extreme_limit=1e10,noise_abs=1e5,local_size=10,
                      local_extreme_deviation=5,local_noise_deviation=5):
    """Find the bad channels in a (single-Stokes) FITS cube, by searching
    for channels with unusual statistical properties (far too noisy, extreme values)
    Assumes cube is in HDU 0. Assumes the frequency axis is the last
    non-degenerate axis (generally axis 3 or 4). Works with either 3 or 4D 
    cubes (with freq in either 3rd or 4th axis).
    
    Inputs: 
        cubename (str): path (absolute, or relative to run directory) to FITS cube
        extreme_limit (float, default 1e10): value (in data units) that constitutes 
            an extreme (a value too high to be realistic)
        noise_abs (float, default 1e5): value (in data units) of the largest
            acceptable channel noise.
        local_size (int, default 10): Size of box to use when analyzing local
            spectral properties
        local_extreme_deviation (float, default 5): maximum allowed multiple 
            of the local median extreme (also 1/this value is allowed minimum)
        local_noise_deviation (float, default 5): maximum allowed multiple 
            of the local median noise (also 1/this value is allowed minimum)
    
    Returns: tuple of 3 elements
        array of bools. True = good channel, False = bad channel
        array of floats, containing extreme values per channel
        array of floats, containing noise estimates per channel
    """

    hdulist=pf.open(cubename)
    data=hdulist[0].data
    header = hdulist[0].header
    
    #Get frequency axis:
    for key, value in header.items():
        if type(value) is str:
            if value.lower() == "freq":
                n_axis = int(key[-1])
    nDim=header['NAXIS']

    #Array of extreme values (max/min) per channel
    extrema=np.zeros(header['NAXIS'+str(n_axis)])

    #Noise is calculated per channel using MADFM. It's a bit more expensive
    #than a standard deviation, but it's nicely outlier resistant.
    noise=np.zeros(header['NAXIS'+str(n_axis)])

    #Boolean array of channel quality


    for i in range(extrema.size):
        #This is expanded out manually (3 cases for different file configurations),
        #in order to avoid unncessary operations (which might read the whole cube
        #into memory). Generates a 2D channel map to extrats stats for.
        #Working on individual slices (rather than spectra) should make for
        #the fastest data access.
        if nDim == 4 and n_axis == 4:
            chan=data[i,0,:,:]
        elif nDim == 4 and n_axis == 3:
            chan=data[0,i,:,:]
        elif nDim == 3:
            chan=data[i,:,:]
        else:
            raise Exception('File dimensions do not seem correct? Cannot slice cube.')

        extrema[i]=np.nanmax(np.abs(chan))
        noise[i]=MAD(chan)

    hdulist.close()

    channel_list=test_channels(extrema,noise,extreme_limit,noise_abs,local_size,
                      local_extreme_deviation,local_noise_deviation)


    return channel_list,extrema,noise



def test_channels(extrema,noise,extreme_limit=1e10,noise_abs=1e5,local_size=10,
                      local_extreme_deviation=5,local_noise_deviation=5):
    """Perform tests of channel quality.
    Inputs: 
        cubename (str): path (absolute, or relative to run directory) to FITS cube
        extreme_limit (float, default 1e10): value (in data units) that constitutes 
            an extreme (a value too high to be realistic)
        noise_abs (float, default 1e5): value (in data units) of the largest
            acceptable channel noise.
        local_size (int, default 10): Size of box to use when analyzing local
            spectral properties
        local_extreme_deviation (float, default 5): maximum allowed multiple 
            of the local median extreme (also 1/this value is allowed minimum)
        local_noise_deviation (float, default 5): maximum allowed multiple 
            of the local median noise (also 1/this value is allowed minimum)
    
    Returns: array of bools. True = good channel, False = bad channel
    
    """
    channel_list=np.ones(extrema.size,dtype=bool)
    
    #Determine the local median noise of the local channels (edges will, by 
    #necessity, have fewer local channels.)
    local_median_noise = np.zeros_like(noise)
    local_median_extrema = np.zeros_like(noise)
    for i in range(local_median_noise.size):
        local_median_noise[i] = np.nanmedian(noise[max([0,i-local_size//2]):
                                                   min(local_median_noise.size,i+local_size//2)])
        local_median_extrema[i] = np.nanmedian(extrema[max([0,i-local_size//2]):
                                                   min(local_median_noise.size,i+local_size//2)])
    

    #Remove channels with extremely large values:
    channel_list=np.where(extrema < extreme_limit,channel_list,False)
    #Remove channels with extremely large noise:
    channel_list=np.where(noise < noise_abs,channel_list,False)
    #Remove channels with anomalously high noise compared to neighbouring channels:
    channel_list=np.where(noise < local_median_noise*local_noise_deviation,
                          channel_list,False)
    #Remove channels with anomalously low noise compared to neighbours:
    channel_list=np.where(noise > local_median_noise/local_noise_deviation,
                          channel_list,False)

    #Remove channels with anomalously high extremes compared to neighbouring channels:
    channel_list=np.where(extrema < local_median_extrema*local_extreme_deviation,
                          channel_list,False)
    #Remove channels with anomalously low extremes compared to neighbours:
    channel_list=np.where(extrema > local_median_extrema/local_extreme_deviation,
                          channel_list,False)

    return channel_list




def check_all_cubes(cubename,extreme_limit=1e10,noise_abs=1e5,local_size=10,
                      local_extreme_deviation=5,local_noise_deviation=5):
    """Finds band channels for all Stokes cubes present (including V if file 
    exists), and combines (if a channel is bad in one Stokes, it's set bad for
    all). Uses find_bad_channels() on each Stokes cube; see that function for
                           description of method.
    Inputs: 
        cubename (str): path+filename to a Stokes cube. It is a assumed to have
            a ASKAP-style filename with a '.i.' to denote Stokes parameter.
            Any Stokes parameter can be provided; the other files will be 
            inferred to have the same names other than the Stokes parameter.
    
    """

    generic_filename=cubename.replace('.i.','._.').replace('.q.','._.').replace('.u.','._.').replace('.v.','._.')

    print('Starting Stokes I characterization.')
    I_channels,I_extrema,I_noise=find_bad_channels(generic_filename.replace('._.','.i.'),
                                 extreme_limit,noise_abs,local_size,
                                 local_extreme_deviation,local_noise_deviation)
    print('Starting Stokes Q characterization.')
    Q_channels,Q_extrema,Q_noise=find_bad_channels(generic_filename.replace('._.','.q.'),
                                 extreme_limit,noise_abs,local_size,
                                 local_extreme_deviation,local_noise_deviation)
    print('Starting Stokes U characterization.')
    U_channels,U_extrema,U_noise=find_bad_channels(generic_filename.replace('._.','.u.'),
                                 extreme_limit,noise_abs,local_size,
                                 local_extreme_deviation,local_noise_deviation)
    
    #If V is present, include it. Otherwise, assume V has no problems.
    if os.path.exists(generic_filename.replace('._.','.v.')):
        print('Starting Stokes V characterization.')
        V_channels,V_extrema,V_noise=find_bad_channels(generic_filename.replace('._.','.v.'),
                                 extreme_limit,noise_abs,local_size,
                                 local_extreme_deviation,local_noise_deviation)
    else:
        V_channels=np.ones_like(I_channels)
        V_extrema=np.zeros_like(I_channels)
        V_noise=np.zeros_like(I_channels)

    #Since True = good, we require all Stokes have True for a channel to be good
    combined_channels=np.all([I_channels,Q_channels,U_channels,V_channels],axis=0)
    
    
    
    
    return combined_channels, np.array([I_extrema,I_noise,Q_extrema,
                                        Q_noise,U_extrema,U_noise,V_extrema,V_noise])





def flag_cube(cubename,channel_flags):
    """Applies flags to a FITS cube. This is applied in-place; no new files
    are made. Bad channels are replaced entirely with NaNs.
    
    Inputs:
        cubename (str): filename/path to the FITS cube to be flagged.
            Must be 3 or 4 dimensional.
        channel_flags (array-like): Boolean array of channel quality.
            True = good, False = bad/to-be-flagged
            
    Returns: none
    """
    hdulist=pf.open(cubename,mode='update')
    data=hdulist[0].data
    header = hdulist[0].header
    
    #Get frequency axis:
    for key, value in header.items():
        if type(value) is str:
            if value.lower() == "freq":
                n_axis = int(key[-1])
    nDim=header['NAXIS']

    for i in range(channel_flags.size):
        if channel_flags[i] == True:
            continue #skip good channels
        
        #This is expanded out manually (3 cases for different file configurations),
        #in order to ensure correct slicing.
        if nDim == 4 and n_axis == 4:
            data[i,0,:,:]=np.nan
        elif nDim == 4 and n_axis == 3:
            data[0,i,:,:]=np.nan
        elif nDim == 3:
            data[i,:,:]=np.nan
        else:
            raise Exception('File dimensions do not seem correct? Cannot flag cube.')

    hdulist.close()
    
    
    
def flag_all_Stokes(cubename,channel_flags):
    """Applies flag_cube() to a set of Stokes cube (I, Q, U, and V if present).
    
    Inputs:
        cubename (str): filename/path to the FITS cube to be flagged.
            Must be 3 or 4 dimensional. It is a assumed to have
            an ASKAP-style filename with a '.i.' to denote Stokes parameter.
            Any Stokes parameter can be provided; the other files will be 
            inferred to have the same names other than the Stokes parameter.
        channel_flags (array-like): Boolean array of channel quality.
            True = good, False = bad/to-be-flagged
        
    """
    generic_filename=cubename.replace('.i.','._.').replace('.q.','._.').replace('.u.','._.').replace('.v.','._.')

    print('Starting Stokes I flagging.')
    flag_cube(generic_filename.replace('._.','.i.'),channel_flags)
    print('Starting Stokes Q flagging.')
    flag_cube(generic_filename.replace('._.','.q.'),channel_flags)
    print('Starting Stokes U flagging.')
    flag_cube(generic_filename.replace('._.','.u.'),channel_flags)
    
    #If V is present, include it. Otherwise, assume V has no problems.
    if os.path.exists(generic_filename.replace('._.','.v.')):
        print('Starting Stokes V flagging.')
        flag_cube(generic_filename.replace('._.','.v.'),channel_flags)






def command_line():
    """When invoked from the command line, reads the arguments and processes
    the supplied cubes. Default behaviour is to only compute the bad channels;
    to flag them you must also set the appropriate flag (this is a safety feature)
    to prevent accidentally overwriting the cubes.
    Note that flagging the cubes is done 'in-place' with the supplied cubes;
    NO NEW CUBES ARE GENERATED, so make sure you want the cubes modified before
    selecting that option.
    
    """

    # Help string to be shown using the -h option
    descStr = """
    Characterizes channels and determines bad channels, and/or flags bad 
    channels. Designed for ASKAP cubes, assumes each Stokes has ".i." type 
    filename, and assumes Stokes IQU are all present with the same file names.
    If Stokes V is present, it will detect it and use it as well.
    Will produce a text file with the channel stats.
    Default behaviour is to only characterize and classify bad channels;
    to flag them you must also set the -f flag (this is a safety feature
    to prevent accidentally overwriting the cubes).
    Note that flagging the cubes is done 'in-place' with the supplied cubes;
    NO NEW CUBES ARE GENERATED, so make sure you want the cubes modified before
    selecting that option.
    If the stats have already been computed, only the flagging done by setting 
    the -f and -s flags.
    """


    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,#epilog=epilog_text,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("cubefile", metavar="cube.fits", 
                        help="FITS cube (any Stokes)")
    parser.add_argument("statsfile", metavar="stats.dat", 
                        help="Output channel statistics file.")
    parser.add_argument("-f", dest="applyflags", action="store_true",
                        help="Flag bad channels (set to NaNs) (modifies supplied cubes!) [False].")
    parser.add_argument("-s", dest="skipstats", action="store_true",
                        help="Skip stats computation (used supplied stats file) [False].")

    parser.add_argument("-e", dest="extreme_limit", type=float, default=1e10,
                        help="Intensity threshold for extreme values. [1e10]")
    parser.add_argument("-n", dest="noise_limit", type=float, default=1e5,
                        help="Noise threshold for extreme values. [1e5]")
    parser.add_argument("-l", dest="local", nargs=3, action='store',
                        type=float, default=[10, 5, 5],
                        help="Number of channels for local estimation, threshold multiple of local noise, threshold multiple for local extrema. [10 5 5]")
    


    args = parser.parse_args()


    if not os.path.exists(args.cubefile):
        raise Exception('Cube file could not be found!')

    extreme_limit=args.extreme_limit
    noise_limit = args.noise_limit
    local_size = int(args.local[0])
    local_noise_threshold = args.local[1]
    local_extreme_threshold = args.local[2]

    if not args.skipstats: #Calculate stats
        channel_flags,channel_stats=check_all_cubes(args.cubefile,extreme_limit,
                                        noise_limit,local_size,
                                        local_noise_threshold,local_extreme_threshold)
        
        outdata=np.vstack((channel_flags,channel_stats))
        np.savetxt(args.statsfile,outdata)
    else: #Or read them from file:
        if not os.path.exists(args.statsfile):
            raise Exception('Stats file could not be found; required as input with -s flag set!')
        indata=np.loadtxt(args.statsfile)
        channel_flags=indata[0]
    
    if args.applyflags:
        flag_all_Stokes()



def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:
    median(abs(a - median(a))) / c
    c = 0.6745 is the constant to convert from MAD to std
    
    This function blatently ripped from RM-Tools; I don't want to make it a
    dependency, so I'm just copying the function.
    """
    
    a = np.ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = np.ma.median(a)
        m = np.ma.median(np.ma.fabs(a - d) / c)
    else:
        d = np.ma.median(a, axis=axis)
        if (axis is not None) and (axis > 0):
            aswp = np.ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = np.ma.median(np.ma.fabs(aswp - d) / c, axis=axis)

    return m



if __name__ == '__main__':
    command_line()
