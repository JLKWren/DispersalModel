#!/usr/bin/python

#matplotlib needs numpy, and possibly the dev files for freetype
import matplotlib.pyplot as plt
import sys
import os
import argparse
from os.path import exists
from math import sqrt, log, atan, cos, sin
from random import seed, random


# just to make sure we randomize the random function
seed()


# uses new code to replace bathy check with simple query of current field
# uses newer bathy file from Johanna, does 5km radius, and skips days not of interest.
# mod from goby1b, this does Ed's simulation, note the out of bounds settings were obsolete and wrong, no harm I think but cluttered screen bounds
# mod from 1A, this is for MacPro desktop, check mem limits...
# this does part of the second year...
# modified from post-crash version of hycom4x.tru which came from MyBook...
# This does goby simulation with streams source/sink
# ============
# tried 15km radius instead of "default" 25km
# mod from hycom3.tru to look at additional years and use less diffusivity.
# note the trajectory file has an additional record to indicate settled or not
# Note that dimensioning and arraystart/arrayend values need updating if doing simulation for other time frames



# DLS:
# In many of the for loops, when we define a range, we are adding 1 to the end of the range.
# This is due to python not including the final value in the range [n, x) .  To match what the original  
# Basic code does, we need this +1.  Other work arounds could be devised, but this was the simplest.


# GLOBAL STATIC VALUES 
la1 = -0.000003437984
la2 = 0.0004652059
la3 = -0.001507974
la4 = 110.572
lo1 = 0.00006610779
lo2 = -0.02041057
lo3 = 0.064812
lo4 = 111.1078
TPI = 8 * atan(1)
TIME_CONVERSION = 60 * 60 * 24
diffco = 0.00175
sqrt_diffco = sqrt(diffco)
EVEL_NVEL_ERR = -99


def release(xlon, xlat, startday, uval, vval, duration):
    xflag1 = False
    xflag2 = False
    for day in xrange(1, duration + 1):
        if xflag1 == False:
            currentday = startday + day
            xlon, xlat = disperse(xlon, xlat, currentday, uval, vval)
            xflag1 = checkbounds(xlon, xlat)
            if day < duration and xflag1 == False:
                xflag2 = checkdepth(xlon, xlat, currentday, uval, vval)
        if xflag2 == True:
            xlon, xlat = tweak(xlon, xlat, uval, vval, currentday)
            xflag2 = False # DLS: needs to be reset, since tweak in the TRU would normally set it to False, since it is global    
    return xlon, xlat


def calculateLonLatCM(avglat):
    avglat2 = avglat**2
    avglat3 = avglat**3
    zkmlat = (la1 * avglat3) + (la2 * avglat2) + (la3 * avglat) + la4
    zkmlon = (lo1 * avglat3) + (lo2 * avglat2) + (lo3 * avglat) + lo4
    zlatkm = 1.0 / zkmlat
    zlonkm = 1.0 / zkmlon
    #zlatm = zlatkm / 1000.0
    #zlonm = zlonkm / 1000.0
    #zlatcm = zlatm / 100.0
    #zloncm = zlonm / 100.0
    zlatcm = zlatkm / 100000.0
    zloncm = zlonkm / 100000.0
    return zloncm, zlatcm


def disperse(xlon, xlat, currentday, uval, vval):
    hycomu = 0
    hycomv = 0
    zloncm, zlatcm = calculateLonLatCM(abs(xlat))

    nnlon = int(round( ((xlon - 194.0) / 0.04), 0)) + 1
    nnlat = int(round( ((xlat - 16.0) / 0.04), 0))  + 1
    if (1 <= nnlon <= 401) and (1 <= nnlat <= 251):
        hycomu = uval[currentday][nnlon][nnlat]
        hycomv = vval[currentday][nnlon][nnlat]
        if hycomu != EVEL_NVEL_ERR:
            hycomu *= 100.0
        else:
            hycomu = 0
        if hycomv != EVEL_NVEL_ERR:
            hycomv *= 100.0
        else:
            hycomv = 0
    X1 = nonZeroRandom()
    X2 = nonZeroRandom()
    TPIX2 = TPI * X2
    sqrtlogX1 = sqrt(-2.0 * log(X1))
    XX1 = cos(TPIX2) * sqrtlogX1
    XX2 = sin(TPIX2) * sqrtlogX1

    degu = zloncm * TIME_CONVERSION * hycomu
    degv = zlatcm * TIME_CONVERSION * hycomv
    degu += XX1 * sqrt_diffco
    degv += XX2 * sqrt_diffco
    xlon += degu
    xlat += degv
    return xlon, xlat


def checkbounds(xlon, xlat):
    if xlat < 16.0 or xlat > 26.0 or xlon < 194.0 or xlon > 210.0:
        return  True
    return False


def checkdepth(xlon, xlat, currentday, uval, vval):
    nnlon = int(round( ((xlon - 194.0) / 0.04), 0)) + 1
    nnlat = int(round( ((xlat - 16.0) / 0.04), 0)) + 1
    if 1 <= nnlon <= 401 and 1 <= nnlat <= 251:
        hycomu = uval[currentday][nnlon][nnlat]
        hycomv = vval[currentday][nnlon][nnlat]
        if hycomu == EVEL_NVEL_ERR or hycomv == EVEL_NVEL_ERR:
            return  True
    return False


def tweak(xlon, xlat, uval, vval, currentday):
    newxlon = 0.0
    newxlat = 0.0
    zloncm, zlatcm = calculateLonLatCM(abs(xlat))
    deguconst = zloncm  * TIME_CONVERSION
    degvconst = zlatcm  * TIME_CONVERSION

    xflag2 = True   
    # DLS
    # This was originally recursive, but switched it to iterative.
    # Only difference between recursive calls seemed to be 
    # X1, X2, degu and degv.
    while xflag2:
        degu = deguconst * (50.0 - (random() * 100.0))
        degv = degvconst * (50.0 - (random() * 100.0))
        X1 = nonZeroRandom()
        X2 = nonZeroRandom()
        TPIX2 = TPI * X2
        sqrtlogX1 = sqrt(-2.0 * log(X1))
        XX1 = cos(TPIX2) * sqrtlogX1
        XX2 = sin(TPIX2) * sqrtlogX1
        degu += XX1 * sqrt_diffco
        degv += XX2 * sqrt_diffco
        newxlon = xlon + degu
        newxlat = xlat + degv
        xflag2 = checkdepth(newxlon, newxlat, currentday, uval, vval)
    return newxlon, newxlat


def settle(endlon, endlat, uval, vval, habilon, habilat):
    # DLS
    # dmin was originally 5000, but check at end requires less than 5.  
    # When  xflag3 is false, dmindex is not needed in parent function
    dmin = 6  
    dmindex = None
    for ijk in xrange(1, 302 + 1):
        avglat = abs((endlat + habilat[ijk]) / 2.0)
        avglat2 = avglat**2
        avglat3 = avglat**3
        
        zkmlat = (la1 * avglat3) + (la2 * avglat2) + (la3 * avglat) + la4
        zkmlon = (lo1 * avglat3) + (lo2 * avglat2) + (lo3 * avglat) + lo4

        latdiff = (habilat[ijk] - endlat) * zkmlat
        londiff = (habilon[ijk] - endlon) * zkmlon
        tmp = sqrt( (latdiff**2) + ( londiff**2) )
        if(tmp < dmin):
            dmin = tmp
            dmindex = ijk
    xflag3 = False
    if dmin < 5:
        xflag3 = True
    return xflag3, dmindex


def nonZeroRandom():
    # (0, 1)
    x = random()
    while x == 0:
        x = random()
    return x


def splitScatter(infile):
    x_axis = []
    y_axis = []
    with open(infile, 'r') as f:
        for ijk in f:
            y, x = map(float, ijk.split(","))
            x_axis.append(x)
            y_axis.append(y)
    return x_axis, y_axis


def plotFigure(infile1, infile2, title, save, output):
    MHI_X, MHI_Y = splitScatter(infile1)
    HIFCZ_X, HIFCZ_Y = splitScatter(infile2)
    fig_x = MHI_X + HIFCZ_X
    fig_y = MHI_Y + HIFCZ_Y
    colors = ['b']*len(MHI_X) + ['r']*len(HIFCZ_X)

    fig = plt.figure(1)
    splot1 = plt.subplot(211, title = title)
    splot1.axes.get_xaxis().set_visible(False)
    splot1.axes.get_yaxis().set_visible(False)
    # http://stackoverflow.com/questions/2176424/hiding-axis-text-in-matplotlib-plots
    # if we want to allow labels
    #splot1.axes.get_xaxis().set_ticks([])
    #splot1.axes.get_yaxis().set_ticks([])
    plt.scatter(fig_x, fig_y, s = 5, lw = 0, c = colors,  edgecolor = colors, marker = '.')
    if(save):
        print "Saved plot as : '%s'"%( "%s.pdf"%(output) )
        plt.savefig("%s.pdf"%(output),  format = "pdf")
    plt.show()


###########
###########


def main(args):
    """
    Main body of the converted program.  For the most part, it should be a direct translation of the basic program.
    We have gone ahead and tried to make as much as possible non-static.  This will hopefully make it simpler 
    to modify the parameters, without having to modifying the code.
    """

    uval = dict([ (y, [ [0.0]*(251 + 1) for x in xrange(401 + 1)], ) for y in xrange(args.arraystart, args.arrayend + 1) ])
    vval = dict([ (y, [ [0.0]*(251 + 1) for x in xrange(401 + 1)], ) for y in xrange(args.arraystart, args.arrayend + 1) ])
    habilon = [0.0]*(302 + 1) # 302
    habilat = [0.0]*(302 + 1) # 302
    island = [0.0]*(302 + 1) # 302 

    print "Reading coastal fringe location file..."
    with open (args.locale) as file3:
        for line, ijk in zip(file3, xrange(1, 302 + 1)):
            habilon[ijk], habilat[ijk], island[ijk] = map(float, line.split(","))

    print "Reading and plotting Island data and EEZ data..."
    plotFigure(args.island, args.eez, "Island & EEZ Data", args.plot, args.output)

    myFile = {}
    print "Reading daily HYCOM currents..."
    with open(args.currents) as file1:
        for line, ijk in zip(file1, xrange(122, 1095 + 1)):
            #line = file1.readline()
            xfile = line.split(",")[0]
            if arraystart <= ijk <= arrayend:
                myFile[ijk] = xfile

    for ijk in xrange(args.arraystart, args.arrayend + 1):
        ufile = os.path.join(args.data, "%s_evel.dat"%(myFile[ijk]))
        vfile = os.path.join(args.data, "%s_nvel.dat"%(myFile[ijk]))
        # DLS
        # old program would try to use non-existant files if allowed.
        # not sure how it should really be handled, so for now print an error message,
        # but continue forward.
        if(not exists(ufile)):
            print >> sys.stderr, "FILE ERROR: %s does not exist"%(ufile)
            continue
        if(not exists(vfile)):
            print >> sys.stderr, "FILE ERROR: %s does not exist"%(vfile)
            continue 
        with open(ufile) as file3:
            with open(vfile) as file4:
                for j in xrange(251, 0, -1):
                    for i in xrange(1, 401 + 1):
                        uval[ijk][i][j] = float(file3.readline())
                        vval[ijk][i][j] =  float(file4.readline())
    del myFile

    with open(args.output, 'w') as file2:
        for startsite in xrange(1, 302 + 1):
            print startsite
            startlon = habilon[startsite]
            startlat = habilat[startsite]
            for startday in xrange(122, 580 + 1):
                if (122 <= startday <= 212) or (486 <= startday <= 577) or (851 <= startday <= 942):
                    for nrun in xrange(args.runs): 
                        endlon, endlat = release(startlon, startlat, startday, uval, vval, args.duration)
                        xflag3, dmindex = settle(endlon, endlat, uval, vval, habilon, habilat)
                        if xflag3 and dmindex:
                            print >> file2,  "%s\t%s\t%s\t%s\t%s" % (startsite, startday, dmindex, island[startsite], island[dmindex])


def parseArgs(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument("-l", "--locale", required = True, help ="Coastal fringe location file", type = str)
    parser.add_argument("-i", "--island", required = True, help ="Island data file", type = str)
    parser.add_argument("-e", "--eez", required = True, help ="EEZ data file", type = str)
    parser.add_argument("-c", "--currents", required = True, help ="Daily HYCOM currents file csv", type = str)
    parser.add_argument("-d", "--data", default = os.path.join(".","data"), help ="Directory path containing *_evel.dat and *_nvel.dat files", type = str)
    parser.add_argument("-o", "--output", required = True, help ="Output file", type = str)

    parser.add_argument("-t", "--duration", default = 45, help ="How many days", type = int)
    parser.add_argument("-s", "--arraystart", default = 122, help ="Start of analysis period", type = int)
    parser.add_argument("-n", "--arrayend", default = 626, help ="End of analysis period", type = int)
    parser.add_argument("-r", "--runs", default = 1000, help ="Originally known as ntotal", type = int)

    parser.add_argument("-p", "--plot", action = "store_true", help = "Save plot to pdf")
    args = parser.parse_args()
    return args
   

if __name__=="__main__":
    main(parseArgs(sys.argv))
