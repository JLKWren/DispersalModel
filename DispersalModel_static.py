#!/usr/bin/python

#matplotlib needs numpy, and possibly the dev files for freetype
import matplotlib.pyplot as plt
import sys
import os
from os.path import exists
from math import sqrt, log, atan, cos, sin
from random import seed, random


# just to make sure we randomize the random function
seed()


# GLOBAL STATIC VALUES 
# Set these for all runns
reefs_csv = "HIreefsNew.csv"
currents_csv = "files_HA8km.csv"
currents_u_suffix = "_100m_u.txt"
currents_v_suffix = "_100m_v.txt"
currents_data_dir = "hawaii8km_/" #"hawaii8km_/"
output_prefix_daily = "Dispersal_python_8kmHA_test_site_"
output_prefix_total = "Dispersal_python_8kmHA_test_total.txt"
#PLD = (15,30,45,60,75)  # This is for looping over PLD's. Comment out if only running one PLD 
duration = 5
arraystart = 122
arrayend = 130
ntotal = 5
releasesites = 687
resolution= 0.08
londim = 438
latdim = 251
lonmin = 175
lonmax = 210
latmin = 15
latmax = 35

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

#----------------------------------- START OF SUBROUTINES ---------------------------------------

#### RELEASE FUNCTION STARTS THE TRANSPORT SIMULATION PROCESS #### 
def release(xlon, xlat, startday, uval, vval):
    xflag1 = False
    xflag2 = False
    global currentday
    for day in xrange(1, duration + 1):
        if xflag1 == False:
            currentday = startday + day
            xlon, xlat = disperse(xlon, xlat, currentday, uval, vval)
            xflag1 = checkbounds(xlon, xlat)
            if day < duration and xflag1 == False:
                xflag2 = checkdepth(xlon, xlat, currentday, uval, vval)
        if xflag2 == True:
            xlon, xlat = tweak(xlon, xlat, uval, vval, currentday)
        print >> file5,  "%s\t%s\t%s\t%s\t%s" % (startsite, day, currentday, xlon, xlat)
    return xlon, xlat


#### FUNCTION CALCULATING LAT AND LONG ####
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


#### DISPERSE FUNCTION TAKES A PARTICLE POSITION AND MOVES IT SPATIALLY 1 TIME STEP ####
def disperse(xlon, xlat, currentday, uval, vval):
    hycomu = 0
    hycomv = 0
    zloncm, zlatcm = calculateLonLatCM(abs(xlat))
    #space
    nnlon = int(round( ((xlon - lonmin) / resolution), 0)) + 1
    nnlat = int(round( ((xlat - latmin) / resolution), 0))  + 1
    if (1 <= nnlon <= londim) and (1 <= nnlat <= latdim):
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
    #space
    degu = zloncm * TIME_CONVERSION * hycomu
    degv = zlatcm * TIME_CONVERSION * hycomv
    degu += XX1 * sqrt_diffco
    degv += XX2 * sqrt_diffco
    xlon += degu
    xlat += degv
    return xlon, xlat


#### FUNCTION CHECKING FOR BEING OUT OF BOUNDS AND RETURNS FLAG ####
def checkbounds(xlon, xlat):
    if xlat < latmin or xlat > latmax or xlon < lonmin or xlon > lonmax:
        return  True
    return False


#### FUNCTION CHECKING DEPTH OF PARTICLE AND RETURNS FLAG ####
def checkdepth(xlon, xlat, currentday, uval, vval):
    nnlon = int(round( ((xlon - lonmin) / resolution), 0)) + 1
    nnlat = int(round( ((xlat - latmin) / resolution), 0)) + 1
    if 1 <= nnlon <= londim and 1 <= nnlat <= latdim:
        hycomu = uval[currentday][nnlon][nnlat]
        hycomv = vval[currentday][nnlon][nnlat]
        if hycomu == EVEL_NVEL_ERR or hycomv == EVEL_NVEL_ERR:
            #print >> file9,  "%s\t%s\t%s\t%s\t%s" % (xlon, xlat, hycomu, hycomv, startday)  # added this for trouble shooting. 
            #print xlon, xlat, uval, vval 
            #print "is this working?"
            return  True
    return False


#### TWEAK FUNCTION TAKES A PARTICLE POSITION AND ADJUSTS IT RANDOMLY ####
# This function used to be recursive in TrueBasic, but is now itterative in python
def tweak(xlon, xlat, uval, vval, currentday):
    newxlon = 0.0
    newxlat = 0.0
    zloncm, zlatcm = calculateLonLatCM(abs(xlat))
    deguconst = zloncm  * TIME_CONVERSION
    degvconst = zlatcm  * TIME_CONVERSION
    #space
    xflag2 = True   
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


#### SETTLE FUNCTION CHECKS FOR PROXIMITY TO DEFINED HABITAT ####
def settle(endlon, endlat, uval, vval, habilon, habilat):
    dmin = 6  
    dmindex = None
    for ijk in xrange(1, releasesites + 1):
        avglat = abs((endlat + habilat[ijk]) / 2.0)
        avglat2 = avglat**2
        avglat3 = avglat**3
        # space
        zkmlat = (la1 * avglat3) + (la2 * avglat2) + (la3 * avglat) + la4
        zkmlon = (lo1 * avglat3) + (lo2 * avglat2) + (lo3 * avglat) + lo4
        # space
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

#### RANDOMIZING FUNCTION ####
def nonZeroRandom():
    # (0, 1)
    x = random()
    while x == 0:
        x = random()
    return x

#### FUNCTION ####
def splitScatter(infile):
    x_axis = []
    y_axis = []
    with open(infile, 'r') as f:
        for ijk in f:
            y, x = ijk.rstrip().split(",")
            x_axis.append(x)
            y_axis.append(y)
    return x_axis, y_axis


#### FUNCTION PLOTTING FIGURE OF HAWAII AND THE HIFCZ ####
# this function is not necessary and can be commented out without problem
def plotFigure(infile1, infile2, title, output):
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
    plt.savefig("%s.pdf"%(output),  format = "pdf")
    plt.show()

#--------------------------------------------------- END OF SUBROUTINES ---------------------------------------



#### SETTING DIMENSIONS FOR THE SIMULATION I THINK ####
myFile = {}
uval = dict([ (y, [ [0.0]*(latdim + 1) for x in xrange(londim + 1)], ) for y in xrange(arraystart, arrayend + 1) ])
vval = dict([ (y, [ [0.0]*(latdim + 1) for x in xrange(londim + 1)], ) for y in xrange(arraystart, arrayend + 1) ])
habilon = [0.0]*(releasesites + 1) 
habilat = [0.0]*(releasesites + 1) 
island = [0.0]*(releasesites + 1)  
propreef = [0.0]*(releasesites + 1)

#### LOADING SETTLEMENT AND RELEASE HABITAT FROM FILE ####
print "Reading habitat file..."
with open (reefs_csv, "r") as file3:
    for ijk in xrange(1, releasesites + 1):
        habilat[ijk], habilon[ijk], propreef[ijk], island[ijk] = map(float, file3.readline().rstrip().split(","))

#print "Reading Island data and EEZ data..."
'''Commented this out because I don't need it right now, but would like for it to stay in the code for future use'''
#plotFigure("MHI.csv", "HIFCZ.csv","Island & EEZ Data", "output")

#### READING IN CURRENT SOLUTIONS ####
# This reads in the 365 days of HYCOM currents from Yanli as parsed to ascii by Evan and regridded by Don
# Note that ALL array indexing is now on a julian style from 1=1/1/2009
# Simple count index starts at 1 for 5/2/2009 (122 julian) 
print "Reading daily HYCOM currents..."
with open(currents_csv, "r") as file1:
    for ijk in xrange(121, 2038 + 1):
        line = file1.readline()
        xfile = line.split(",")[0]
        if arraystart <= ijk <= arrayend:
            myFile[ijk] = xfile

for ijk in xrange(arraystart, arrayend + 1):
    print ijk
    ufile = currents_data_dir + myFile[ijk] + currents_u_suffix
    vfile = currents_data_dir + myFile[ijk] + currents_v_suffix
    if(not exists(ufile)):
        print >> sys.stderr, "FILE ERROR: %s does not exist"%(ufile)
        continue
    if(not exists(vfile)):
        print >> sys.stderr, "FILE ERROR: %s does not exist"%(vfile)
        continue 
    with open(ufile, "r") as file3:
        with open(vfile, "r") as file4:
            #for j in xrange(latdim, 0, -1):
            for j in xrange(1, latdim + 1): 
                for i in xrange(1, londim + 1):
                    uval[ijk][i][j] = float(file3.readline())
                    vval[ijk][i][j] =  float(file4.readline())
    

#### STARTS THE DIFFERENT FUCTIONS AND THE TRANSPORT SIMULATION ####
# This is tha part that actually initiates all the subroutines and initiates the simulation
# It is also responsible for opening output files and closing them after the simulation is done
outputfile = output_prefix_total
with open(outputfile, 'w') as file2 :
    #for idx, duration in enumerate(PLD):
        for startsite in xrange(1, releasesites + 1):
            outputfile5= output_prefix_daily + str(startsite) + "_daily_1.txt"
            with open(outputfile5, 'w') as file5 :              
                print startsite
                #print duration
                startlon = habilon[startsite]
                startlat = habilat[startsite]
                for startday in xrange(arraystart, arrayend - duration + 1):
                    for nrun in xrange(ntotal): 
                        endlon, endlat = release(startlon, startlat, startday, uval, vval)
                        xflag3, dmindex = settle(endlon, endlat, uval, vval, habilon, habilat)
                        if xflag3 and dmindex:
                            test = 0
                            print >> file2,  "%s\t%s\t%s\t%s\t%s" % (startsite, startday, dmindex, island[startsite], island[dmindex])
                            print >> file5,  "%s\t%s\t%s\t%s\t%s" % (startsite, currentday, test, endlon, endlat)



