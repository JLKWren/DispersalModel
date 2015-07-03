#!/usr/bin/env python

#matplotlib needs numpy, and possibly the dev files for freetype
#import matplotlib.pyplot as plt
import sys
import os
from os.path import exists
from math import sqrt, log, atan, cos, sin
from random import seed, random
import sqlite3
import argparse
from datetime import datetime


######################################
######################################
### Fixed variables, do not change ###
######################################
######################################
TPI = 8 * atan(1)
TIME_CONVERSION = 60 * 60 * 24
diffco = 0.00175
sqrt_diffco = sqrt(diffco)
EVEL_NVEL_ERR = -99
######################################
######################################


#### RELEASE FUNCTION STARTS THE TRANSPORT SIMULATION PROCESS#### 
def release(args, xlon, xlat, startday, con, habilon, habilat, startsite, daily):
    xflag1 = False
    xflag2 = False
    xflag3 = False
    currentday = 0
    dmindex = None
    prevxlon = None
    prevxlat = None
    '''add if loop if smaller than 51 days, keep going, if larger then check with Settle after every day.'''
    for day in xrange(1, args.maxPLD + 1):
        if day < args.minPLD:
            xflag3 = False
        else: #elif day >= args.minPLD:
            # DLS - Settle does not change if we have the same lon/lat... 
            # simple storage of the previous lat/lon should reduce unneeded computation calls for test show about a 50% call reduction
            if xlon != prevxlon or xlat != prevxlat:
                xflag3, dmindex = settle(args, xlon, xlat, habilon, habilat)
                prevxlon = xlon
                prevxlat = xlat
        if xflag3 == False:
            if xflag1 == False:
                currentday = startday + day
                xlon, xlat = disperse(args, xlon, xlat, currentday, con)
                xflag1 = checkbounds(args, xlon, xlat)
                if xflag1 == False:
                    xflag2 = checkdepth(args, xlon, xlat, currentday, con)
            if xflag2 == True:
                xlon, xlat = tweak(args, xlon, xlat, con, currentday) 
        print >> daily,  "%s\t%s\t%s\t%s\t%s" % (startsite, day, currentday, xlon, xlat)
    return xflag3, dmindex, xlon, xlat, currentday


def computeKMLatLon(avglat):
    """
    Pull out the calculation for the zkmlat and zkmlon since it is used in at least 2 locations.
    Function call overhead isn't really an issue yet, so manual inlining is not required
    """
    la1 = -0.000003437984
    la2 = 0.0004652059
    la3 = -0.001507974
    la4 = 110.572
    lo1 = 0.00006610779
    lo2 = -0.02041057
    lo3 = 0.064812
    lo4 = 111.1078
    avglat2 = avglat**2
    avglat3 = avglat**3
    zkmlat = (la1 * avglat3) + (la2 * avglat2) + (la3 * avglat) + la4
    zkmlon = (lo1 * avglat3) + (lo2 * avglat2) + (lo3 * avglat) + lo4
    return zkmlat, zkmlon


def computeNNlonlat(xlon, lonmin, xlat, latmin, resolution):
    """
    Much like the zkmlat,zkmlon calculation, this set of code is used in at lesat two locations.
    Pulled into a function to make it simpler to modify if needed.
    """
    nnlon = int(round( ((xlon - lonmin) / resolution), 0)) + 1
    nnlat = int(round( ((xlat - latmin) / resolution), 0))  + 1
    return nnlon, nnlat


#### FUNCTION CALCULATING LAT AND LONG ####
def calculateLonLatCM(avglat):
    zkmlat, zkmlon = computeKMLatLon(avglat)
    zlatkm = 1.0 / zkmlat
    zlonkm = 1.0 / zkmlon
    zlatcm = zlatkm / 100000.0
    zloncm = zlonkm / 100000.0
    return zloncm, zlatcm


def getVU(con, day, lon, lat):
    """
    For a given day, lat, lon, query the datbase for the v and u values (REAL)
    """
    return con.execute("""SELECT v, u FROM uvvals WHERE day=? AND lon=? AND lat=? LIMIT 1""", (day, lon, lat,)).fetchone()


#### DISPERSE FUNCTION TAKES A PARTICLE POSITION AND MOVES IT SPATIALLY 1 TIME STEP ####
def disperse(args, xlon, xlat, currentday, con):
    hycomu = 0
    hycomv = 0
    zloncm, zlatcm = calculateLonLatCM(abs(xlat))
    #space
    nnlon, nnlat = computeNNlonlat(xlon, args.lonmin, xlat, args.latmin, args.resolution)
    if (1 <= nnlon <= args.londim) and (1 <= nnlat <= args.latdim):
        hycomv, hycomu = getVU(con, currentday, nnlon, nnlat)
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
def checkbounds(args, xlon, xlat):
    #(args.latmin <= xlat <= args.latmax) or (args.lonmin <= xlon <= args.lonmax)
    if xlat < args.latmin or xlat > args.latmax or xlon < args.lonmin or xlon > args.lonmax:
        return  True
    return False


#### FUNCTION CHECKING DEPTH OF PARTICLE AND RETURNS FLAG ####
def checkdepth(args, xlon, xlat, currentday, con):
    nnlon, nnlat = computeNNlonlat(xlon, args.lonmin, xlat, args.latmin, args.resolution)
    if 1 <= nnlon <= args.londim and 1 <= nnlat <= args.latdim:
        return EVEL_NVEL_ERR in  getVU(con, currentday, nnlon, nnlat)
    return False


#### TWEAK FUNCTION TAKES A PARTICLE POSITION AND ADJUSTS IT RANDOMLY ####
# This function used to be recursive in TrueBasic, but is now itterative in python
def tweak(args, xlon, xlat, con, currentday):
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
        xflag2 = checkdepth(args, newxlon, newxlat, currentday, con)
    return newxlon, newxlat


#### SETTLE FUNCTION CHECKS FOR PROXIMITY TO DEFINED HABITAT ####
def settle(args, endlon, endlat, habilon, habilat):
    dmin = 6  
    dmindex = None

    # la1 = -0.000003437984
    # la2 = 0.0004652059
    # la3 = -0.001507974
    # la4 = 110.572
    # lo1 = 0.00006610779
    # lo2 = -0.02041057
    # lo3 = 0.064812
    # lo4 = 111.1078

    for ijk in xrange(1, args.releasesites + 1):
        hlat = habilat[ijk]
        hlon = habilon[ijk]
        avglat = abs((endlat + hlat) / 2.0)

        # avglat2 = avglat**2
        # avglat3 = avglat**3
        # zkmlat = (la1 * avglat3) + (la2 * avglat2) + (la3 * avglat) + la4
        # zkmlon = (lo1 * avglat3) + (lo2 * avglat2) + (lo3 * avglat) + lo4        
        # DLS - The call overhead can be brutal.. Need to look into this more.. 
        # inlining would be ideal, but that would req use of something like pypy
        zkmlat, zkmlon = computeKMLatLon(avglat)
        latdiff = (hlat - endlat) * zkmlat
        londiff = (hlon - endlon) * zkmlon
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


# #### FUNCTION ####
# def splitScatter(infile):
#     x_axis = []
#     y_axis = []
#     with open(infile, 'r') as f:
#         for ijk in f:
#             y, x = ijk.rstrip().split(",")
#             x_axis.append(x)
#             y_axis.append(y)
#     return x_axis, y_axis


# #### FUNCTION PLOTTING FIGURE OF HAWAII AND THE HIFCZ ####
# # this function is not necessary and can be commented out without problem
# def plotFigure(infile1, infile2, title, output):
#     MHI_X, MHI_Y = splitScatter(infile1)
#     HIFCZ_X, HIFCZ_Y = splitScatter(infile2)
#     fig_x = MHI_X + HIFCZ_X
#     fig_y = MHI_Y + HIFCZ_Y
#     colors = ['b']*len(MHI_X) + ['r']*len(HIFCZ_X)
#     fig = plt.figure(1)
#     splot1 = plt.subplot(211, title = title)
#     splot1.axes.get_xaxis().set_visible(False)
#     splot1.axes.get_yaxis().set_visible(False)
#     # http://stackoverflow.com/questions/2176424/hiding-axis-text-in-matplotlib-plots
#     # if we want to allow labels
#     #splot1.axes.get_xaxis().set_ticks([])
#     #splot1.axes.get_yaxis().set_ticks([])
#     plt.scatter(fig_x, fig_y, s = 5, lw = 0, c = colors,  edgecolor = colors, marker = '.')
#     plt.savefig("%s.pdf"%(output),  format = "pdf")
#     plt.show()


def getUandVFile(pth, prefix, u_suffix, v_suffix):
    """
    create a path to the u and v files.  Also check that they exist.  If they don't throw out a warning, and return None so that a parent
    function can deal with the problem as it deems fit.
    """

    ufile = os.path.join(pth, prefix + u_suffix)
    vfile = os.path.join(pth, prefix + v_suffix)
    if(not exists(ufile)):
        print >> sys.stderr, "FILE ERROR: %s does not exist"%(ufile)
        ufile = None
    if(not exists(vfile)):
        print >> sys.stderr, "FILE ERROR: %s does not exist"%(vfile)
        vfile = None
    return ufile, vfile


def createDatabase(args):
    # remove and old sqlite database, since we dont want it around
    removeDatabase(args)
    con = sqlite3.connect(args.database_loc, check_same_thread = False)
    try:
        con.execute("""CREATE TABLE IF NOT EXISTS uvvals( day INTEGER NOT NULL, 
                                                      lat INTEGER NOT NULL, 
                                                      lon INTEGER NOT NULL, 
                                                      u REAL NOT NULL, 
                                                      v REAL NOT NULL, 
                                                      PRIMARY KEY(day, lat, lon)) WITHOUT ROWID;""")     
    except:
        print >> sys.stderr, "sqlite version is older than 3.8.2.  Using tables with rowid"
        con.execute("""CREATE TABLE IF NOT EXISTS uvvals( day INTEGER NOT NULL, 
                                                      lat INTEGER NOT NULL, 
                                                      lon INTEGER NOT NULL, 
                                                      u REAL NOT NULL, 
                                                      v REAL NOT NULL, 
                                                      PRIMARY KEY(day, lat, lon));""")     

    # dangerous in the case we care about the database in that it could get corrupt with power loss. 
    con.execute("""PRAGMA synchronous = OFF;""")
    con.commit()
    return con

def removeDatabase(args):
    # remove and old sqlite database, since we dont want it around  
    try:
        os.remove(args.database_loc)
    except:
        pass


def readinUandVfiles(con, args):
    """
    The currents_csv file already contains the julian date value.  This means we do not need to hard code
    a look that covers whatever range we desire.  Instead, we can read it line by line, not storing it in memory
    and filter to only the range that we want, using the arraystart/arrayend values compared to the juliand date
    """
    print "Reading daily HYCOM currents..."
    with open(args.currents_csv, "r") as currentFile:
        for line in currentFile:
            xfile, index, julian = line.split(",")
            julian = int(julian)
            if args.arraystart > julian or args.arrayend < julian:
                continue

            ufile, vfile = getUandVFile(args.currents_data_dir, xfile, args.currents_u_suffix, args.currents_v_suffix)
            # If the ufile or vfile doesn't not exist, we can't really continue adding it to the database..
            # Therefore, skip over it and move on to the next file.
            if not ufile or not vfile:
                continue

            print julian
            with open(ufile, "r") as ufilep:
                with open(vfile, "r") as vfilep:
                    vudat = []
                    #for j in xrange(latdim, 0, -1):
                    # prepare for a bulk insert into the database
                    for j in xrange(1, args.latdim + 1): 
                        for i in xrange(1, args.londim + 1):
                            vudat.append( (julian, j, i, float(vfilep.readline()), float(ufilep.readline()), ) )
                    # need to verify the practical limits of executemany.. in which case we might need to fall back to plain execute()
                    con.executemany("""INSERT INTO uvvals(day, lat, lon, v, u) VALUES(?, ?, ?, ?,?)""", vudat)
    con.commit()
    print "Done."


#@profile
def main():
    args = parseArgs(sys.argv)
    if args.rng_seed != None:
        seed(args.rng_seed)

    #### SETTING DIMENSIONS FOR THE SIMULATION I THINK ####
    # DLS - attempted to add this to the sqliteDB, but didn't really help speed, and the amount of memory savings is minimal..
    habilon = [0.0]*(args.releasesites + 1)
    habilat = [0.0]*(args.releasesites + 1)
    #island = [0.0]*(args.releasesites + 1)  ## DLS - Not used?
    #propreef = [0.0]*(args.releasesites + 1) ## DLS - Not used?

    #### LOADING SETTLEMENT AND RELEASE HABITAT FROM FILE ####
    print "Reading habitat file..."
    with open (args.reefs_csv, "r") as file3:
        for ijk in xrange(1, args.releasesites + 1):
            habilat[ijk], habilon[ijk], propreef, island = map(float, file3.readline().rstrip().split(","))
            #habilat[ijk], habilon[ijk], propreef[ijk], island[ijk] = map(float, file3.readline().rstrip().split(","))
            

    '''Commented this out because I don't need it right now, but would like for it to stay in the code for future use'''
    #print "Reading Island data and EEZ data..."  
    #plotFigure("MHI.csv", "HIFCZ.csv","Island & EEZ Data", "output")

    # create a new DB with tables and return the connector
    con = createDatabase(args)
    # read in all our input files into the database
    readinUandVfiles(con, args)
    con.execute("""PRAGMA synchronous = ON;""")
    con.commit()
    #### STARTS THE DIFFERENT FUCTIONS AND THE TRANSPORT SIMULATION ####
    # This is tha part that actually initiates all the subroutines and initiates the simulation
    # It is also responsible for opening output files and closing them after the simulation is done
    start = args.arraystart
    end = args.arrayend - args.maxPLD + 1
    for startsite in xrange(1, args.releasesites + 1):
        outputfile_endpoint = None
        if args.output_prefix_total:
            outputfile_endpoint = args.output_prefix_total + str(startsite) + ".txt" # Total settlemetn file. Comment out if only want daily.
        outputfile = args.output_prefix_daily + str(startsite) + ".txt"
        with open(outputfile, 'w') as daily:
            totout = None
            if outputfile_endpoint:
                totout = open(outputfile_endpoint, "w")
            #stime = datetime.now()
            print startsite
            startlon = habilon[startsite]
            startlat = habilat[startsite]
            for startday in xrange(start, end):
                for nrun in xrange(args.ntotal): 
                    xflag3, dmindex, xlon, xlat, currentday = release(args, startlon, startlat, startday, con, habilon, habilat, startsite, daily)
                    if xflag3 and dmindex:
                        print >> daily,  "%s\t%s\t%s\t%s\t%s" % (startsite, dmindex, currentday, xlon, xlat)
                    if xflag3 == False:
                        print >> daily,  "%s\t%s\t%s\t%s\t%s" % (startsite, 0, currentday, xlon, xlat)
            #etime = datetime.now()
            #print etime - stime
            if totout:
                totout.close()
    con.close()
    removeDatabase(args)


def parseArgs(argv):
    """
    commandline parser.  Should cover all the previously static variables that required changes.  
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument("-l", "--reefs_csv", required = True, help ="Habitat location file", type = str) 

    parser.add_argument("-c", "--currents_csv", required = True, help ="Daily HYCOM currents file csv", type = str)
    parser.add_argument("-d", "--currents_data_dir", default = os.path.join(".","data"), help ="Directory path containing *_u and *_v files", type = str)
    parser.add_argument("-u", "--currents_u_suffix", required = True , help ="u file suffix", type = str)
    parser.add_argument("-v", "--currents_v_suffix", required = True, help ="v file suffix", type = str)

    parser.add_argument("-o", "--output_prefix_daily", required = True, help ="daily output path prefix (example: ./out/Dispersal_python_8kmHA_test_site_ )", type = str)
    parser.add_argument("-e", "--output_prefix_total", default = None, help ="total output path prefix (example: ./out/Dispersal_python_8kmHA_test_site_total )", type = str)

    parser.add_argument("-t", "--minPLD", default = 2, help ="PLD in days", type = int)
    parser.add_argument("-m", "--maxPLD", default = 5, help ="PLD in days", type = int)

    parser.add_argument("-s", "--arraystart", default = 122, help ="Start of analysis period", type = int)
    parser.add_argument("-n", "--arrayend", default = 626, help ="End of analysis period", type = int)

    parser.add_argument("-r", "--ntotal", default = 100, help ="Particles released per site per day", type = int)
    parser.add_argument("-z", "--releasesites", default = 687, help ="Number of habitat pixels", type = int)

    parser.add_argument("-a", "--londim", required = True, help ="Longitude dimension", type = int)
    parser.add_argument("-b", "--latdim", required = True, help ="Latitude dimension", type = int)

    parser.add_argument("-f", "--lonmin", required = True, help ="Longitude min", type = int)
    parser.add_argument("-g", "--lonmax", required = True, help ="Longitude max", type = int)
    parser.add_argument("-j", "--latmin", required = True, help ="Latitude min", type = int)
    parser.add_argument("-k", "--latmax", required = True, help ="Latitude max", type = int)
    parser.add_argument("-q", "--resolution", required = True, help ="Resolution of current file in km", type = float)

    parser.add_argument("-w", "--database_loc", default = os.path.join(".", "model.db3"), help ="Path at which to build the sqlite database.", type = str)    
    parser.add_argument("-x", "--rng_seed", default = None, help ="Set a static seed value for random", type = int)    

    #parser.add_argument("-i", "--island_csv", required = True, help ="Island data file", type = str) # delete/comment
    #parser.add_argument("-e", "--eez", required = True, help ="EEZ data file", type = str) # delete/commnt
    #parser.add_argument("-p", "--plot", action = "store_true", help = "Save plot to pdf") #delete/comment

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # just to make sure we randomize the random function
    main()
