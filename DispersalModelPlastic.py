#!/usr/bin/env python

#matplotlib needs numpy, and possibly the dev files for freetype
#from meliae import scanner
#scanner.dump_all_objects("memmap.out" ) 


#import matplotlib.pyplot as plt
import sys
import os
from os.path import exists
from math import sqrt, log, atan, cos, sin
from random import seed, random
import sqlite3
import argparse
from datetime import datetime, timedelta
#import dateutil.parser
#from tmepfile import NamedTemporaryFile
from multiprocessing import Pool


######################################
######################################
### Fixed variables, do not change ###
######################################
######################################
TPI = 8 * atan(1)
#TIME_CONVERSION = 60 * 60 * 24
diffco = 0.00175
sqrt_diffco = sqrt(diffco)
EVEL_NVEL_ERR = -99
######################################
######################################


#### RELEASE FUNCTION STARTS THE TRANSPORT SIMULATION PROCESS#### 
def release(args, xlon, xlat, startstep, con, habilon, habilat, startsite, all_output, startdate, timeslice):
    xflag1 = False
    xflag2 = False
    xflag3 = False
    currentdate = startdate
    currentstep = 0
    dmindex = None
    prevxlon = None
    prevxlat = None
    '''add if loop if smaller than 51 days, keep going, if larger then check with Settle after every day.'''
    for step in xrange(1, args.maxPLD + 1):
        currentdate = currentdate + timeslice
        if step < args.minPLD:
            xflag3 = False
        else: #elif step >= args.minPLD:
            # DLS - Settle does not change if we have the same lon/lat... 
            # simple storage of the previous lat/lon should reduce unneeded computation calls.  Small test shows about a 50% call reduction.
            if xlon != prevxlon or xlat != prevxlat:
                xflag3, dmindex = settle(args, xlon, xlat, habilon, habilat)
                prevxlon = xlon
                prevxlat = xlat
        if xflag3 == False:
            if xflag1 == False:
                currentstep = startstep + step
                xlon, xlat = disperse(args, xlon, xlat, currentstep, con)
                xflag1 = checkbounds(args, xlon, xlat)
                if xflag1 == False:
                    xflag2 = checkdepth(args, xlon, xlat, currentstep, con)
            if xflag2 == True:
                xlon, xlat = tweak(args, xlon, xlat, con, currentstep) 
        
        print >> all_output,  "%s\t%s\t%s\t%s\t%s" % (startsite, step, currentdate.isoformat(), xlon, xlat)
    return xflag3, dmindex, xlon, xlat, currentdate


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


def getVU(con, step, lon, lat):
    """
    For a given step, lat, lon, query the datbase for the v and u values (REAL)
    """
    return con.execute("""SELECT v, u FROM uvvals WHERE step=? AND lon=? AND lat=? LIMIT 1""", (step, lon, lat,)).fetchone()


#### DISPERSE FUNCTION TAKES A PARTICLE POSITION AND MOVES IT SPATIALLY 1 TIME STEP ####
def disperse(args, xlon, xlat, currentstep, con):
    hycomu = 0
    hycomv = 0
    zloncm, zlatcm = calculateLonLatCM(abs(xlat))
    #space
    nnlon, nnlat = computeNNlonlat(xlon, args.lonmin, xlat, args.latmin, args.resolution)
    if (1 <= nnlon <= args.londim) and (1 <= nnlat <= args.latdim):
        hycomv, hycomu = getVU(con, currentstep, nnlon, nnlat)
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
    degu = zloncm * args.TIME_CONVERSION * hycomu
    degv = zlatcm * args.TIME_CONVERSION * hycomv
    degu += XX1 * sqrt_diffco
    degv += XX2 * sqrt_diffco
    xlon += degu
    xlat += degv
    return xlon, xlat


#### FUNCTION CHECKING FOR BEING OUT OF BOUNDS AND RETURNS FLAG ####
def checkbounds(args, xlon, xlat):
    if xlat < args.latmin or xlat > args.latmax or xlon < args.lonmin or xlon > args.lonmax:
        return  True
    return False


#### FUNCTION CHECKING DEPTH OF PARTICLE AND RETURNS FLAG ####
def checkdepth(args, xlon, xlat, currentstep, con):
    nnlon, nnlat = computeNNlonlat(xlon, args.lonmin, xlat, args.latmin, args.resolution)
    if 1 <= nnlon <= args.londim and 1 <= nnlat <= args.latdim:
        return EVEL_NVEL_ERR in  getVU(con, currentstep, nnlon, nnlat)
    return False


#### TWEAK FUNCTION TAKES A PARTICLE POSITION AND ADJUSTS IT RANDOMLY ####
# This function used to be recursive in TrueBasic, but is now itterative in python
def tweak(args, xlon, xlat, con, currentstep):
    newxlon = 0.0
    newxlat = 0.0
    zloncm, zlatcm = calculateLonLatCM(abs(xlat))
    deguconst = zloncm  * args.TIME_CONVERSION
    degvconst = zlatcm  * args.TIME_CONVERSION
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
        xflag2 = checkdepth(args, newxlon, newxlat, currentstep, con)
    return newxlon, newxlat


#### SETTLE FUNCTION CHECKS FOR PROXIMITY TO DEFINED HABITAT ####
def settle(args, endlon, endlat, habilon, habilat):
    dmin = 6  
    dmindex = None

    # testing the idea of local access is quicker than global..
    absf = abs
    sqrtf = sqrt
    xrangef = xrange
    releasesites =  args.releasesites
    computeKMLatLonf =  computeKMLatLon

    # la1 = -0.000003437984
    # la2 = 0.0004652059
    # la3 = -0.001507974
    # la4 = 110.572
    # lo1 = 0.00006610779
    # lo2 = -0.02041057
    # lo3 = 0.064812
    # lo4 = 111.1078

    for ijk in xrangef(1, releasesites + 1):
        hlat = habilat[ijk]
        hlon = habilon[ijk]
        avglat = absf((endlat + hlat) / 2.0)

        # avglat2 = avglat**2
        # avglat3 = avglat**3
        # zkmlat = (la1 * avglat3) + (la2 * avglat2) + (la3 * avglat) + la4
        # zkmlon = (lo1 * avglat3) + (lo2 * avglat2) + (lo3 * avglat) + lo4        
        # DLS - The call overhead can be brutal.. Need to look into this more.. 
        # inlining would be ideal, but that would req use of something like pypy
        zkmlat, zkmlon = computeKMLatLonf(avglat)
        latdiff = (hlat - endlat) * zkmlat
        londiff = (hlon - endlon) * zkmlon
        tmp = sqrtf( (latdiff**2) + ( londiff**2) )
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
    """                                                                         
    Remove any database that matches the one we plan to use, and then create a new database. 
    In the case, we do not have a recent version of sqlite, we need to fall back to a version that has a ROWID on the rows. 
    """

    # remove and old sqlite database, since we dont want it around
    removeDatabase(args)
    con = sqlite3.connect(args.database_loc, check_same_thread = False)
    try:
        con.execute("""CREATE TABLE IF NOT EXISTS uvvals( step INTEGER NOT NULL, 
                                                      lat INTEGER NOT NULL, 
                                                      lon INTEGER NOT NULL, 
                                                      u REAL NOT NULL, 
                                                      v REAL NOT NULL, 
                                                      PRIMARY KEY(step, lat, lon)) WITHOUT ROWID;""")     
    except:
        print >> sys.stderr, "sqlite version is older than 3.8.2.  Using tables with rowid"
        con.execute("""CREATE TABLE IF NOT EXISTS uvvals( step INTEGER NOT NULL, 
                                                      lat INTEGER NOT NULL, 
                                                      lon INTEGER NOT NULL, 
                                                      u REAL NOT NULL, 
                                                      v REAL NOT NULL,
                                                      PRIMARY KEY(step, lat, lon));""")     
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
    # dangerous in the case we care about the database in that it could get corrupt with power loss. 
    con.execute("""PRAGMA synchronous = OFF;""")

    #### READING IN CURRENT SOLUTIONS ####
    # This reads in the 365 days of HYCOM currents from Yanli as parsed to ascii by Evan and regridded by Don
    # Note that ALL array indexing is now on a julian style from 1=1/1/2009
    # Simple count index starts at 1 for 5/2/2009 (122 julian) 
    startdatestamp = datetime(year = 2009, month = 1, day = 1)
    step = 0
    startDate = None
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
            if not startDate:
                startDate = startdatestamp + timedelta(days = (julian - 1))
            cdate = startdatestamp + timedelta(days = (julian - 1))
            cdate = cdate.isoformat()
            print cdate
            with open(ufile, "r") as ufilep:
                with open(vfile, "r") as vfilep:
                    vudat = []
                    #for j in xrange(latdim, 0, -1):
                    # prepare for a bulk insert into the database
                    for j in xrange(1, args.latdim + 1): 
                        for i in xrange(1, args.londim + 1):
                            vudat.append( (step, j, i, float(vfilep.readline()), float(ufilep.readline())) )
                    # need to verify the practical limits of executemany.. in which case we might need to fall back to plain execute()
                    con.executemany("""INSERT INTO uvvals(step, lat, lon, v, u) VALUES(?, ?, ?, ?, ?)""", vudat)
            
            step += args.total_steps
    if args.total_steps != 1:
        print "Begin Interpolation"
        con.executemany("""INSERT INTO uvvals(step, lat, lon, v, u) VALUES(?, ?, ?, ?, ?)""", stepGenerator(con, args, step) )
        print "\nEnd Interpolation"
    con.commit()
    args.arrayend = ((args.arrayend - args.arraystart) * args.total_steps) + 1
    args.arraystart = 0
    #con.execute("""VACUUM;""")    
    con.execute("""PRAGMA query_only = 1;""")
    con.execute("""PRAGMA synchronous = OFF;""")
    con.execute("""PRAGMA journal_mode = OFF;""")
    con.execute("""PRAGMA cache_size = 20000;""")
    con.execute("""PRAGMA locking_mode = EXCLUSIVE;""")
    con.commit()
    print "Done."
    return startDate

def stepGenerator(con, args, step):
    """ 
    A generator function that we use to interpolate the current changes during a day in a given time step. 
    """

    cur = con.cursor()
    for i, j in zip(xrange(0, step, args.total_steps), xrange(args.total_steps, step, args.total_steps)):
        sys.stderr.write(".")
        for row in list(cur.execute("""SELECT a.lat as clat, b.lon as clon, a.v AS cv, a.u AS cu, b.v AS nv, b.u AS nu FROM uvvals AS a JOIN uvvals AS b ON (a.lat = b.lat AND a.lon = b.lon) WHERE a.step =? AND b.step = ?;""", (i, j,))):
            lat = row[0]
            lon = row[1]
            cv = row[2]
            cu = row[3]
            nv = row[4]
            nu = row[5]
            vstep = ( (nv - cv) / float(args.total_steps))
            ustep = ( (nu - cu) / float(args.total_steps))
            timeslice = timedelta(hours = 24.0 / float(args.total_steps) )
            for x in xrange(1, args.total_steps):                
                cv += vstep
                cu += ustep
                yield (i+x, lat, lon, cv, cu, )


#@profile
def main():
    args = parseArgs(sys.argv)
    if args.rng_seed != None:
        seed(args.rng_seed)

    #### SETTING DIMENSIONS FOR THE SIMULATION I THINK ####
    # DLS - attempted to add this to the sqliteDB, but didn't really help speed, and the amount of memory savings is minimal..
    habilon = [0.0]*(args.releasesites + 1)
    habilat = [0.0]*(args.releasesites + 1)
    island = [0.0]*(args.releasesites + 1)  ## DLS - Not used?
    propreef = [0.0]*(args.releasesites + 1) ## DLS - Not used?

    #### LOADING SETTLEMENT AND RELEASE HABITAT FROM FILE ####
    print "Reading habitat file..."
    with open (args.reefs_csv, "r") as file3:
        for ijk in xrange(1, args.releasesites + 1):
            habilat[ijk], habilon[ijk], propreef[ijk], island[ijk] = map(float, file3.readline().rstrip().split(","))
            #habilat[ijk], habilon[ijk], propreef[ijk], island[ijk] = map(float, file3.readline().rstrip().split(","))

    '''Commented this out because I don't need it right now, but would like for it to stay in the code for future use'''
    #print "Reading Island data and EEZ data..."  
    #plotFigure("MHI.csv", "HIFCZ.csv","Island & EEZ Data", "output")

    # create a new DB with tables and return the connector
    con = createDatabase(args)
    # read in all our input files into the database
    offset = readinUandVfiles(con, args)
    
    #### STARTS THE DIFFERENT FUCTIONS AND THE TRANSPORT SIMULATION ####
    # This is tha part that actually initiates all the subroutines and initiates the simulation
    # It is also responsible for opening output files and closing them after the simulation is done
    start = args.arraystart
    end = args.arrayend - args.maxPLD #+ 1
    timeslice = timedelta(hours = 24.0 / float(args.total_steps) )
    for startsite in xrange(1, args.releasesites + 1):
        outputfile_endpoint = None
        if args.output_prefix_total:
            outputfile_endpoint = args.output_prefix_total + str(startsite) + ".txt" # Total settlemetn file. Comment out if only want daily.
        outputfile = args.output_prefix_all + str(startsite) + ".txt"
        with open(outputfile, 'w') as all_output:
            totout = None
            if outputfile_endpoint:
                totout = open(outputfile_endpoint, "w")

            startlon = habilon[startsite]
            startlat = habilat[startsite]

            #stime = datetime.now()
            startstamp = offset
            for step in xrange(end):
                for nrun in xrange(args.ntotal):
                    xflag3, dmindex, xlon, xlat, currentdate = release(args, startlon, startlat, step, con, habilon, habilat, startsite, all_output, startstamp, timeslice)
                    if xflag3 and dmindex:
                        print >> all_output,  "%s\t%s\t%s\t%s\t%s" % (startsite, dmindex, currentdate.isoformat(), xlon, xlat)
                        if totout:
                            print >> totout,  "%s\t%s\t%s\t%s\t%s" % (startsite, startstamp, dmindex, island[startsite], island[dmindex])
                    if xflag3 == False:
                        print >> all_output,  "%s\t%s\t%s\t%s\t%s" % (startsite, 0, currentdate.isoformat(), xlon, xlat)
                startstamp = startstamp + timeslice
            sys.stderr.write(".")
            #etime = datetime.now()
            #print etime - stime
            if totout:
                totout.close()
    print "\nDone"
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

    parser.add_argument("-o", "--output_prefix_all", required = True, help ="all timestep output path prefix (example: ./out/Dispersal_python_8kmHA_test_site_ )", type = str)
    parser.add_argument("-e", "--output_prefix_total", default = None, help ="total output path prefix (example: ./out/Dispersal_python_8kmHA_test_site_total )", type = str)

    parser.add_argument("-t", "--minPLD", default = 2, help ="min PLD in steps", type = int)
    parser.add_argument("-m", "--maxPLD", default = 5, help ="max PLD in steps", type = int)

    parser.add_argument("-s", "--arraystart", default = 122, help ="Start of analysis period (days)", type = int)
    parser.add_argument("-n", "--arrayend", default = 626, help ="End of analysis period (days)", type = int)

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

    parser.add_argument("-i", "--total_steps", default = 1, help ="How many steps in a given day >=1", type = int) # delete/comment
    #parser.add_argument("-e", "--eez", required = True, help ="EEZ data file", type = str) # delete/commnt
    #parser.add_argument("-p", "--plot", action = "store_true", help = "Save plot to pdf") #delete/comment
    parser.set_defaults(TIME_CONVERSION=60*60*24)
    args = parser.parse_args()
    args.TIME_CONVERSION = 60.0 * 60.0 * (24.0/args.total_steps)
    print args.TIME_CONVERSION
    return args


if __name__ == '__main__':
    # just to make sure we randomize the random function
    main()
