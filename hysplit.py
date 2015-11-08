"""
General Python class for creating HYSPLIT4 back trajectories.
Created on Sat Nov  1 21:11:10 2014

@author: Von P. Walden
         Washington State University
         Laboratory for Atmospheric Research
"""
class HYSPLIT4:
    """This class defines attributes and methods for calculating back trajectories
    		using NOAA's HYSPLIT4 model.
        
        Example usage #1:
            # Retrieve reanalysis data from NOAA for the ICECAPS time period.
            import hysplit
            hy = hysplit.HYSPLIT4('100501','140601')
            hy.retrieveReanalysisDataFromNOAA()
            
        Example usage #2:
            # Calculate a back trajectory for a given date.
            import hysplit
            hy = hysplit.HYSPLIT4('121218','121220')
            hy.runBackTrajectory(72.83, -38.46, [50, 500, 1000, 3000])  # Summit Station at 50 m, 500 m, 1000 m, and 3000 m.
            hy.plotBackTrajectory()

    """
    def __init__(self, beginningDate, endingDate):
        """INITIAL: Initializes the HYSPLIT processing object with the
                    desired dates to process.
                    
                    Written by Von P. Walden
                                 1 Nov 2014
        """
        import sys
        from socket import gethostname
        from datetime import datetime, timedelta
        from dateutil.relativedelta import relativedelta
        
        # Specify the beginning and ending dates to process.
#        beginningDate = datetime.strptime(sys.argv[1],'%y%m%d')
#       endingDate    = datetime.strptime(sys.argv[2],'%y%m%d')
        self.beginningDate = datetime.strptime(beginningDate,'%y%m%d')
        self.endingDate    = datetime.strptime(endingDate,   '%y%m%d')
        self.deltaHours    = timedelta(hours=6)
        self.deltaMonths   = relativedelta(months=1)
        if (self.beginningDate > self.endingDate):
            print('Beginning date is greater than the ending date.')
            sys.exit()
        
        # Set up the necessary directories.
        self.hostname = gethostname()
        if self.hostname.rfind('aeolus')>=0 or self.hostname.rfind('compute')>=0:
            # For aeolus.wsu.edu  (Linux)
            self.directory = {'data':  '/data/vonw/hysplit4/reanalysis/',
                              'traj':  '/data/vonw/hysplit4/backTrajectories/',
                              'code':  '/home/vonw/work/projects/icecaps/python/hysplit4/',
                              'plot':  '/data/vonw/hysplit4/backTrajectories/plots/'}
        elif self.hostname.rfind('sila')>=0 or self.hostname.rfind('nuia')>=0:
            # For sila.cee.wsu.edu  (iMac)
            self.directory = {'data':  '/Users/vonw/data/hysplit4/reanalysis/',
                              'traj':  '/Users/vonw/data/hysplit4/backTrajectories/',
                              'code':  '/Users/vonw/software/ICECAPS/backTrajectories/',
                              'plot':  '/Users/vonw/data/hysplit4/backTrajectories/plots/'}
        else:
            print('Your local computer is unrecognized by hysplit.HYSPLIT4.  Talk to Von!')
            sys.exit()
        
        return
    
    def retrieveReanalysisDataFromNOAA(self):
        """Automates the retrieval of reanalysis data from NOAAs
            data archive. Currently one can only retrieve data from the
            NCEP/NCAR reanalyses.
                    
                    Written by  Von P. Walden
                                 1 Nov 2014
        """
        import os
        from subprocess import call
        
        os.chdir(self.directory['data'])
        CDC = 'ftp://arlftp.arlhq.noaa.gov/pub/archives/reanalysis/'
        
        date  = self.beginningDate 
        while (date <= self.endingDate):
            f = 'RP'+date.strftime('%Y%m')+'.gbl'
            print('Retrieving NCEP/NCAR reanalyses data: '+f)
            call(['wget', CDC+f])
            
            date+=self.deltaMonths
        
        return
    
    def runBackTrajectory(self, latitutde, longitude, altitudes):
        """Creates HYSPLIT4 input file based on (lat,lon) and altitude, 
            then runs the back trajectory.  Note that "altitudes" can be a list
            of altitudes for running multiple altitudes for a given location.
                    
                    Written by  Von P. Walden
                                 1 Nov 2014
                    Updated     29 Sep 2014 - Added ability to specify
                                                lat, lon and alt.
        """
        import os
        from subprocess import call
        
        # Navigate to the "trajectory" directory.
        os.chdir(self.directory['traj'])
        
        date = self.beginningDate
        while (date <= self.endingDate):
            # Define a unique date string for the trajectory.
            dstr   = date.strftime('%Y%m%d%H')
            
            # Based on the desired date, determine the three nearest reanalysis data files.
            month1 = date.strftime('%Y%m')
            month2 = (date-self.deltaMonths).strftime('%Y%m')
            month3 = (date+self.deltaMonths).strftime('%Y%m')
            
            # Create the SETUP input file for HYSPLIT4
            f = open(self.directory['traj']+'SETUP.'+dstr,'w')
            f.write('&SETUP\n')
            f.write('KMSL=0,\n')
            f.write('tm_tpot=1,\n')
            f.write('tm_tamb=1,\n')
            f.write('tm_rain=1,\n')
            f.write('tm_mixd=1,\n')
            f.write('tm_relh=1,\n')
            f.write('tm_terr=1,\n')
            f.write('tm_dswf=1,\n')
            f.write('/ \n')
            f.close()
            
            # Create the CONTROL input file for HYSPLIT4
            f = open(self.directory['traj']+'CONTROL.'+dstr,'w')
            f.write(date.strftime('%y %m %d %H')+'\n')
            f.write('%d\n' % len(alt))
            for altitude in altitudes:
                f.write('%9f %10f %d\n' % (latitude, longtitude, altitude))
            f.write('-336\n')
            f.write('0	\n')
            f.write('13000\n')
            f.write('3\n')
            f.write(self.directory['data']+'\n')
            f.write('RP'+month1+'.gbl\n')
            f.write(self.directory['data']+'\n')
            f.write('RP'+month2+'.gbl\n')
            f.write(self.directory['data']+'\n')
            f.write('RP'+month3+'.gbl\n')
            f.write(self.directory['traj']+'\n')
            f.write('Summit'+date.strftime('%Y%m%d%H')+'.trj\n')
            f.close()
            
            if self.hostname.rfind('aeolus')>=0 or self.hostname.rfind('compute')>=0:
                # For aeolus.wsu.edu  (Linux)
                call(['/home/vonw/hysplit4/exec/hyts_std', dstr])
            elif self.hostname.rfind('sila')>=0 or self.hostname.rfind('nuia')>=0:
                # For sila.cee.wsu.edu  (iMac)
                call(['/Applications/Hysplit4/exec/hyts_std', dstr])
            
            date+=self.deltaHours
        
        return
    
    def plotBackTrajectory(self):
        """Generate a plot of HYSPLIT4 backtrajectories.
                    
                    Written by  Von P. Walden
                                 1 Nov 2014
        """
        from datetime import datetime
        from mpl_toolkits.basemap import Basemap
        import numpy as np
        import matplotlib.pyplot as plt
        
        date = self.beginningDate
        while (date <= self.endingDate):
            # Define a unique date string for the trajectory.
            dstr   = date.strftime('%Y%m%d%H')
            print('Creating plot for: Summit'+dstr+'.trj')
            
            # Read in the trajectory data file.
            backtraj = np.loadtxt(self.directory['traj']+'Summit'+dstr+'.trj',skiprows=9)
            
            # Store each altitude in a separate array.
            dn1  = []
            lat1 = np.array([])
            lon1 = np.array([]) 
            alt1 = np.array([]) 
            for line in backtraj[::3]:
                dn1.append(datetime(int(line[2])+2000,int(line[3]),int(line[4]),int(line[5]),int(line[5]),int(line[7])))
                lat1 = np.append(lat1, line[9])
                lon1 = np.append(lon1, line[10])
                alt1 = np.append(alt1, line[11])
            
            dn2  = []
            lat2 = np.array([])
            lon2 = np.array([]) 
            alt2 = np.array([]) 
            for line in backtraj[1::3]:
                dn2.append(datetime(int(line[2])+2000,int(line[3]),int(line[4]),int(line[5]),int(line[5]),int(line[7])))
                lat2 = np.append(lat2, line[9])
                lon2 = np.append(lon2, line[10])
                alt2 = np.append(alt2, line[11])
            
            dn3  = []
            lat3 = np.array([])
            lon3 = np.array([]) 
            alt3 = np.array([]) 
            for line in backtraj[2::3]:
                dn3.append(datetime(int(line[2])+2000,int(line[3]),int(line[4]),int(line[5]),int(line[5]),int(line[7])))
                lat3 = np.append(lat3, line[9])
                lon3 = np.append(lon3, line[10])
                alt3 = np.append(alt3, line[11])
            
            dn4  = []
            lat4 = np.array([])
            lon4 = np.array([]) 
            alt4 = np.array([]) 
            for line in backtraj[3::3]:
                dn4.append(datetime(int(line[2])+2000,int(line[3]),int(line[4]),int(line[5]),int(line[5]),int(line[7])))
                lat4 = np.append(lat4, line[9])
                lon4 = np.append(lon4, line[10])
                alt4 = np.append(alt4, line[11])
            
            # Create a nifty map.
            plt.figure()
            m = Basemap(resolution='l',projection='ortho',lat_0=90.,lon_0=0.)
            m.drawmapboundary()
            m.drawcoastlines()
            parallels = np.arange(-80.,90,10.)
            meridians = np.arange(0.,360.,10.)
            m.drawparallels(parallels)
            m.drawmeridians(meridians)
            
            x1, y1 = m(lon1,lat1)
            m.scatter(x1, y1, c=alt1/1000., linewidths=0, vmin=0., vmax=8., s=20)
            #m.plot(x1, y1, 'w', linewidth=0.7)
            plt.text(x1[-1],y1[-1],'50 meters',color='r')
            x2, y2 = m(lon2,lat2)
            m.scatter(x2, y2, c=alt2/1000., linewidths=0, vmin=0., vmax=8., s=20)
            #m.plot(x2, y2, 'w', linewidth=1)
            plt.text(x2[-1],y2[-1],'500 meters',color='r')
            x3, y3 = m(lon3,lat3)
            m.scatter(x3, y3, c=alt3/1000., linewidths=0, vmin=0., vmax=8., s=20)
            #m.plot(x3, y3, 'w', linewidth=1)
            plt.text(x3[-1],y3[-1],'1000 meters',color='r')
            x4, y4 = m(lon4,lat4)
            m.scatter(x4, y4, c=alt4/1000., linewidths=0, vmin=0., vmax=8., s=20)
            #m.plot(x3, y3, 'w', linewidth=1)
            plt.text(x3[-1],y3[-1],'3000 meters',color='r')
            cb=plt.colorbar(shrink=0.7)
            cb.set_label('Altitude (km)')
            plt.title('Back Trajectories from Summit, Greenland: '+dstr)
            
            plt.savefig(self.directory['plot']+'Summit'+dstr+'.png')
            plt.close('all')
            
            date+=self.deltaHours
        
        return
