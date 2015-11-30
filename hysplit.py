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
            from datetime import datetime
            dates = [datetime(2010,5,1), datetime(2014,6,1)]
            hy = hysplit.HYSPLIT4(dates)
            hy.retrieveReanalysisDataFromNOAA()
            
        Example usage #2:
            # Calculate a back trajectory for a given date.
            import hysplit
            from datetime import datetime
            dates  = [datetime(2015,1,1,0,0), datetime(2015,1,1,6,0), datetime(2015,1,1,12,0), datetime(2015,1,1,18,0)]
            length = 336
            lat    = +72.83
            lon    = -38.46
            alts   = [50, 500, 1000, 3000]
            hy = hysplit.HYSPLIT4(dates,length,lat,lon,alts,'Summit')
            hy.runBackTrajectory()  # Summit Station at 50 m, 500 m, 1000 m, and 3000 m.
            hy.plotBackTrajectory()

    """
    def __init__(self, dates, length = 48, latitude=0., longitude=0., altitudes=0., descriptor='backTrajectory'):
        """INITIAL: Initializes the HYSPLIT processing object with the
                    desired dates to process.

            Inputs:
                dates      - list of datetime objects to create back trajectories for.
                length     - length of back trajectory in hours.
                latitude   - latitude of the start position for back trajectory.
                longitude  - longitude of the start position for back trajectory.
                altitude   - series of altitudes to start back trajectories from.
                descriptor - text that describes back trajectories; used for filenames.

                    Written by Von P. Walden
                                 1 Nov 2014
                    Updated on   8 Nov 2015
                                20 Nov 2015 - Changed the way dates are inputed; list of
                                                datetime objects instead of beginning and ending.
                                29 Nov 2015 - Can now input latitude, longitude and altitudes.
        """
        import sys
        from socket import gethostname
        from datetime import datetime, timedelta
        from dateutil.relativedelta import relativedelta
        
        # Specify the dates to process.
        self.dates         = dates
        self.length        = length
        self.latitude      = latitude
        self.longitude     = longitude
        self.altitudes     = altitudes
        self.descriptor    = descriptor
        self.deltaMonths   = timedelta(months=1)
        
        # Set up the necessary directories.
        self.hostname = gethostname()
        if self.hostname.rfind('aeolus')>=0 or self.hostname.rfind('compute')>=0:
            # For aeolus.wsu.edu  (Linux)
            self.directory = {'data':  '/data/vonw/hysplit4/reanalysis/',
                              'traj':  '/data/vonw/hysplit4/backTrajectories/',
                              'code':  '/home/vonw/work/software/backtrajectories/',
                              'plot':  '/data/vonw/hysplit4/backTrajectories/plots/'}
        elif self.hostname.rfind('sila')>=0 or self.hostname.rfind('nuia')>=0:
            # For sila.cee.wsu.edu  (iMac)
            self.directory = {'data':  '/Users/vonw/data/hysplit4/reanalysis/',
                              'traj':  '/Users/vonw/data/hysplit4/backTrajectories/',
                              'code':  '/Users/vonw/software/backTrajectories/',
                              'plot':  '/Users/vonw/data/hysplit4/backTrajectories/plots/'}
        else:
            print('Whoa, Partner... Your local computer is unrecognized by hysplit.HYSPLIT4.  Talk to Von!')
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
        
        for date in dates:
            f = 'RP'+date.strftime('%Y%m')+'.gbl'
            print('Retrieving NCEP/NCAR reanalyses data: '+f)
            call(['wget', CDC+f])
        
        return
    
    def runBackTrajectory(self):
        """Creates HYSPLIT4 input file based on the __init__ function.
                    
                    Written by  Von P. Walden
                                 1 Nov 2014
                    Updated     29 Sep 2014 - Added ability to specify
                                                lat, lon and alt.
                                20 Nov 2015 - Added capability to specify the 
                                                number of hours to calculate the
                                                back trajectory for.
                                29 Nov 2015 - Cleaned things up by putting all
                                                initializations into __init__.
        """
        import os
        from subprocess import call
        
        # Navigate to the "trajectory" directory.
        os.chdir(self.directory['traj'])
        
        for date in dates:
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
            f.write(date.strftime('%y %m %d %H %M')+'\n')
            f.write('%d\n' % len(alt))
            for altitude in altitudes:
                f.write('%9f %10f %d\n' % (latitude, longitude, altitude))
            f.write(str(int(-length))+'\n')
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
            f.write(self.descriptor+date.strftime('%Y%m%d%H')+'.trj\n')
            f.close()
            
            if self.hostname.rfind('aeolus')>=0 or self.hostname.rfind('compute')>=0:
                # For aeolus.wsu.edu  (Linux)
                call(['/home/vonw/hysplit4/exec/hyts_std', dstr])
            elif self.hostname.rfind('sila')>=0 or self.hostname.rfind('nuia')>=0:
                # For sila.cee.wsu.edu  (iMac)
                call(['/Applications/Hysplit4/exec/hyts_std', dstr])
        
        return
    
    def plotBackTrajectory(self):
        """Generate a plot of HYSPLIT4 back trajectories.
                    
                    Written by  Von P. Walden
                                 1 Nov 2014
        """
        from datetime import datetime
        from mpl_toolkits.basemap import Basemap
        import numpy as np
        import matplotlib.pyplot as plt
        
        for date in dates:
            # Define a unique date string for the trajectory.
            dstr   = date.strftime('%Y%m%d%H')
            print('Creating plot for: '+self.descriptor+dstr+'.trj')
            
            # Read in the trajectory data file.
            numalts  = len(self.altitudes)
            backtraj = np.loadtxt(self.directory['traj']+self.descriptor+dstr+'.trj',skiprows=6+numalts)

            # Create the back trajectory plot.
            plt.figure()
            m = Basemap(resolution='l',projection='ortho',lat_0=90.,lon_0=0.)
            m.drawmapboundary()
            m.drawcoastlines()
            parallels = np.arange(-80.,90,10.)
            meridians = np.arange(0.,360.,10.)
            m.drawparallels(parallels)
            m.drawmeridians(meridians)
            
            # Plot each altitude.
            for ind in range(numalts):
                lat = backtraj[ind::numalts,9]
                lon = backtraj[ind::numalts,10]
                alt = backtraj[ind::numalts,11]

                x, y = m(lon,lat)
                m.scatter(x, y, c=alt/1000., linewidths=0, vmin=0., vmax=8., s=20)
                plt.text(x[-1],y[-1],str(altitudes[ind])+' meters',color='r')
            
            cb=plt.colorbar(shrink=0.7)
            cb.set_label('Altitude (km)')
            plt.title('Back Trajectories for '+self.descriptor+': '+dstr)
            
            plt.savefig(self.directory['plot']+self.descriptor+dstr+'.png')
            plt.close('all')
        
        return
