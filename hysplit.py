"""
General Python class for creating HYSPLIT back trajectories.
Created on Sat Nov  1 21:11:10 2014

@author: Von P. Walden
         Washington State University
         Laboratory for Atmospheric Research
"""
class HYSPLIT:
    """This class defines attributes and methods for calculating back trajectories
    		using NOAA's HYSPLIT model.
        
        Example usage #1:
            # Retrieve reanalysis data from NOAA for the ICECAPS time period.
            import hysplit
            import pandas as pd
            dates = pd.date_range('2019-12-01','2020-03-01', freq='1M')   # Dec 2019, Jan 2020 and Feb 2020
            hy = hysplit.HYSPLIT(dates)
            hy.retrieveReanalysisDataFromNOAA()
            
        Example usage #2:
            # Calculate a back trajectory for a given date.
            import hysplit
            import pandas as pd
            dates  = pd.date_range('2020-01-09', '2020-01-17 18:00', freq='3H')   # Every 3 hours for 9 days in Jan 2020  
            length = 120
            lat    = +46.7
            lon    = -117.2
            alts   = [500, 1000, 1500]
            hy = hysplit.HYSPLIT(dates,length,lat,lon,alts,'Pullman')
            hy.runBackTrajectory()  # Pullman, WA at 500 m, 1000 m and 1500 m.
    """
    def __init__(self, dates, length = 48, latitude=0., longitude=0., altitudes=0., descriptor='backTrajectory'):
        """INITIAL: Initializes the HYSPLIT processing object with the
                    desired dates to process.

            Inputs:
                dates      - list of datetime objects to create back trajectories for; MUST BE A LIST!!
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
        from dateutil.relativedelta import relativedelta
        
        # Specify the dates to process.
        self.dates         = dates
        self.length        = length
        self.latitude      = latitude
        self.longitude     = longitude
        self.altitudes     = altitudes
        self.descriptor    = descriptor
        self.deltaMonths   = relativedelta(months=1)
        
        # Set up the necessary directories.
        self.hostname = gethostname()
        if self.hostname.rfind('aeolus')>=0 or self.hostname.rfind('compute')>=0:
            # For aeolus.wsu.edu  (Linux)
            self.directory = {'data':  '/data/vonw/hysplit4/reanalysis/',
                              'traj':  '/data/vonw/hysplit4/backTrajectories/',
                              'code':  '/home/vonw/hysplit4/exec/',
                              'plot':  '/data/vonw/hysplit4/backTrajectories/plots/'}
        elif self.hostname.rfind('sila')>=0 or self.hostname.rfind('nuia')>=0:
            # For sila.paccar.wsu.edu  (iMac)
            self.directory = {'data':  '/Users/vonw/data/hysplit/reanalysis/',
                              'traj':  '/Users/vonw/data/hysplit/backTrajectories/',
                              'code':  '/Users/vonw/hysplit/exec/',
                              'plot':  '/Users/vonw/data/hysplit/backTrajectories/plots/'}
        else:
            print('Whoa, Partner... Your local computer is unrecognized by hysplit.HYSPLIT.  Talk to Von!')
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
        CDC = 'ftp://arlftp.arlhq.noaa.gov/archives/reanalysis/'
        
        for date in self.dates:
            f = 'RP'+date.strftime('%Y%m')+'.gbl'
            print('Retrieving NCEP/NCAR reanalyses data: '+f)
            call(['wget', CDC+f])
        
        return
    
    def retrieveMostRecent7daysReanalysisDataFromNOAA(self):
        """Automates the retrieval of reanalysis data from NOAAs
            data archive. Currently one can only retrieve data from the
            NCEP/NCAR reanalyses.

        Usage:
            import hysplit
            dates = [datetime.datetime.now()]
            hy    = hysplit.HYSPLIT(dates)
            hy.retrieveMostRecent7daysReanalysisDataFromNOAA()
                    
                    Written by  Von P. Walden
                                 1 Nov 2014
        """
        import os
        from subprocess import call
        
        os.chdir(self.directory['data'])
        CDC = 'ftp://arlftp.arlhq.noaa.gov/pub/archives/gdas1/'
        
        for date in self.dates:
            f = 'current7days.t00z'
            print('Retrieving GDAS reanalyses data: '+f)
            call(['wget', CDC+f])
        
        return
    
    def runBackTrajectory(self):
        """Creates HYSPLIT input file based on the __init__ function, then
            runs the specified backtrajectory.
                    
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
        
        for date in self.dates:
            # Define a unique date string for the trajectory.
            dstr   = date.strftime('%Y%m%d%H')
            
            # Based on the desired date, determine the three nearest reanalysis data files.
            month1 = date.strftime('%Y%m')
            month2 = (date-self.deltaMonths).strftime('%Y%m')
            month3 = (date+self.deltaMonths).strftime('%Y%m')
            
            # Create the SETUP input file for HYSPLIT
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
            
            # Create the CONTROL input file for HYSPLIT
            f = open(self.directory['traj']+'CONTROL.'+dstr,'w')
            f.write(date.strftime('%y %m %d %H %M')+'\n')
            f.write('%d\n' % len(self.altitudes))
            for altitude in self.altitudes:
                f.write('%9f %10f %d\n' % (self.latitude, self.longitude, altitude))
            f.write(str(int(-self.length))+'\n')
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
                call([self.directory['code']+'hyts_std', dstr])
            elif self.hostname.rfind('sila')>=0 or self.hostname.rfind('nuia')>=0:
                # For sila.cee.wsu.edu  (iMac)
                call([self.directory['code']+'hyts_std', dstr])
        
        return
    
    def plotBackTrajectory(self, colorScaleVariable='altitude'):
        """Generate a plot of HYSPLIT back trajectories.        
                    Written by  Von P. Walden
                                 1 Nov 2014
                    Updated:    11 Aug 2020 - Now uses cartopy and geopandas instead of basemap.
        """
        import pandas as pd
        import hvplot
        import hvplot.pandas  # noqa
        
        for date in self.dates:
            # Define a unique date string for the trajectory.
            dstr   = date.strftime('%Y%m%d%H')
            print('Creating HTML plot for: '+self.descriptor+dstr+'.trj')
            
            # Read in the trajectory data file.
            numalts  = len(self.altitudes)
            df = pd.read_csv(self.directory['traj']+self.descriptor+dstr+'.trj', 
                             skiprows=6+numalts, 
                             delimiter=r"\s+", 
                             names=['trajectory', 
                                    'run', 
                                    'year', 
                                    'month', 
                                    'day', 
                                    'hour', 
                                    'minute', 
                                    'seconds', 
                                    'time', 
                                    'latitude', 
                                    'longitude', 
                                    'altitude', 
                                    'pressure', 
                                    'potential temperature', 
                                    'air temperature', 
                                    'rainfall', 
                                    'mix depth', 
                                    'relative humidity', 
                                    'terrain above msl', 
                                    'solar flux'])
            
            plot = df.hvplot.points('longitude', 'latitude', 
                             geo=True, 
                             c=df[colorScaleVariable], 
                             clabel=colorScaleVariable,
                             cmap='viridis', 
                             tiles='CartoEco', 
                             title='Back Trajectories for ' + self.descriptor + ' (' + colorScaleVariable + ')' + ': ' + dstr)
            
            hvplot.save(plot, self.directory['plot'] + self.descriptor + '_' + colorScaleVariable + '_' + dstr + '.html')
        
        return
