import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import pytz
import hysplit

lance = Dataset('10minute.nc')
lat   = lance.variables['latitude'][...]
lon   = lance.variables['longitude'][...]
utime = lance.variables['unix_time'][...]
time  = np.array([datetime.fromtimestamp(t,tz=pytz.UTC) for t in utime])

stime = np.array([datetime(2000+int(f.strip()[-10:-8]),int(f.strip()[-8:-6]),int(f.strip()[-6:-4]),int(f.strip()[-2:]),0,0,0,pytz.UTC) for f in open('N-ICEsondeList.dat')])

for t in stime[:10]:
	ind = np.min(np.where(t<=time))
	print(t,lat[ind],lon[ind])

	hy = hysplit.HYSPLIT4([t],24*7,lat[ind],lon[ind],[0,22,50,100],'N-ICE_nearSfc_')
	hy.runBackTrajectory()
	hy.plotBackTrajectory()
