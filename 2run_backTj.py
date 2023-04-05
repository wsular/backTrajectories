#Inputs for hysplit back trajectory 
# for Wenatchee and Yakima wineries
#%%
from datetime import datetime
import hysplit
import pandas as pd
import hvplot
import hvplot.pandas
import xarray as xr
from datetime import datetime, timedelta
import numpy as np

#%%


def addCoordsToACONC(date):
    from datetime import datetime, timedelta
    import xarray as xr

    yyyy = date.strftime('%Y')
    mm = date.strftime('%m')
    yyyymmdd = date.strftime('%Y%m%d')
    yyyymmddhh = date.strftime('%Y%m%d'+'00')
    grid_path = '/home/anacarla/airpact5_saved/grid/' 
    aconc3d = xr.open_dataset('/mnt/AIRNAS1/airpact/saved/' +
                            yyyy + '/' + mm + '/aconc/' + 'combined3d_' + yyyymmdd + '.ncf', engine='netcdf4' )
    aconc = xr.open_dataset('/mnt/AIRNAS1/airpact/saved/' +
                            yyyy + '/' + mm + '/aconc/' + 'combined_' + yyyymmdd + '.ncf', engine='netcdf4' )
    gridcro2d = xr.open_dataset(grid_path + 'GRIDCRO2D', engine='netcdf4')

    metcro3d = xr.open_dataset( '/mnt/AIRNAS1/airpact/AIRRUN/' + yyyy + '/' + yyyymmddhh + '/MCIP37/' + 'METCRO3D',  engine='netcdf4')
    re_aconc3d = xr.open_dataset('/mnt/AIRNAS1/airpact/AIRRUN/'+yyyy+'/'+ yyyymmdd + '00/CCTM/'+ 'ACONC_' + yyyymmdd +'.ncf', engine='netcdf4')

    # ....Assign coordinates to aconc file
    times = [date + timedelta(hours=int(hour)) for hour in aconc.TSTEP.values]
    ds3d = aconc3d.copy()
    ds3d = ds3d.drop('TFLAG')
    ds3d = ds3d.assign_coords({'time': (["TSTEP"], times),
                           'pressure': (["TSTEP", "LAY", "ROW", "COL"], metcro3d.PRES[:-1].values),
                           'latitude': (["ROW", "COL"], gridcro2d.LAT[0, 0].values),
                           'longitude': (["ROW", "COL"], gridcro2d.LON[0, 0].values)})
    re_ds3d = re_aconc3d.copy()
    re_ds3d = re_ds3d.drop('TFLAG')
    re_ds3d = re_ds3d.assign_coords({'time': (["TSTEP"], times),
                           'pressure': (["TSTEP", "LAY", "ROW", "COL"], metcro3d.PRES[:-1].values),
                           'latitude': (["ROW", "COL"], gridcro2d.LAT[0, 0].values),
                           'longitude': (["ROW", "COL"], gridcro2d.LON[0, 0].values)}) 
    ds = aconc.copy()
    ds = ds.drop('TFLAG')
    ds = ds.assign_coords({'time': (["TSTEP"], times),
                           'layer': (["LAY"], aconc.LAY[:].values),
                           'latitude': (["ROW", "COL"], gridcro2d.LAT[0, 0].values),
                           'longitude': (["ROW", "COL"], gridcro2d.LON[0, 0].values)})

    return ds,ds3d,re_ds3d


#%%
def find_WRF_pixel(latvar,lonvar,lat0,lon0):
    # Read latitude and longitude from file into numpy arrays
    # Renamed findWRFpixel from original function, naive_fast, written by Vikram Ravi.
    latvals = latvar[:]
    lonvals = lonvar[:]
    dist_sq = (latvals-lat0)**2 + (lonvals-lon0)**2
    minindex_flattened = dist_sq.argmin()  # 1D index of min element
    iy_min, ix_min = np.unravel_index(minindex_flattened, latvals.shape)

    return int(iy_min),int(ix_min)


#%%
def findACONCindex(time, pressure, lat, lon):
    # ....Find time
    it = np.where(aconc.time == np.datetime64(time))[0][0]
    print(it)
    # ....Find nearest gridcell
    ilat, ilon = find_WRF_pixel(aconc.latitude, aconc.longitude, lat, lon)
    # ....Find nearest pressure level
    dP = aconc3d.pressure[it, :, ilat, ilon].values
    #ipr = np.where(np.abs(dP) == np.abs(dP).min())[0][0]
    diff = np.subtract((np.abs(dP)/100), pressure)

    ipr = np.where(np.abs(diff) == np.amin(np.abs(diff)))[0][0]


    return it, ipr, ilat, ilon
#print(lats)
#%%
#%% Vertical Regriding following steps from: https://www.epa.gov/hesc/how-rsig-regrids-data



#%%
#Reading backtrajectories data
#path='/mnt/vonw/data/hysplit/backTrajectories/'
path= '/home/anacarla/hysplit/traj/'
#descriptor = 'Wenatcheev2'
descriptor = 'Yakima'
#dates  = pd.date_range('2020-09-10', '2020-09-15 00:00', freq='6H')
date = datetime(2020,9,13,00)
altitudes   = [50, 100, 500]
names=['trajectory','run','year','month','day','hour','minute','seconds','time','latitude','longitude','altitude','pressure','potential temperature','air temperature','rainfall','mix depth','relative humidity','terrain above msl','solar flux']
#for date in dates:
dstr= date.strftime('%Y%m%d%H')
       #print(dstr)
numalts  = len(altitudes)
df = pd.read_csv(path+descriptor+dstr+'.trj',skiprows=6+numalts,delimiter=r"\s+",names= names)
df['year'] = '20' + df['year'].astype(str)
df['Datetime']= pd.to_datetime(df.rename(columns={'day':'Day'})[['Day','month','year','hour']])
df.set_index(['Datetime','trajectory'], inplace=True)
df['AP_O3']=np.nan
df['AP_PM']=np.nan
df['AP_Pres']=np.nan
df['RRAP_O3'] = np.nan
df['RRAP_PM'] = np.nan
df['RRAP_Pres'] = np.nan
df['RRAP_OH'] = np.nan
df['RRAP_NO3']=np.nan ##nitrate ppmV
df['AP_ALD2']=np.nan ##acetaldehyde MW:44	Explicit
df['AP_CO']=np.nan ##carbon monoxide MW:28	Explicit
df['AP_ETH']=np.nan ##ethene MW:28	Explicit
df['AP_FORM']=np.nan ##formaldehyde MW:30	Explicit
df['AP_H2O2']=np.nan ##hydrogen peroxide MW:34	Explicit
df['AP_HNO3']=np.nan ##nitric acid MW:63	Explicit
df['AP_HONO']=np.nan ##nitrous acid MW:47	Explicit
df['AP_ISOP']=np.nan ##isoprene MW:68.1	Explicit
df['AP_IOLE']=np.nan ##internal alkene bond MW:56.1	Lumped
df['AP_N2O5']=np.nan ##dinitrogen pentoxide MW:108	Explicit
df['AP_NH3']=np.nan ##ammonia MW:17	Explicit
df['AP_NO']=np.nan ##nitric oxide MW:30	Explicit
df['AP_NO2']=np.nan ##nitrogen dioxide MW:46	Explicit
df['AP_NOX']=np.nan ##1000.0*(NO[1]+NO2[1] in ppbV 
df['AP_ANO3_PPB']=np.nan ##Aeorosol Nitrate MW:62	Explicit
df['AP_PAN']=np.nan ##peroxyacylnitrate MW:121	Explicit
df['AP_PANX']=np.nan ##peroxyacylnitrate with 3 or more carbons MW:135	Lumped
df['AP_SO2']=np.nan ##sulfur dioxide MW:64	Explicit
df['AP_SULF']=np.nan ##sulfuric acid (gaseous) MW:98	Explicit
df['AP_TERP']=np.nan ##monoterpenes MW:136.2	Lumped
df['AP_TOL']=np.nan ##toluene and other monoalkyl aromatics MW:92.1	Lumped
df['AP_VOC']=np.nan ##VOC:1000.0*(PAR[1]+2.0*ETH[1]+2.0*ETOH[1]+2.0*OLE[1]+7.0*TOL[1]+8.0*XYL[1]+FORM[1]+2
df['AP_XLY']=np.nan ##XYLENE
#%%
yak_df=pd.DataFrame()
dates=pd.date_range(start='09/08/2020', end='09/13/2020')
yak_df['Datetime']=dates
yak_o3=[]
yak_no3=[]
yak_oh=[]
yak_pm=[]
for d in dates:
    aconc,aconc3d,re_aconc3d = addCoordsToACONC(d.to_pydatetime())
    it, ipr, ilat, ilon = findACONCindex(d.to_pydatetime(), 930, 46.600, -120.500)
    print(it,ipr,ilat,ilon)
    yak_o3.append(aconc3d.O3[it,0,ilat,ilon].values)
    yak_no3.append(re_aconc3d.NO3[it,0,ilat,ilon].values)
    yak_oh.append(re_aconc3d.OH[it,0,ilat,ilon].values)
    yak_pm.append(aconc3d.PMIJ[it,0,ilat,ilon].values)
yak_df['AP_O3']=yak_o3
yak_df['AP_NO3']=yak_no3
yak_df['AP_OH']=yak_oh
yak_df['AP_PM']=yak_pm
yak_df['Intgr_NO3']=yak_df['AP_NO3'].cumsum()
yak_df['Intgr_OH']=yak_df['AP_OH'].cumsum()



#%%
for ind,rows in df.iterrows():
    #for i in range(0,len(df.i)):
    #if rows['trajectory'] == 1:
    print("Taking the conc values for:", ind[0].to_pydatetime(),"for trajecotry: ", ind[1])
    aconc,aconc3d,re_aconc3d = addCoordsToACONC(ind[0].to_pydatetime())
    it, ipr, ilat, ilon = findACONCindex(ind[0].to_pydatetime(), rows['pressure'], rows['latitude'], rows['longitude'])
    print(ipr)
    df.loc[ind,'AP_O3'] = aconc3d.O3[it,ipr,ilat,ilon].values
    df.loc[ind,'AP_PM'] = aconc3d.PMIJ[it,ipr,ilat,ilon].values 
    df.loc[ind,'AP_Pres']=aconc3d.pressure[it,ipr,ilat,ilon].values
    df.loc[ind,'RRAP_O3'] = re_aconc3d.O3[it,ipr,ilat,ilon].values
    df.loc[ind,'RRAP_NO3'] = re_aconc3d.NO3[it,ipr,ilat,ilon].values
    #df.loc[ind,'RRAP_PM'] = re_aconc3d.PMIJ[it,ipr,ilat,ilon].values 
    df.loc[ind,'RRAP_Pres']=re_aconc3d.pressure[it,ipr,ilat,ilon].values 
    df.loc[ind,'RRAP_OH']=re_aconc3d.OH[it,ipr,ilat,ilon].values
    df.loc[ind,'AP_ALD2']=aconc.ALD2[it,0,ilat,ilon].values
    df.loc[ind,'AP_CO']=aconc.CO[it,0,ilat,ilon].values
    df.loc[ind,'AP_ETH']=aconc.ETH[it,0,ilat,ilon].values
    df.loc[ind,'AP_FORM']=aconc.FORM[it,0,ilat,ilon].values
    df.loc[ind,'AP_H2O2']=aconc.H2O2[it,0,ilat,ilon].values
    df.loc[ind,'AP_HNO3']=aconc.HNO3[it,0,ilat,ilon].values
    df.loc[ind,'AP_HONO']=aconc.HONO[it,0,ilat,ilon].values
    df.loc[ind,'AP_ISOP']=aconc.ISOP[it,0,ilat,ilon].values
    df.loc[ind,'AP_IOLE']=aconc.IOLE[it,0,ilat,ilon].values
    df.loc[ind,'AP_N2O5']=aconc.N2O5[it,0,ilat,ilon].values
    df.loc[ind,'AP_NH3']=aconc.NH3[it,0,ilat,ilon].values
    df.loc[ind,'AP_NO']=aconc.NO[it,0,ilat,ilon].values
    df.loc[ind,'AP_NO2']=aconc.NO2[it,0,ilat,ilon].values
    df.loc[ind,'AP_NOX']=aconc.NOX[it,0,ilat,ilon].values
    df.loc[ind,'AP_ANO3']=aconc.ANO3_PPB[it,0,ilat,ilon].values
    df.loc[ind,'AP_PAN']=aconc.PAN[it,0,ilat,ilon].values
    df.loc[ind,'AP_PANX']=aconc.PANX[it,0,ilat,ilon].values
    df.loc[ind,'AP_SO2']=aconc.SO2[it,0,ilat,ilon].values
    df.loc[ind,'AP_SULF']=aconc.SULF[it,0,ilat,ilon].values
    df.loc[ind,'AP_TERP']=aconc.TERP[it,0,ilat,ilon].values
    df.loc[ind,'AP_TOL']=aconc.TOL[it,0,ilat,ilon].values
    df.loc[ind,'AP_VOC']=aconc.VOC[it,0,ilat,ilon].values
    df.loc[ind,'AP_XYL']=aconc.XYL[it,0,ilat,ilon].values
#%%
df.to_csv("/home/anacarla/hysplit/traj/dframe_2020_ap_hysp_Yakima.csv")
#%% Plotting
#from bokeh.models.formatters import DatetimeTickFormatter

pm=df.hvplot(x='Datetime', y=['AP_PM'], value_label=('ug/m^3'),logy=True, by='trajectory', legend=True)
o3=df.hvplot(x='Datetime', y=['AP_O3'], value_label=('ppb'), by='trajectory', legend=True)
oh = df.hvplot(x='Datetime', y=['RRAP_OH'], value_label=('ppmV'), by='trajectory', legend=True,logy=True) 
alt=df.hvplot(x='Datetime', y=['altitude'], value_label=('Meter'), by='trajectory', legend=True)
no3=df.hvplot(x='Datetime', y=['RRAP_NO3'], value_label=('ppmV'), by='trajectory', legend=True, logy=True)
#%%
#df_surf=df.reset_index()
#df_surf=df.drop_duplicates(subset=['AP_VOC','AP_TOL','AP_NOX','AP_ISOP','AP_ANO3'],keep='first')
df_surf=df[['AP_VOC','AP_TOL','AP_NOX','AP_ISOP','AP_ANO3','AP_SO2','RRAP_OH','RRAP_NO3']]
#df.iloc[df.index.get_level_values('A') == 1]
df_surf=df_surf.loc[df_surf.index.get_level_values('trajectory')==1]
gases=df_surf.hvplot(x='Datetime', y=['AP_VOC','AP_NOX','RRAP_OH','RRAP_NO3'], value_label=('ppbV'),logy=True)
#%%
cum_no3=yak_df.hvplot(x='Datetime',y=['Intgr_NO3'],value_label=('ppmV'),title='Integrated NO3',logy=True)
cum_oh=yak_df.hvplot(x='Datetime',y=['Intgr_OH'],value_label=('ppmV'),title='Integrated OH')
yak_o3_plot=yak_df.hvplot(x='Datetime', y=['AP_O3'], value_label=('ppb'), legend=True)
yak_pm_plot=yak_df.hvplot(x='Datetime', y=['AP_PM'], value_label='ug/m^3', legend=True)
#%%
df.hvplot(x='Datetime', y=['AP_Pres'], value_label=('Pa'))
#df.hvplot(x='Datetime', y = ['AP_PM', 'AP_O3'], value_label=('ug/m^3','ppbV'),subplots = True, shared_axes=False)
#%%
traj=df.hvplot('latitude','longitude',
                        geo=True,
                        c='trajectory',
                        cmap='Category10',
                        tiles='StamenTerrain',
                        title='Back Trajecotry for ' + descriptor + ' winery ' + ': ' + dstr)
# %%
pm_geo= df.hvplot.points('longitude', 'latitude', 
                             geo=True, 
                             c=df['AP_PM'], 
                             clabel= 'AP PM (ug/m^3)',
                             cmap='viridis', 
                             logz=True,
                             tiles='StamenTerrain', 
                             title='Back Trajectories for ' + descriptor + ' (' + 'Airpact PM' + ')' + ': ' + dstr)
hvplot.save(pm_geo, '/home/anacarla/hysplit/traj/plots/' + descriptor + '_' + 'AP_PM' + '_' + dstr + '.html')
# %%
o3_geo = df.hvplot.points('longitude', 'latitude', 
                             geo=True, 
                             c=df['AP_O3'], 
                             clabel= 'AP O3 (ppb)',
                             cmap='viridis', 
                             tiles='StamenTerrain', 
                             marker='circle',
                             title='Back Trajectories for ' + descriptor + ' (' + 'Airpact O3' + ')' + ': ' + dstr)

# %%
oh_geo= df.hvplot.points('longitude', 'latitude', 
                             geo=True, 
                             c=df['RRAP_OH'], 
                             clabel= 'AP OH (ppm)',
                             cmap='viridis', 
                             logz=True,
                             tiles='StamenTerrain', 
                             title='Back Trajectories for ' + descriptor + ' (' + 'Airpact OH' + ')' + ': ' + dstr)
# %%
no3_geo= df.hvplot.points('longitude', 'latitude', 
                             geo=True, 
                             c=df['RRAP_NO3'], 
                             clabel= 'AP NO3 (ppm)',
                             cmap='viridis', 
                             logz=True,
                             tiles='StamenTerrain', 
                             title='Back Trajectories for ' + descriptor + ' (' + 'Airpact NO3' + ')' + ': ' + dstr)

# %%
