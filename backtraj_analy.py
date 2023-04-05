#%%
from datetime import datetime
import hysplit
import pandas as pd
import hvplot
import hvplot.pandas
import xarray as xr
from datetime import datetime, timedelta
import numpy as np
import folium
from folium.plugins import TimestampedGeoJson, HeatMap, MarkerCluster
from matplotlib import cm
from branca.colormap import StepColormap
from branca.colormap import LinearColormap
from matplotlib.dates import DateFormatter
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
#Reading backtrajectories data
#path='/mnt/vonw/data/hysplit/backTrajectories/'
path= '/home/anacarla/hysplit/traj/'
#descriptor = 'Wenatcheev2'
descriptor = 'Yakima'
dates  = pd.date_range('2020-09-10', '2020-09-15 00:00', freq='6H')
#date = datetime(2020,9,13,00)
altitudes   = [50, 100, 500]
names=['trajectory','run','year','month','day','hour','minute','seconds','time','latitude','longitude','altitude','pressure','potential temperature','air temperature','rainfall','mix depth','relative humidity','terrain above msl','solar flux']
#for date in dates:
dstr= dates.strftime('%Y%m%d%H')
       #print(dstr)
numalts  = len(altitudes)
traj_50={}
traj_100={}
traj_500={}
for d in dstr:
   df = pd.read_csv(path+descriptor+d+'.trj',skiprows=6+numalts,delimiter=r"\s+",names= names)
   df['year'] = '20' + df['year'].astype(str)
   df['Datetime']= pd.to_datetime(df.rename(columns={'day':'Day'})[['Day','month','year','hour']])
   #df.set_index(['Datetime'], inplace=True)
   df['trajectory_id']=d
   df=df[:147] #Taaking only 48h of the trajectories
   df_50=df[df['trajectory']==1]
   df_100=df[df['trajectory']==2]
   df_500=df[df['trajectory']==3]
   #df=df[(df['date']>datetime.date(2016,1,1))
   #trj=mpd.TrajectoryCollection(traj_yak,d,'altitude',)
   traj_50[d]=df_50
   traj_100[d]=df_100
   traj_500[d]=df_500

# %%
# Define the location of the map center
map_center = [46.600, -120.500] ## The initial location of the trajectoy.

# Create a map object
m = folium.Map(location=map_center, zoom_start=6)

time_interval = timedelta(hours=1)  # Change to your desired time interval
colors = ['red', 'blue', 'green', 'purple', 'orange']
# Calculate the minimum and maximum start dates across all trajectories
start_dates = sorted([df['Datetime'].iloc[0] for df in traj_50.values()])
min_start_date = min(start_dates)
max_start_date = max(start_dates)
listofd=list(start_dates)

# Create a colormap based on the start dates of each trajectory
colormap = LinearColormap(colors=['blue','green', 'yellow', 'red'], vmin=min_start_date.timestamp(), vmax=max_start_date.timestamp())
#colormap = StepColormap(colors, vmin=min(start_dates), vmax=max(start_dates), caption='Arrival Date')
# Create a GeoJSON FeatureCollection with one feature for each trajectory
# Plot each data frame on the map
for key, df in traj_50.items():

    # Create a feature group for the data frame
    fg = folium.FeatureGroup(name=key)

    # Create a line for the data frame
    start_date = df['Datetime'].iloc[0]
    end_date = df['Datetime'].iloc[-1] 
    delta_date = end_date - start_date
    #dates = list(df['Datetime'])
    #start_date = min(dates)
    #end_date = max(dates)
    #delta_date = end_date - start_date
    #line_color = colors[i % len(colors)]
    line_color = colormap(start_date.timestamp())
    locations = df[['latitude', 'longitude']].values.tolist()
    line = folium.PolyLine(locations, color=line_color, weight=1, opacity=0.7)
    line.add_to(fg)

    # Add markers to the line
    for j in range(int(delta_date / time_interval) + 1):
        marker_time = start_date + j * time_interval
        marker_index = (marker_time - start_date).total_seconds() / (time_interval.total_seconds())
        if j == 0:
            continue
        point_a = locations[int(marker_index - 1)]
        point_b = locations[int(marker_index)]
        line = folium.PolyLine([point_a, point_b], color=line_color, weight=1, opacity=0.7)
        line.add_to(fg)
        marker_lat, marker_lon = point_b
        popup = f"{key} ({marker_lat}, {marker_lon})"
        marker = folium.CircleMarker(location=[marker_lat, marker_lon], radius=0.4, color=line_color, fill=True, fill_color=line_color, fill_opacity=0.5, popup=popup)
        marker.add_to(fg)

            # Add the feature group to the map
    fg.add_to(m)
# # Create a GeoJson object for each trajectory
# geojsons = []
# for name, df in traj_50.items():
#     features = []
#     for _, row in df.iterrows():
#         feature = {
#             'type': 'Feature',
#             'geometry': {
#                 'type': 'Point',
#                 'coordinates': [row['longitude'], row['latitude']]
#             },
#             'properties': {
#                 'time': row['Datetime'].strftime('%Y-%m-%dT%H:%M:%S'),
#                 'icon': 'circle',
#                 'iconstyle': {
#                     'color': color_map(row['Datetime'].iloc[0]),
#                     'fillOpacity': 0.8,
#                     'radius': 5
#                 }
#             }
#         }
#         features.append(feature)

#     # Add the GeoJson object to the list
#     geojson = TimestampedGeoJson({'type': 'FeatureCollection', 'features': features}, name=name)
#     geojsons.append(geojson)

# # Create a MarkerCluster object to display the trajectory points
# marker_cluster = MarkerCluster().add_to(m)

# # Add the GeoJson objects to the map
# for geojson in geojsons:
#     geojson.add_to(marker_cluster)

# # Add the color map legend to the map
# color_map.add_to(m)

# # Position the legend in the upper left corner of the map
# m.get_root().html.add_child(color_map.caption)
# color_map.caption.get_root().set_style('position: fixed; top: 10px; left: 10px; z-index: 9999;')

# # Display the map
# m

# Add the colorbar legend to the map
date_format = DateFormatter('%Y-%m-%d')
colormap.caption = 'Arrival Date'
colormap.add_to(m)
m.add_child(colormap)
#m.get_child(0).add_child(date_format)

# Add a layer control to the map
folium.LayerControl().add_to(m)
# Display the map
m

# %%
m.save("yak_traj_50_2020.png")
# %%
