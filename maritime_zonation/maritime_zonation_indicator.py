#!/usr/bin/env python
# -*- coding: utf-8

"""
Script Name: maritime_zonation_indicator.py
Description: The following script shows how the list of coordinates (EPSG:4326) could be distributed into maritime zones
Author: Polina Guseva
Date: 05.05.2023

Usage: for everyone, enjoy

Dependencies: pip install -r requirements.txt; python3, see requirements.txt (geopandas==0.12.2, matplotlib==3.7.1, numpy==1.22.4, pandas==1.5.3, shapely==2.0.1)

Notes: Due to the copyright claims we can not provide the reference shapefiles of maritime zones. Please, insert a path to a folder with zips of maritime zones shapefiles

The pipeline of the script:
- upload reference shapefiles of all maritime zones (Flanders Marine Institute (2023): MarineRegions.org. Available online at www.marineregions.org. Consulted on 2023-03-05 (CC BY 4.0));
- upload a table with ur samples where at least longitute and latitude should be stated;
- for each sample get a fit into one of maritime zones (otherwise the exception is raised);
- show nice plots of distribution;
- save a new table with the maritime zone column.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point



## part 1: upload reference shapefiles ##

folder_path: str = input('Please insert a path to a folder with zipped shapefiles of maritime zones: ')
# for example, './maritime_zones_shapefiles'
file_list = os.listdir(folder_path)
print('In this folder are:\n', *file_list)

# unzip each zip per a maritime zone
os.system(f'for file in {folder_path}/*.zip; ' + 'do unzip "$file" -d "${file%.*}"; done')

# combine all shapefiles into one pandas dataframe
# WARNING this script is tailored to a specific set. Please, modify it for yours
df_mar_zones = gpd.GeoDataFrame()

for file_name in file_list:
    if not file_name.endswith('.zip'):

        if file_name in os.listdir(f'{folder_path}/{file_name}/'):
            path_folder_zone = f'{folder_path}/{file_name}/{file_name}/'
        else:
            path_folder_zone = f'{folder_path}/{file_name}/'

        file_shp = list(filter(lambda x: x.endswith('.shp') and 'boundaries' not in x, os.listdir(path_folder_zone)))
        print(file_name, file_shp)
        df_zone = gpd.read_file(path_folder_zone + file_shp[0])

        if file_shp[0].startswith('High_Seas'):
            out_list = [i for i in list(list(df_zone['geometry'])[0].geoms)]
            df_zone = pd.concat([df_zone] * len(out_list), ignore_index=True)
            df_zone.geometry = out_list
            df_zone.area_km2 = df_zone.to_crs(epsg=3857).area / 1e6
            df_zone = df_zone.rename(columns={'name': 'POL_TYPE'})

        df_mar_zones = pd.concat([df_mar_zones, df_zone])

df_mar_zones.POL_TYPE = df_mar_zones.POL_TYPE.apply(lambda x: '200NM' if x == '200 NM' else x)

# let's visualise what we have
df_mar_zones.plot(column='POL_TYPE', cmap='tab20',
                  legend=True).get_legend().set_bbox_to_anchor((1.5, 1))
plt.title('The map of availible maritime zonations')
plt.show()

print('\nTHE FIRST PART IS FINISHED: reference shapefiles uploaded\n')



## part 2: upload a table with ur samples ##
file_path: str = input('Provide a url or a path to an excel table with ur samples: ')
print('\nPlease, recheck the fitting for ur specific data!!!!\n')

df_data = pd.read_excel(file_path, header=2).dropna(axis=1, how='all')

if 'Longitude' not in df_data.columns or 'Latitude' not in df_data.columns:
    raise ValueError('Please recheck that Longitude and Latitude columns are in ur table')

print('\nTHE SECOND PART IS FINISHED: a table with ur samples uploaded\n')



## part 3: indicating maritime zone for each sample ##

# Create a list of points as (longitude, latitude) tuples
points_list = df_data.loc[:, ['Longitude', 'Latitude']].to_numpy()

# Convert the points list into a geopandas dataframe
points_df = gpd.GeoDataFrame(geometry=[Point(lonlat) for lonlat in points_list],
                             crs=df_mar_zones.crs)

# Check if each point falls within any of the polygons using a spatial join and group by
points_in_polygons = gpd.sjoin(points_df, df_mar_zones, op='within')
points_within_polygons = points_in_polygons.groupby('index_right').any().values

# check the presence in any of EEZ
presence_in_zones = [i in points_in_polygons.geometry for i in points_df.geometry]

# Extract the MAR_ZONE value for each point from the joined GeoPandas dataframe
labels = [points_in_polygons.loc[points_in_polygons.index == i, 'MAR_ZONE'].iloc[0] for i in range(len(points_list))]

assert sum(presence_in_zones) / len(presence_in_zones) == 1.0, 'Some samples fall into several zones!'
# if the check is passed we can add labels to the table with ur sample coordinates
points_df['MAR_ZONE'] = labels

print('\nTHE THIRD PART IS FINISHED: a maritime zone is identified per each sample\n')



## part 4: visualisation ##

# calculate the number of samples per each zone
points_dict = points_df.MAR_ZONE.value_counts().to_dict()

# slightly reshape it for visualisation
label_list = [f'{group}: {points_dict[group]}' for group in points_df.MAR_ZONE.unique()]
points_dict2 = dict(zip(points_df.MAR_ZONE.unique(), label_list))
points_df2 = points_df.copy()
points_df2.MAR_ZONE = points_df2.MAR_ZONE.apply(lambda x: points_dict2[x])

fig, ax = plt.subplots(1, 1)
# background with grey maritime zones
df_mar_zones.plot(column='MAR_ZONE', cmap='Greys_r', ax=ax, legend=False)

# actual points
ax2 = ax.twinx()
points_df2.plot(column='MAR_ZONE', cmap='rainbow', legend=True,
                ax=ax2).get_legend().set_bbox_to_anchor((1.5, 0.9))
ax2.get_legend().set_title('Maritime zone\nof geo points:')
ax2.set_yticks([])

ax.set_ylim(-90, 90)
ax2.set_ylim(-90, 90)
ax.set_ylabel('latitude')
ax.set_xlabel('longitude')
plt.title('Geodata for my favourite samples', x=0.5)
plt.show()



## part 5: saving the table with a maritime zonation column ##

df_data['maritime_zone'] = points_df.MAR_ZONE
df_data.to_csv(os.path.splitext(file_path)[0] + '_maritime_zones.csv', index=False)

print('\nTHE PROGRAM IS FINISHED: please, find a new table near an old one with _maritime_zones.csv at th end\n')
