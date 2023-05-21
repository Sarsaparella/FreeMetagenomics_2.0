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
import numba
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd

## part 1: upload reference shapefiles ##

folder_path: str = input('Please insert a path to a folder with zipped shapefiles of maritime zones: ')
# for example, './maritime_zones_shapefiles'
print('In this folder are:\n', *os.listdir(folder_path), sep='\n')

# unzip each zip per a maritime zone
os.system(f'for file in {folder_path}/*.zip; ' + 'do unzip "$file" -d "${file%.*}"; done')

# combine all shapefiles into one pandas dataframe
# WARNING this script is tailored to a specific set. Please, modify it for yours
df_mar_zones = gpd.GeoDataFrame()
file_list = os.listdir(folder_path)

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
plt.title('The map of available maritime zonations')
plt.show()

print('\nTHE FIRST PART IS FINISHED: reference shapefiles uploaded\n')

## part 2: upload a table with ur samples ##

file_path: str = input('Provide a url or a path to an excel table with ur samples: ')
print('\nPlease, recheck the fitting for ur specific data!!!!\n')

df_data = pd.read_excel(file_path, header=2).dropna(axis=1, how='all')

if 'Longitude' not in df_data.columns or 'Latitude' not in df_data.columns:
    raise ValueError('Please recheck that Longitude and Latitude columns are in ur table')

# Create a list of points as (longitude, latitude) tuples
points_list = df_data.loc[:, ['Longitude', 'Latitude']].to_numpy()

print('\nTHE SECOND PART IS FINISHED: a table with ur samples uploaded\n')


## part 3: indicating maritime zone for each sample ##

# some functions to detect whatever a point is inside a polygon or not
@numba.njit
def is_inside_postgis(polygon, point):
    """
    Determines if a point is inside a polygon using the PostGIS algorithm.

    Parameters:
    -----------
    polygon : list
        A list of (longitude, latitude) tuples representing the vertices of the polygon.
    point : tuple
        A (longitude, latitude) tuple representing the point to check.

    Returns:
    --------
    int
        If the point is inside the polygon, returns True.
        If the point is on the boundary of the polygon, returns 2.
        Otherwise, returns False.
    """
    length = len(polygon)
    intersections = 0

    dx2 = point[0] - polygon[0][0]
    dy2 = point[1] - polygon[0][1]
    jj = 1

    while jj < length:
        dx = dx2
        dy = dy2
        dx2 = point[0] - polygon[jj][0]
        dy2 = point[1] - polygon[jj][1]

        F = (dx - dx2) * dy - dx * (dy - dy2)
        if 0.0 == F and dx * dx2 <= 0 and dy * dy2 <= 0:
            return 2

        if (dy >= 0 > dy2) or (dy2 >= 0 > dy):
            if F > 0:
                intersections += 1
            elif F < 0:
                intersections -= 1

        jj += 1

    return intersections != 0


@numba.njit(parallel=True)
def is_inside_postgis_parallel(points, polygon):
    """
    Determines if a list of points is inside a polygon using the PostGIS algorithm,
    using parallel processing.

    Parameters:
    -----------
    points : list
        A list of (longitude, latitude) tuples representing the points to check.
    polygon : list
        A list of (longitude, latitude) tuples representing the vertices of the polygon.

    Returns:
    --------
    numpy.ndarray
        A boolean array of length `len(points)`, where each element is True
        if the corresponding point is inside the polygon,
        and False otherwise.
    """
    ln = len(points)
    D = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        D[i] = is_inside_postgis(polygon, points[i])
    return D


def sample_polygon_relation(zone_poly_row, row_number):
    """
    Given a GeoPandas DataFrame row with a polygon, checks
    if any points in the global list points_list are inside
    the polygon using the function is_inside_postgis_parallel.

    If there are any matches, the function appends
    the corresponding Metagenomes names and the name of the POL_TYPE of the row
    to the global lists list_names and list_polys, respectively.

    Parameters:
    zone_poly_row (GeoPandas row): A GeoPandas DataFrame row containing a polygon.
    row_number (int): The index of the row in the DataFrame.

    Returns:
    None
    """
    global list_names, list_polys, points_list
    global df_data, df_mar_zones

    polygon = np.array(zone_poly_row.exterior.coords)
    list_of_matches = is_inside_postgis_parallel(points_list, polygon)
    sum_of_matches = sum(list_of_matches)
    if sum_of_matches > 0:  # if there is at least one match
        for one_match in df_data[list_of_matches].Metagenomes:
            list_names += [one_match]
            list_polys += [df_mar_zones.POL_TYPE.iloc[row_number]]


# apply each polygon to a list of points and save overlapping events
list_names = []
list_polys = []

for row_number in range(len(df_mar_zones)):  # each polygon
    zone_poly_row = df_mar_zones.geometry.iloc[row_number]

    try:  # in case of a polygon
        sample_polygon_relation(zone_poly_row, row_number)

    except AttributeError:  # in case of a multipolygon
        list_of_polygons = list(zone_poly_row.geoms)  # deconstruct it to polygons
        for polygon in list_of_polygons:
            sample_polygon_relation(polygon, row_number)

# create a new df with samples and maritime zones
df_new = pd.DataFrame([list_names, list_polys], index=['sample', 'mar_zone'])
df_new = df_new.T.groupby('sample').agg({'mar_zone': lambda x: list(x)})
df_new = df_new.reset_index()

# if the order of samples are the same as in the original file
# let's add a column with detected maritime zones!
if all(df_data.Metagenomes == df_new['sample']):
    print('They are equal!')
    df_data['maritime_zones'] = df_new['mar_zone']

print('\nTHE THIRD PART IS FINISHED: a maritime zone is identified per each sample\n')



## part 4: visualisation ##

# to our new dataframe let's add coordinates
df_new['Longitude'] = df_data['Longitude']
df_new['Latitude'] = df_data['Latitude']

# let's replace each sample with geographically controversial position
df_new['mar_zone'] = df_new['mar_zone'].apply(lambda x: x[0] if len(x) == 1 else 'controversial')
# overall, how many categories right now?
count_dict = {element: df_new['mar_zone'].str.count(element).sum() for element in df_new['mar_zone'].unique()}
# prepare a legend
label_list = [f'{group}: {count_dict[group]}' for group in count_dict.keys()]
points_legend = dict(zip(count_dict.keys(), label_list))

fig, ax = plt.subplots(1, 1, figsize=(10, 5))

# draw a background with maritime zones
df_mar_zones.plot(column='POL_TYPE', cmap='Greys_r', ax=ax,
                  legend=False)
ax.set_ylim(-90, 95)
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')

# add dotes as samples on top
ax2 = ax.twinx()
category_to_number = dict(zip(df_new['mar_zone'].unique(),
                              range(len(df_new['mar_zone'].unique()))))
colors = df_new['mar_zone'].map(category_to_number)
ax2.scatter(x=df_new.Longitude, y=df_new.Latitude, c=colors, cmap='rainbow')
handles, _ = ax2.collections[0].legend_elements()
ax2.legend(handles, points_legend.values(), bbox_to_anchor=(1.3, 0.9),
           title='Maritime zone\nof geo points:')
ax2.set_yticks([])
ax2.set_yticklabels([])
ax2.set_ylim(-90, 95)
ax2.grid()

plt.suptitle('Geodata of TARA_Oceans from doi: 10.1038/s41564-018-0176-9', y=0.95)
plt.show()



## part 5: saving the table with a maritime zonation column ##

df_data.to_csv(os.path.splitext(file_path)[0] + '_maritime_zones.csv', index=False)

print('\nTHE PROGRAM IS FINISHED: please, find a new table near an old one with _maritime_zones.csv at th end\n')
