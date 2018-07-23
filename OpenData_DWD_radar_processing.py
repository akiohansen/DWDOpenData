#!/usr/bin/env python
# -*- coding: utf-8 -*-
# OpenData DWD FX Radar processing tool
# Author: Akio Hansen
# v0: First test version to read in FX binary data and create first static plot
# v1: Included automatic download of latest image and further improvements
# v2: Splittet program into single functions for better code structure
# --> not yet finished, has to be done

# Load necessary python modules
import numpy as np
from datetime import datetime
from sys import argv

# Define user inputs
temp_folder = '/work/bm0834/u233139/DWDRadar/temp/'               # Static link for debug options
#infile_fx  = '/work/bm0834/u233139/DWDRadar/FXLATEST_000_MF002'  # Static link for debug options
#infile_fx   = temp_folder+"FXLATEST_000_MF002"                    # Latest image update option
infile_fx  = sys.argv[1]                                         # Option for argument input
plot_path  = sys.argv[2]                                         # Option for argument input
#plot_path   = '/work/bm0834/u233139/DWDRadar/plots/'              # Static link for debug options


# Define DWD FX Radar download function
def download_dwd_radar(latest):
    "This function downloads latest DWD radar file."

    ## Download latest radar image to temporary directory
    # Check if temporary folder exists otherwise create one
    import os
    import urllib
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    # Download file
    if latest == 1:
        latest_fx_file = urllib.URLopener()
        latest_fx_file.retrieve("https://opendata.dwd.de/weather/radar/composit/fx/FXLATEST_000_MF002", temp_folder+"FXLATEST_000_MF002")
    else:
        print "Error, no other function implemented yet."
        return 1

    # If everything runs without errors, return 0
    return 0


# Define DWD FX Radar processing function
def process_dwd_radar(infile_fx):
    "This function process DWD radar file."

    ### Read and process DWD FX data ###
    ## Documentation about DWD OpenData FX format is found at:
    ## https://www.dwd.de/DE/leistungen/radarprodukte/formatbeschreibung_fxdaten.pdf

    ### Read header information and binary data
    with open(infile_fx, 'r') as f:
        f.seek(0, 0)                                # Set reading cursor to start of file
        file_data   = f.read()                      # Read all lines and store to object
        h_line      = file_data.split('\x03', 1)[0] # Read header information till first ETX
        fx_bin_data = file_data.split('\x03', 1)[1] # Binary data is stored after header ETX

    # Process header information and save variables
    prod_type = h_line[0:2]           # Product type (FX)
    meas_time = h_line[2:8]           # Measurement time UTC (DDhhmm)
    meas_date = h_line[13:17]         # Measurement date (MMYY)
    byte_len  = np.int(h_line[19:26]) # Product length in bytes
    intv_min  = np.int(h_line[51:55]) # Intervall length in minutes
    dim_x     = np.int(h_line[57:61]) # Number of Pixels (x-direction)
    dim_y     = np.int(h_line[62:66]) # Number of Pixels (y-direction)

    # Convert header time information to python time object
    datetime_meas = datetime.strptime(meas_date+'-'+meas_time, '%m%y-%d%H%M')

    # Read binary data from buffer variable and clean buffer
    fx_data = np.frombuffer(fx_bin_data, dtype='<i2').copy()
    del file_data

    # Set missing values according to file format specs.
    fx_data[fx_data > 4095] = 0
    # Data handling for optimized postprocessing
    fx_data = np.asfarray(fx_data)          # Convert data to float
    fx_data = fx_data.reshape(dim_x, dim_y) # Reshape vector to array
    # Convert pixel to dBz values
    fx_data = fx_data/10.0
    fx_data = (fx_data/2.0) - 32.5
    # Set all values below to NaN, non physical values
    fx_data[fx_data < 0.0] = np.nan

    # If everything runs without errors, return dictonary with variables
    return {'dim_x': dim_x, 'dim_y': dim_y, 'datetime_meas': datetime_meas, 'fx_data': fx_data}


# Define DWD FX Radar radar grid function
def create_radar_grid(dim_x, dim_y):
    "This function downloads latest DWD radar file."

    ### Geolocation of pixels
    # Creating array for positions
    grid_pos    = np.empty((dim_x, dim_y,))
    grid_pos[:] = np.nan

    # Reference for Coordinates Calculation see:
    # https://www.dwd.de/DE/leistungen/radolan/radolan_info/radolan_radvor_op_komposit_format_pdf.pdf?__blob=publicationFile&v=7
    #  page 13 and following

    # Set constants
    earthRadius = 6370.04 # km

    junctionNorth = 60.0 # N
    junctionEast  = 10.0 # E

    # Equidistant cartesian coordinates
    dx_dist = 1.0 # km
    dy_dist = 1.0 # km

    # Coordinates of corners of cartesian grid
    lowerleft_x = -523.4622
    lowerleft_y = -4658.645

    lowerright_x = 376.5378
    lowerright_y = -4658.645

    upperright_x = 376.5378
    upperright_y = -3758.645

    upperleft_x = -523.4622
    upperleft_y = -3758.645

    # Create x and y equidistant vectors
    x_cord_cart = np.arange(lowerleft_x,lowerright_x,dx_dist)
    y_cord_cart = np.arange(lowerleft_y,upperleft_y,dy_dist)

    # Convert cartesian coordinates to lat lon coordinates
    x_vec_cord = np.rad2deg( np.arctan(-x_cord_cart[0:900]/y_cord_cart[0]) + np.deg2rad( junctionEast ) ) # see p.13, formula 1.4a

    y_vec_cord = np.rad2deg( np.arcsin(( \
        ((earthRadius**2) * ((1+np.sin(np.deg2rad(junctionNorth)))**2) - ( (x_cord_cart[0:900]**2) + (y_cord_cart[0:900]**2))) / \
        ((earthRadius**2) * ((1+np.sin(np.deg2rad(junctionNorth)))**2) + ( (x_cord_cart[0:900]**2) + (y_cord_cart[0:900]**2))) )) )

    # Create meshgrid to have coordinates matrix
    xm_ll, ym_ll = np.meshgrid(x_vec_cord, y_vec_cord)

    # If everything runs without errors, return dictonary with grid variables
    return {'xm_ll': xm_ll, 'ym_ll': ym_ll}


# Define DWD FX Radar plot function function
def plot_dwd_radar(xm_ll, ym_ll, fx_data, datetime_meas, plot_path, temp_folder):
    "This function gets DWD radar and grid information and creates a plot."

    # Import map data for plotting
    ### Plot DWD Open Radardata on a map
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    # Create formatted time string
    time_string = datetime.strftime(datetime_meas, '%d/%m/%Y %H:%M UTC')
    time_file   = datetime.strftime(datetime_meas, '%d%m%Y_%H%M')

    # Create new figure
    fig = plt.figure(num=None, facecolor='w', edgecolor='k')

    plt.contourf(xm_ll, ym_ll, fx_data, interpolation='none')
    plt.title("Rain at: "+time_string, fontweight="bold", size=18)

    map = Basemap(llcrnrlon=3.5,llcrnrlat=46.5,urcrnrlon=15.5,urcrnrlat=55.5,
                 resolution='i')
    map.drawcoastlines()
    map.drawcountries()
    lat_ticks = [55.0,53.0,51.0,49.0,47.0]
    map.drawparallels(lat_ticks, labels=[1,0,0,0], linewidth=0.0)
    lon_ticks = [4.0,6.0,8.0,10.0,12.0,14.0]
    map.drawmeridians(lon_ticks, labels=[0,0,0,1], linewidth=0.0)

    # Insert colorbar
    cbar = plt.colorbar()
    cbar.set_label('Radar reflectivity (dBz)')

    # Save plot to file
    plt.savefig(plot_path + 'DWD_OpenData_FX_5min_precip_'+ time_file + '.png', dpi=100, format='png')
    plt.close()
    #plt.show()

    # Clean temporary directory
    import os
    os.remove(temp_folder+"FXLATEST_000_MF002")

    # If everything runs without errors, return 0
    return 0


# Main control function
def main():
    # Define parameters
    print "Starting DWD OpenData Radar FX plot generation..."

    # Download latest Radar file
    radar_download = download_dwd_radar(1) # 1 = latest radar image

    # Process Radar file
    radar_process  = process_dwd_radar(infile_fx)

    # Generate DWD RADOLAN grid
    #print(radar_process["dim_x"], radar_process["dim_y"]) # Debug option for grid size
    radar_grid = create_radar_grid(radar_process["dim_x"], radar_process["dim_y"])

    # Generate latest DWD FX Radar plot
    radar_plot = plot_dwd_radar(radar_grid["xm_ll"], radar_grid["ym_ll"], radar_process["fx_data"], \
                    radar_process["datetime_meas"], plot_path, temp_folder)

    # End of main

# Call main routine
main()
print "Finished processing DWD OpenData for radar."
