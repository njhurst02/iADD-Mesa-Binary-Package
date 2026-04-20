#!/usr/bin/env python
# coding: utf-8

# In[38]:

##########################################################################################################   
    
#Module for reading out data from MESA binary simulations, intended for use in tandem with MESA_Binary from iadd.astro.illinois.edu.
#However, the methods used here could easily be modified for single star evolution. 
#For questions, comments, and concerns about this code, please contact njhurst2@illinois.edu

#IMPORTANT NOTES:
#To streamline the process of aligning timesteps, this code assumes that history_interval and profile_interval are
#equal for star 1, star 2, and the binary data.  

#For the kippenhahn plot, you must edit profile_columns.list to include gradr (radiative gradient) and grada (adiabatic gradient). 
#This is done by simply removing the ! before their entries.

##For the abundance plot to work, you must edit profiles_columns.list to include the elements/isotopes you want. 

##########################################################################################################   
    
#Modules

import os
import itertools
import numpy as np
import mesaPlot as mp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

from PyAstronomy import pyasl
from tulips import get_isotopes
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from moviepy.editor import ImageSequenceClip
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Rectangle, Patch
from mesa_reader import MesaProfileIndex, MesaData, MesaLogDir

##########################################################################################################   

#Helper functions

def get_Binary_Data(dir, PM):

    if not PM:
        
        try: #Read in data, with error handling for file reading issues.
            h1 = MesaData(f"{dir}/binary_logs1/history.data") #Primary or donor star
            h2 = MesaData(f"{dir}/binary_logs2/history.data") #Secondary or companion star
            b = MesaData(f"{dir}/binary_history.data") #Binary data
        except (OSError, IOError, FileNotFoundError) as e: #Handle file reading errors, such as missing files or incorrect paths.
            print(f"Error: Could not read the data files. Please check your paths.")
            print(f"Details: {e}")
            return None
        except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
            print(f"An unexpected error occurred: {e}")
            return None
        
        #Star 1 details
        M1 = h1.star_mass
        R1 = 10.**h1.log_R
        T1 = 10.**h1.log_Teff
        Lum1 = 10.**h1.log_L
            
        #Star 2 details
        M2 = h2.star_mass
        R2 = 10.**h2.log_R
        T2 = 10.**h2.log_Teff
        Lum2 = 10.**h2.log_L

        #Binary details
        b_Age = b.age
        b_a = b.binary_separation
        b_period = b.period_days

        #Interpolate Masses, Age, Period, and Separation. They have different lengths than the Radii, Luminoisites, and Temperatures due to how the binary data is recorded.
        x_old = np.linspace(0, 1, len(b_Age))
        x_new = np.linspace(0, 1, len(R1))

        interp_Age = interp1d(x_old, b_Age, kind='cubic')
        interp_a = interp1d(x_old, b_a, kind='cubic')
        interp_period = interp1d(x_old, b_period, kind='cubic')

        Age = interp_Age(x_new)
        a = interp_a(x_new)
        period = interp_period(x_new)
        Total_M = M1 + M2 

        return Age, a, M1, M2, Total_M, R1, R2, T1, T2, Lum1, Lum2, period
    
    if PM:
        
        try: #Read in data, with error handling for file reading issues.
            h1 = MesaData(f"{dir}/binary_logs1/history.data") #Primary or donor star
            b = MesaData(f"{dir}/binary_history.data") #Binary data
        except (OSError, IOError, FileNotFoundError) as e: #Handle file reading errors, such as missing files or incorrect paths.
            print(f"Error: Could not read the data files. Check your paths.")
            print(f"Details: {e}")
            return None
        except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
            print(f"An unexpected error occurred: {e}")
            return None
        
        #Star 1 details
        b_M1 = b.star_1_mass
        R1 = 10.**h1.log_R
        T1 = 10.**h1.log_Teff
        Lum1 = 10.**h1.log_L
            
        #Star 2 details
        b_M2 = b.star_2_mass

        #Binary details
        b_Age = b.age
        b_a = b.binary_separation
        b_period = b.period_days

        #Interpolate Masses, Age, Period, and Separation. They have different lengths than the Radii, Luminoisites, and Temperatures due to how the binary data is recorded.
        x_old = np.linspace(0, 1, len(b_Age))
        x_new = np.linspace(0, 1, len(R1))

        interp_Age = interp1d(x_old, b_Age, kind='cubic')
        interp_a = interp1d(x_old, b_a, kind='cubic')
        interp_M1 = interp1d(x_old, b_M1, kind='cubic')
        interp_M2 = interp1d(x_old, b_M2, kind='cubic')
        interp_period = interp1d(x_old, b_period, kind='cubic')

        Age = interp_Age(x_new)
        a = interp_a(x_new)
        M1 = interp_M1(x_new)
        M2 = interp_M2(x_new)
        period = interp_period(x_new)
        Total_M = M1 + M2 

        return Age, a, M1, M2, Total_M, R1, T1, Lum1, period
    
    else:
        print('PM must be given as a boolean; either True or False.')
        return None

#Functions for Roche lobe geometry plot

#Mass ratio function
def Get_Mass_Ratio(m1, m2, i):
    q = m2[i]/m1[i]
    return q
    
def Get_Lagrange_Points(q):
    # Positions (and potentials) of Lagrange points
    l1, l1pot = pyasl.get_lagrange_1(q)
    l2, l2pot = pyasl.get_lagrange_2(q)
    l3, l3pot = pyasl.get_lagrange_3(q)
    l4, l5 = pyasl.get_lagrange_4(), pyasl.get_lagrange_5()
    return l1, l2, l3, l4, l5
    
def Get_Lagrange_Potential(q):
    # Positions (and potentials) of Lagrange points
    l1, l1pot = pyasl.get_lagrange_1(q)
    l2, l2pot = pyasl.get_lagrange_2(q)
    l3, l3pot = pyasl.get_lagrange_3(q)
    l4, l5 = pyasl.get_lagrange_4(), pyasl.get_lagrange_5()
    
    # Potential for l4 and l5 aren't included in the methods above. Get them separately.
    l4pot = pyasl.rochepot_dl(l4[0], l4[1], l4[2], q)
    l5pot = pyasl.rochepot_dl(l5[0], l5[1], l5[2], q)
    return l1pot, l2pot, l3pot, l4pot, l5pot
    
#Getting potential meshgrid
def Create_Potential_Grid(l2, l3, l4, l5, q):
    x, y = np.linspace(-2.5, 2.5, 300), np.linspace(-2.5, 2.5, 300)
    xx, yy = np.meshgrid(x, y)
    z = 0
    
    # Get dimensional values of Roche potential
    p = pyasl.rochepot_dl(xx, yy, z, q)
    
    #Extent vector just the span of the space
    ext_vect = [-2.5, 2.5, -2.5, 2.5]
    
    return p, ext_vect

#Functions used across multiple plots

#Generate a label for the age of the star(s)
def Label_Star_Age(age_y):
    label = ''
    if age_y >= 1e6 and age_y < 1e9: #Millions of years
        label = 'Star Age: ' + str(round(age_y / 1e6, 3)) + ' Myr'
    elif age_y >= 1e9: #Billions of years
        label = 'Star Age: ' + str(round(age_y / 1e9, 3)) + ' Gyr'
    else: #Years
        label = 'Star Age: ' + str(round(age_y, 3)) + ' yr'
    return label

#Temperatures and corresponding color function, values taken from wikipedia
def Color(T, i): 
    if T[i] > 33000:
        return "#94b4fc"  # Blue
    elif 10000 < T[i] <= 33000:
        return "#a4c4fc"  # Deep bluish white
    elif 7300 < T[i] <= 10000:
        return "#d4e4fc"  # Bluish white
    elif 6000 < T[i] <= 7300:
        return "#fbf4fc"  # White
    elif 5300 < T[i] <= 6000:
        return "#fcece4"  # Yellowish white
    elif 3900 < T[i] <= 5300:
        return "#fcdcb4"  # Pale yellowish orange
    elif 2300 < T[i] <= 3900:
        return "#fcb46c"  # Orangish red
    else:
        return "#fc946c"  # Deep red

#Generate a label for the orbital period
def Label_Period(period_y, i):
    period_label = ''
    if period_y[i] < 1e0:
        label = str(round((period_y[i] / 24) , 3)) + ' hours'
    elif 1e0 < period_y[i] < 365:
        label = str(round(period_y[i], 3)) + ' days'
    else:
        label = str(round((period_y[i] / 365), 3)) + 'years'
    return label

#Function for checking if the inputted fps is valid
def fps_check(fps):
    if (type(fps) != int):
        print('fps must be an integer.')
        return False
    elif (type(fps) == int and fps <= 0):
        print('fps must be a positive integer.')
        return False
    return True

def PM_check(PM):
    if not isinstance(PM, bool):
        print('PM must be given as a Boolean; True or False.')
        return False
    return True

##########################################################################################################

#User specified data out reader

#If the data specified is in a history or binary_history file, it will be returned as a simple numpy array. 
#If the data specified is in a LOGS folder, users will have to specify if they want data from one profile (stored as a numpy array) or all profiles (stored as a list of lists).

def get_Data(dir, parameter):

    profiles_parameters = ['zone', 'mass', 'logR', 'logT', 'logRho', 'logP', 'x_mass_fraction_H', 'y_mass_fraction_He', 'z_mass_fraction_metals', 'log_g', 'logE',
                           'logS', 'logPgas', 'eta', 'grada', 'csound', 'eps_nuc', 'pp', 'cno', 'tri_alpha', 'h1', 'he3', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24', 
                           'si28', 's32', 'ar36', 'ca40', 'fe56', 'opacity', 'luminosity', 'gradT', 'gradr', 'tau']
    history_parameters = ['model_number', 'star_age', 'star_mass', 'log_abs_mdot', 'log_Lnuc', 'log_Lneu', 'log_Lneu_nuc', 'he_core_mass', 'co_core_mass', 'one_core_mass',
                          'fe_core_mass', 'envelope_mass', 'dynamic_timescale', 'kh_timescale', 'nuc_timescale', 'log_Teff', 'log_L', 'log_R', 'log_g','log_center_T',
                            'log_center_Rho', ]
    binary_history_parameters = ['age', 'period_days', 'binary_separation', 'star_1_mass','star_2_mass', 'lg_mtransfer_rate', 'lg_mstar_dot_1', 'lg_mstar_dot_2', 
                                 'lg_system_mdot_1', 'lg_system_mdot_2', 'lg_wind_mdot_1', 'lg_wind_mdot_2', 'J_orb', 'J_spin_1', 'Jdot', 'jdot_mb', 'jdot_gr', 'jdot_ml']

    if (type(parameter) != str):
        print('When using get_data, parameter must given as strings.')
        return None
    
    try:
        if parameter in profiles_parameters:
            print("Which star? (1/2)")
            which_star = input()
            if which_star == '1' or which_star == '2':
                print("One profile or all? (One/All)")
                profile_n_input = input()
                if profile_n_input == "One":
                    print("Which profile? Make sure to input an existing profile index.")
                    which_profile = input()
                    profile = MesaData(f"{dir}/binary_logs{which_star}/profile{which_profile}.data")
                    data = profile.data(parameter)
                elif profile_n_input == "All":
                    data = []
                    profiles_index = MesaProfileIndex(f"{dir}binary_logs{which_star}/profiles.index")
                    profiles = profiles_index.profile_numbers
                    for i in profiles:
                        profile = MesaData(f"{dir}/binary_logs{which_star}/profile{i}.data")
                        data.append(profile.data(parameter))
                else:
                    print("Please type 'One' or 'All'.")
                    return None
            else:
                print("Please type '1' or '2'.")
                return None

        if parameter in history_parameters:
            print("Which star? (1/2)")
            history_input = input()
            if history_input == '1' or history_input == '2':
                h = MesaData(f"{dir}/binary_logs{history_input}/history.data")
                data = h.data(parameter)
        
        if parameter in binary_history_parameters:
            h = MesaData(f"{dir}/binary_history.data")
            data = h.data(parameter)

    except (OSError, IOError, FileNotFoundError) as e: #Handle file reading errors, such as missing files or incorrect paths.
        print(f"Error: Could not read out data. Please check your input parameter, star number, or profile selection.")
        print(f"Details: {e}")
        return 
    except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
        print(f"An unexpected error occurred: {e}")
        return 

    return data

##########################################################################################################

#Used for the EOS diagram and temporal locator

#logT and logRho are enabled by default as part of MESA Binary
def get_profile_temp(dir, logs, index):
    
    index += 1
    profile = MesaData(f"{dir}/{logs}/profile{index}.data")
    temp = profile.data('logT')
    
    return temp

def get_profile_density(dir, logs, index):
    
    index += 1
    profile = MesaData(f"{dir}/{logs}/profile{index}.data")
    density = profile.data('logRho')
    
    return density

def get_profile_mass(dir, logs, index):

    index += 1
    profile = MesaData(f"{dir}/{logs}/profile{index}.data")
    mass = profile.header('star_mass')
    
    return mass

def get_profile_age(dir, logs, index):

    index += 1
    profile = MesaData(f"{dir}/{logs}/profile{index}.data")
    age = profile.header('star_age')
    
    return age

##########################################################################################################

#Temporal locator

#Given a time, locates the index of the temporally closest profile. 
#Can either be called to return a MesaData object of that profile or just the index. 
#Must be given as a time in years (yr), millions of years (Myr), or billions of years (Gyr).
#For example, '10 yr' or '30 Myr'.

#Finds the nearest value in an array given an array and a value
def nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def reverse_Label_Star_age(input_age):
    input_age = input_age.strip()
    if "Gyr" in input_age:
        val = float(input_age.replace( "Gyr", ""))
        val = val * 10 ** 9 
    elif "Myr" in input_age:
        val = float(input_age.replace(" Myr", ""))
        val = val * 10 ** 6
    elif "yr" in input_age:
        val = float(input_age.replace(" yr", ""))
    else:
        raise ValueError(f"Could not convert '{input_age}'. Check format (yr, Myr, Gyr).")
    return int(val)

def get_Time(dir, time):
    print("Star 1 or 2 (1/2)?")
    which_star = input()
    if which_star == '1' or which_star == '2':
        profiles_index = MesaProfileIndex(file_name = f"{dir}/binary_logs{which_star}/profiles.index")
        profiles = profiles_index.profile_numbers
        ages = []
        for i in range(len(profiles)):
            ages.append(get_profile_age(dir, f"binary_logs{which_star}", i))
        time = reverse_Label_Star_age(time)
        nearest_time = nearest(ages, time)
        index = ages.index(nearest_time) + 1
        profile = MesaData(f"{dir}/binary_logs{which_star}/profile{index}.data")
        return nearest_time, index, profile
    else:
        print("Invalid star selection; please input 1 or 2.")
        return None

##########################################################################################################

#Stellar Mass vs. Time

#For examples, see documentation.
#fps is the frames/s of the generated movie, ~12 generally provides good temporal resolution. Must be a positive integer.
def plot_Mass_Transfer(dir, fps, PM):

    if not PM_check(PM):
        return

    if not fps_check(fps):
        return

    data = get_Binary_Data(dir, PM)
    if data is None:
        return 
    
    if not PM:
        Age, a_, M1, M2, Total_M, _, _, T1, T2, _, _, _ = data

        #Calculate colors for animation, not done in animation loop to make plotting previous points easier.
        colors1 = [Color(T1, i) for i in range(len(T1))]
        colors2 = [Color(T2, i) for i in range(len(T2))]   
    
    if PM:
        Age, _, M1, M2, Total_M, _, T1, _, _ = data

        colors1 = [Color(T1, i) for i in range(len(T1))]

    MT_xmin = np.min(Age) - (np.min(Age)/10)
    MT_xmax = np.max(Age) + (np.max(Age)/10)
    MT_ymax = (max(np.max(M1), np.max(M2), np.max(Total_M)))+1  
    
    #The amount of frames should correspond to the number of data points
    #Since all data arrays have the same length, which data are used is irrelevant
    n_frames = len(M1)
        
    os.makedirs(f"{dir}_Mass_Transfer_frames", exist_ok=True)
    frames = []

    if (type(fps) != int):
        print('fps must be given as an integer')
        return
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'sans-serif'

    #Generate a frame for each plotted point
    for i in range(n_frames):
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(1, 1, figure=fig)
        fig.patch.set_facecolor('black')

        ax1 = plt.subplot(gs[0:4, 0:5])
        ax1.set_facecolor('black')
        if not PM:
            plot_MT(ax1, Age[:i+1], M1[:i+1], M2[:i+1], Total_M[:i+1], colors1[:i+1], colors2[:i+1], MT_xmin, MT_xmax, MT_ymax, PM)
        if PM:
            plot_MT(ax1, Age[:i+1], M1[:i+1], M2[:i+1], Total_M[:i+1], colors1[:i+1], 'lawngreen', MT_xmin, MT_xmax, MT_ymax, PM)

        filename = f"{dir}_Mass_Transfer_frames/frame_{i:03d}.png"
        plt.savefig(filename, dpi=100)

        frames.append(filename)

    #Stich all frames together at specified fps to generate movie
    clip = ImageSequenceClip(frames, fps=fps)
    clip.write_videofile(f"{dir}_Mass_Transfer.mp4", codec="libx264")

#Actual plotting function.
def plot_MT(ax, Age, M1, M2, Total_M, colors1, colors2, MT_xmin, MT_xmax, MT_ymax, PM):
    fontsize = 19
    ax.text(0.5, 1.01, 'Stellar Mass vs. Time', transform=ax.transAxes, ha='center', color='white', fontsize=fontsize)

    #X axis formatting.
    x = Age
    ax.tick_params(axis='x', colors='white', labelsize=fontsize)
    ax.set_xlabel('Time (Years)', color='white', labelpad=0, fontsize=fontsize)
    ax.set_xscale('log')
    ax.set_xlim(MT_xmin, MT_xmax)
            
    #Y axis formatting.
    ax.tick_params(axis='y', colors='white', labelsize=fontsize)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_ylabel(r'Mass ($M_{\odot}$)', color='white', labelpad=0, fontsize=fontsize)
    ax.set_ylim(0, MT_ymax)

    #Plot formatting
    for spine in ax.spines.values():
        spine.set_color('white')
        ax.tick_params(which='major', length=7, color='white', direction='out', width=1)
        ax.tick_params(which='minor', length=4, color='white', direction='out', width=0.75)

    # Mass scatter plots.
    if len(M1) == 1:
        #Acounting for the first point.
        ax.scatter(x, M1, c=colors1, s=500, label='Donor',
            edgecolors='none', marker='*')
        ax.scatter(x, M2, c=colors2, s=500, label='Companion',
            edgecolors='none', marker='*')
        ax.scatter(x, Total_M, c='violet', s=100, label='Total Mass',
            edgecolors='none')
    else:
        #All previous points with static trail colors.
        ax.scatter(x[:-1], M1[:-1], c='deepskyblue', s=100, label='_nolegend_',
            edgecolors='none')
        ax.scatter(x[:-1], M2[:-1], c='lawngreen', s=100, label='_nolegend_',
            edgecolors='none')
        ax.scatter(x[:-1], Total_M[:-1], c='violet', s=100, label='_nolegend_',
            edgecolors='none')

        if not PM:    
                # Current (last) points with evolving colors.
            ax.scatter(x[-1], M1[-1], c=colors1[-1], s=500, label='Donor',
                edgecolors='none', marker='*')
            ax.scatter(x[-1], M2[-1], c=colors2[-1], s=500, label='Companion',
                edgecolors='none', marker='*')
            ax.scatter(x[-1], Total_M[-1], c='violet', s=100, label='Total Mass',
                edgecolors='none')
        else:
            # Current (last) points with evolving colors.
            ax.scatter(x[-1], M1[-1], c=colors1[-1], s=500, label='Donor',
                edgecolors='none', marker='*')
            ax.scatter(x[-1], M2[-1], c='lawngreen', s=500, label='Companion',
                edgecolors='none', marker='*')
            ax.scatter(x[-1], Total_M[-1], c='violet', s=100, label='Total Mass',
                edgecolors='none')

    # Remove previous legend
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    #To ensure the marker color in the legend is the color of the most recent point
    current_color1 = colors1[-1]
    current_color2 = colors2[-1]

    #To remove edges on legend markers
    legend_handles = [
    Line2D([], [], marker='o', linestyle='None',
        markerfacecolor='deepskyblue',
        markeredgewidth=0,
        markeredgecolor='none',
        markersize=fontsize, label='Donor'),
    Line2D([], [], marker='o', linestyle='None',
        markerfacecolor='lawngreen',
        markeredgewidth=0,
        markeredgecolor='none',
        markersize=fontsize, label='Companion'),
    Line2D([], [], marker='o', linestyle='None',
        markerfacecolor='violet',
        markeredgewidth=0,
        markeredgecolor='none',
        markersize=fontsize, label='Total Mass')
    ]

    ax.legend(handles=legend_handles, loc='lower left', facecolor='black', labelcolor='white', fontsize=fontsize)

    return

##########################################################################################################

#Hertzsprung-Russel Diagram

#Used for making the observers HR diagram
def colorFader(color1, color2, mix): #Taken from https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
    color1 = np.array(mpl.colors.to_rgb(color1))        
    color2 = np.array(mpl.colors.to_rgb(color2))
    return mpl.colors.to_hex((1-mix)*color1 + mix*color2)

#For examples, see documentation.
#fps is the frames/s of the generated movie, ~12 generally provides good temporal resolution. Must be a positive integer.
def plot_Hertzsprung_Russel(dir, fps, observers, PM):

    fontsize=20

    if not PM_check(PM):
        return

    if not fps_check(fps):
        return

    data = get_Binary_Data(dir, PM)
    if data is None:
        return 
    
    if not PM:

        _, _, _, _, _, _, _, T1, T2, Lum1, Lum2, _ = data

        #Calculate colors for animation, not done in animation loop to make plotting previous points easier.
        colors1 = [Color(T1, i) for i in range(len(T1))]
        colors2 = [Color(T2, i) for i in range(len(T2))]

    if PM:
        _, _, _, _, _, _, T1, Lum1, _ = data
        colors1 = [Color(T1, i) for i in range(len(T1))]

    #The amount of frames should correspond to the number of data points
    #Since all arrays have the same length, which data are used is irrelevant
    n_frames = len(T1)

    #Frame generation
    os.makedirs(f"{dir}_Hertzsprung_Russel_frames", exist_ok=True)
    frames = []

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'sans-serif'

    for i in range(n_frames):
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(1, 1, figure=fig)
        fig.patch.set_facecolor('black')


        ax1 = plt.subplot(gs[0:1, 0:1])
        ax1.set_facecolor('black')
        if not PM:
            plot_HR(ax1, T1[:i+1], Lum1[:i+1], T2[:i+1], Lum2[:i+1], colors1[:i+1], colors2[:i+1], observers, PM)

            #To ensure the marker color in the legend is the color of the most recent point
            current_color1 = colors1[-1]
            current_color2 = colors2[-1]

            #To remove edges on legend markers
            legend_handles = [
            Line2D([], [], marker='o', linestyle='None',
                markerfacecolor='lawngreen',
                markeredgewidth=0,
                markeredgecolor='none',
                markersize=20, label='Companion'),
            Line2D([], [], marker='o', linestyle='None',
                markerfacecolor='deepskyblue',
                markeredgewidth=0,
                markeredgecolor='none',
                markersize=20, label='Donor')
            ]

        if PM:
            plot_HR(ax1, T1[:i+1], Lum1[:i+1], None, None, colors1[:i+1], None, observers, PM)

            #To remove edges on legend markers
            legend_handles = [
            Line2D([], [], marker='o', linestyle='None',
                markerfacecolor='deepskyblue',
                markeredgewidth=0,
                markeredgecolor='none',
                markersize=20, label='Donor')
            ]

            #Remove previous legend
        if ax1.get_legend() is not None:
            ax1.get_legend().remove()

        if observers:
            ax1.legend(
                    handles=legend_handles,
                facecolor='black',
                labelcolor='white',
                bbox_to_anchor=(0.365, 0.275), 
                borderaxespad=0,
                bbox_transform=fig.transFigure,
                fontsize=20
            )
        else:
            ax1.legend(
                handles=legend_handles,
                facecolor='black',
                labelcolor='white',
                loc = 'lower left',
                fontsize=20
            )

        filename = f"{dir}_Hertzsprung_Russel_frames/frame_{i:03d}.png"
        plt.savefig(filename, dpi=100)
        plt.close(fig)
        frames.append(filename)

    #Create animation
    clip = ImageSequenceClip(frames, fps=fps)
    clip.write_videofile(f"{dir}_Hertzsprung_Russel.mp4", codec="libx264")

#Actual Plotting Function
def plot_HR(ax, T1, Lum1, T2, Lum2, colors1, colors2, observers, PM):
    
    fontsize=20
    ax.clear() #Clear plot before formatting to avoid the animation laying over previous plots.
    
    #X axis formatting
    ax.set_xlabel('Surface Temperature (K)', color='white', labelpad=0, fontsize=fontsize)
    ax.set_xlim(2000, 45000)
    ax.set_xscale('log')
    ax.xaxis.label.set_color('white')
    ax.invert_xaxis()
    
    #X axis ticks
    ax.set_xticks([40000, 20000, 10000, 5000, 2500])
    ax.set_xticklabels(['40,000', '20,000', '10,000', '5,000', '2,500'])
    ax.tick_params(axis='x', labelsize=fontsize)
    
    #Y axis formatting
    ax.set_ylabel(fr'Luminosity ($L_{{\odot}}$)', color='white', labelpad=0, fontsize=fontsize)
    ax.set_ylim(10**(-5), 10**6)
    ax.set_yscale('log')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='y', labelsize=fontsize)
    
    #Y axis ticks
    ax.set_yticks([10**i for i in range(-5, 7)])
    ax.set_yticklabels(['$10^{-5}$', '$10^{-4}$', '$10^{-3}$','$10^{-2}$', '$10^{-1}$','$1$', '$10$', '$10^2$', '$10^{3}$', '$10^4$', '$10^{5}$', '$10^6$'])
    
    #Plot formatting
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_color('white')
    ax.text(0.5, 1.01, 'HR Diagram', transform=ax.transAxes, ha='center', color='white', fontsize=fontsize)      

    if observers:

        color_list = ["#94b4fc", "#94b4fc", "#a4c4fc", "#d4e4fc", "#fbf4fc", "#fcece4", "#fcdcb4", "#fcb46c", "#fc946c"]
        range_list = np.array([45000, 33000, 10000, 7300, 6000, 5300, 3900, 2300, 2000])
        log_range_list = np.log10(range_list)
        class_list = ['O', 'B', 'A', 'F', 'G', 'K', 'M']           

        #Set new axis parameters for observers properties & making room for the colorbar
        ax.set_ylim(2*10**(-6), 10**6)
        ax.set_xticks([30000, 10000, 6000, 3000])
        ax.set_xticklabels(['30,000', '10,000', '6,000', '3,000'])
        ax.tick_params(axis='x', labelsize=fontsize)

        #Plotting the class dividers
        for i in range(len(range_list) - 1):
            ax.vlines(range_list[i], 2*10**(-6), 1.5*10**(-5), color='white')

        #Plotting the colorbar at the bottom
        for i in range(len(range_list) - 1):
            n = range_list[i] - range_list[i + 1]
            if n > 10000:
                step = 100
            elif 10000 > n  > 1000:
                step = 10
            else:
                step = 3
            for x in range(0, n + 1, step):
                ax.vlines(x + range_list[i + 1], ymin=10**(-6), ymax=10**(-5), color=colorFader(color_list[i + 1], color_list[i], (x/n)), linewidth=step)

        #Plotting the classes
        for i in range(len(class_list)):
            ax.text((10**((log_range_list[i] + log_range_list[i + 1])/2)), 10**(-5.5), class_list[i], fontsize=fontsize, ha='center', color='black')

    if not PM:
        # Plot donor and companion
        if len(T1) == 1:  # Only one data point
            ax.scatter(T1, Lum1, c=colors1[-1], s=500,
                label='Donor', edgecolors='none', marker='*')
            ax.scatter(T2, Lum2, c=colors2[-1], s=500,
                label='Companion', edgecolors='none', marker='*')
        else:
            # All previous points with static trail colors
            ax.scatter(T1[:-1], Lum1[:-1], c='deepskyblue', s=100, label='_nolegend_', edgecolors='none')
            ax.scatter(T2[:-1], Lum2[:-1], c='lawngreen', s=100, label='_nolegend_', edgecolors='none')
            
            # Current (last) points with evolving colors
            ax.scatter(T1[-1], Lum1[-1], c=colors1[-1], s=500,
                    label='Donor', edgecolors='none', marker='*')
            ax.scatter(T2[-1], Lum2[-1], c=colors2[-1], s=500,
                    label='Companion', edgecolors='none', marker='*')
            
    if PM:
            # Plot donor and companion
        if len(T1) == 1:  # Only one data point
            ax.scatter(T1, Lum1, c=colors1[-1], s=500,
                label='Donor', edgecolors='none', marker='*')
        else:
            # All previous points with static trail colors
            ax.scatter(T1[:-1], Lum1[:-1], c='deepskyblue', s=100, label='_nolegend_', edgecolors='none')
            
            # Current (last) points with evolving colors
            ax.scatter(T1[-1], Lum1[-1], c=colors1[-1], s=500,
                    label='Donor', edgecolors='none', marker='*')

    return

##########################################################################################################

#Roche lobe geometry

#For examples, see documentation.
#fps is the frames/s of the generated movie, ~12 generally provides good temporal resolution. Must be a positive integer.
def plot_Roche_Lobe(dir, fps, PM):

    if not PM_check(PM):
        return

    if not fps_check(fps):
        return

    data = get_Binary_Data(dir, PM)
    if data is None:
        return 

    if not PM:

        Age, a, M1, M2, _, R1, R2, T1, T2, _, _, _ = data

        #Calculate colors for animation, not done in animation loop to make plotting previous points easier.
        colors1 = [Color(T1, i) for i in range(len(T1))]
        colors2 = [Color(T2, i) for i in range(len(T2))]

    if PM:
        Age, a, M1, M2, _, R1, T1, _, _ = data
        colors1 = [Color(T1, i) for i in range(len(T1))]

    #The amount of frames should correspond to the number of data points
    #Since all arrays have the same length, which data are used is irrelevant
    n_frames = len(R1)
    
    #Frame generation
    os.makedirs(f"{dir}_Roche_Lobe_frames", exist_ok=True)
    frames = []

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'sans-serif'

    for i in range(n_frames):
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(1, 1, figure=fig)
        fig.patch.set_facecolor('black')
        
        q = Get_Mass_Ratio(M1, M2, i)

        if q <= 1.:
            q = q
            flip_grid = False
            lpoints = list(Get_Lagrange_Points(q))
            lpot = Get_Lagrange_Potential(q)
            p, ext_vect = Create_Potential_Grid(lpoints[1], lpoints[2], lpoints[3], lpoints[4], q)
        else:
            q = 1.0 / q
            flip_grid = True

            lpoints = list(Get_Lagrange_Points(q))
            lpot = Get_Lagrange_Potential(q)
            p, ext_vect = Create_Potential_Grid(lpoints[1], lpoints[2], lpoints[3], lpoints[4], q)

        if flip_grid:
            p = np.fliplr(p)
            ext_vect = [ext_vect[0] + 0.5, ext_vect[1] + 0.5, ext_vect[2], ext_vect[3]]
            lpoints[0] = -lpoints[0] + 0.5
            lpoints[1] = -lpoints[1] + 0.5
            lpoints[2] = -lpoints[2] + 0.5
            lpoints[3] = (-lpoints[3][0] + 0.5, lpoints[3][1])
            lpoints[4] = (-lpoints[4][0] + 0.5, lpoints[4][1])
        else:
            ext_vect = [ext_vect[0] - 0.5, ext_vect[1] - 0.5, ext_vect[2], ext_vect[3]]
            lpoints[0] = lpoints[0] - 0.5
            lpoints[1] = lpoints[1] - 0.5
            lpoints[2] = lpoints[2] - 0.5
            lpoints[3] = (lpoints[3][0] - 0.5, lpoints[3][1])
            lpoints[4] = (lpoints[4][0] - 0.5, lpoints[4][1])

        ax1 = plt.subplot(gs[0:1, 0:1])
        ax1.set_facecolor('black')
        if not PM:
            plot_RL(ax1, lpoints, lpot, p, R1[i], R2[i], a[i], Age[i], ext_vect, q, colors1[i], colors2[i])
        
        if PM:
            plot_RL(ax1, lpoints, lpot, p, R1[i], 0, a[i], Age[i], ext_vect, q, colors1[i], 'white')
            ax1.add_line(mlines.Line2D([0.475, 0.525], [0.05, -0.05], color='white', linewidth=1)) #Add an x to indicate the position of the companion star in PM cases.
            ax1.add_line(mlines.Line2D([0.525, 0.475], [0.05, -0.05], color='white', linewidth=1)) 
        
        filename = f"{dir}_Roche_Lobe_frames/frame_{i:03d}.png"
        plt.savefig(filename, dpi=100)
        plt.close(fig)
        frames.append(filename)

    #Create animation
    clip = ImageSequenceClip(frames, fps=fps)
    clip.write_videofile(f"{dir}_Roche_Lobe.mp4", codec="libx264")
    
    return

#Actual plotting function
def plot_RL(ax, lpoints, lpot, p, R1, R2, a, Age, ext_vect, q, colors1, colors2):
    l1, l2, l3, l4, l5 = lpoints
    
    fontsize = 20
    ax.clear() #Clear plot before formatting to avoid the animation laying over previous plots.
    
    #X axis formatting
    ax.set_xlabel('x/a', color='white', labelpad=0, fontsize=fontsize)
    ax.tick_params(axis='x', colors='white', labelsize=fontsize)
    ax.set_xlim(-2.5, 2.5)
    
    #Y axis formatting
    ax.set_ylabel('y/a', color='white', labelpad=0, fontsize=fontsize)
    ax.tick_params(axis='y', colors='white', labelsize=fontsize)
    ax.set_ylim(-2.5, 2.5)
    
    #Plot formatting
    for spine in ax.spines.values():
        spine.set_color('white')
    ax.tick_params(colors='white')
    ax.set_aspect('equal', 'box')
    ax.title.set_color('white')
    
    #Plotting contours and Lagrange points
    ax.contour(p, [lpot[3]*1.02, lpot[2], lpot[1], lpot[0]], cmap='Reds', extent=ext_vect)
    ax.text(l1, 0, 'L1', ha='center', color='white', fontsize=fontsize-7.5)
    ax.text(l2, 0, 'L2', ha='center', color='white', fontsize=fontsize-7.5)
    ax.text(l3, 0, 'L3', ha='center', color='white', fontsize=fontsize-7.5)
    ax.text(l4[0], l4[1], 'L4', ha='center', color='white', fontsize=fontsize-7.5)
    ax.text(l5[0], l5[1], 'L5', ha='center', color='white', fontsize=fontsize-7.5)
    ax.text(0.5, 1.01, 'Roche Lobe Geometry', transform=ax.transAxes, ha='center', color='white', fontsize=fontsize+5)

    r1eff = pyasl.roche_lobe_radius_eggleton(q, 1)
    r2eff = pyasl.roche_lobe_radius_eggleton(q, 2)

    #Star geometry
    ax.add_artist(plt.Circle((-0.5, 0), min(R1 / a, r1eff), color=colors1))
    ax.add_artist(plt.Circle((0.5, 0), min(R2 / a, r2eff), color=colors2))
    
    #Simulation age label
    ax.text(0.015, 0.98, Label_Star_Age(Age), transform=ax.transAxes,
               fontsize=fontsize, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), color='white')
               
    return

##########################################################################################################

#Abundances plot

#For this code to work, you must edit profiles_columns.list to include the elements you want.

#Auxillary function for removing elements below a threshold value. 0.001 is normally suitable.
def element_remover(dir, abundances, elements, threshold):
    to_remove = []
    
    for element, profiles in abundances.items():
        flat_profiles = np.concatenate([np.atleast_1d(p) for p in profiles])
        
        if np.max(flat_profiles) < threshold:
            to_remove.append(element)

    for element in to_remove:
        if element in elements:
            elements.remove(element)
        abundances.pop(element, None)
    
    return elements, abundances
    
#Function for getting a unique style for each element displayed, frame by frame.
def get_style(element, element_styles, style_cycle):
    if element not in element_styles:
        element_styles[element] = next(style_cycle)
    return element_styles[element]


##For examples, see documentation.

#fps is the frames/s of the generated movie, ~12 generally provides good temporal resolution. Must be a positive interger.
#Threshold value of 0.001 is normally suitable.
def plot_Abundances(dir, fps, threshold, PM):

    if not PM_check(PM):
        return

    if not fps_check(fps):
        return
    
    if threshold < 0:
        print('Threshold must be a non-negative value')
        return
    elif threshold >= 1:
        print('Threshold must be less than 1')
        return

    #Mesa object and profiles for star 1
    try: #Read in data, with error handling for file reading issues
        M1 = mp.MESA()
        M1.loadHistory(f"{dir}/binary_logs1")
        profiles_index_1 = MesaProfileIndex(file_name=f"{dir}/binary_logs1/profiles.index")

        #For star ages
        m1 = MesaLogDir(log_path= f"{dir}/binary_logs1")
        h1 = MesaData(file_name= f"{dir}/binary_logs1/history.data")
    except (OSError, IOError, FileNotFoundError) as e:
        print(f"Error: Could not read the data files. Check your paths.")
        print(f"Details: {e}")
        return None
    except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
        print(f"An unexpected error occurred: {e}")
        return None
    
    profiles_1 = profiles_index_1.profile_numbers
    models_1 = m1.model_numbers

    #List of elements
    elements_1 = get_isotopes(M1)[:-1]

    ages_1 = []
    for i in models_1:
        ages_1.append(h1.data_at_model_number("star_age", i))
    ages_1 = [Label_Star_Age(age) for age in ages_1]

    # Star 1
    abundances_data_1 = {element: [] for element in elements_1}
    masses_data_1 = []

    for j in profiles_1:
        profile = MesaData(f"{dir}/binary_logs1/profile{j}.data")
        for element in elements_1:
            abundances_data_1[element].append(profile.data(element))
        masses_data_1.append(profile.data('mass'))
        
    #Calculate the mass limit from all profiles
    mass_lim_1 = np.max([masses_data_1[0][0], masses_data_1[-1][0]])

    #Remove elements below the threshold to prevent them from being displayed in the legend
    elements_1, abundances_data_1 = element_remover(dir, abundances_data_1, elements_1, threshold)

    n_frames = len(profiles_1)

    if not PM:

        #Mesa object and profiles for star 2
        try: #Read in data, with error handling for file reading issues
            M2 = mp.MESA()
            M2.loadHistory(f"{dir}/binary_logs2")
            profiles_index_2 = MesaProfileIndex(file_name = f"{dir}/binary_logs2/profiles.index")

            #Star 2 ages
            m2 = MesaLogDir(log_path = f"{dir}/binary_logs2")
            h2 = MesaData(file_name = f"{dir}/binary_logs2/history.data")
        except (OSError, IOError, FileNotFoundError) as e:
            print(f"Error: Could not read the data files. Check your paths.")
            print(f"Details: {e}")
            return None
        except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
            print(f"An unexpected error occurred: {e}")
            return None

        profiles_2 = profiles_index_2.profile_numbers
        models_2 = m2.model_numbers

        #List of elements
        elements_2 = get_isotopes(M2)[:-1]

        ages_2 = []
        for i in models_2:
            ages_2.append(h2.data_at_model_number("star_age", i))
        ages_2 = [Label_Star_Age(age) for age in ages_2]

        # Star 2
        abundances_data_2 = {element: [] for element in elements_2}
        masses_data_2 = []

        for k in profiles_2:
            profile = MesaData(f"{dir}/binary_logs2/profile{k}.data")
            for element in elements_2:
                abundances_data_2[element].append(profile.data(element))
            masses_data_2.append(profile.data('mass'))

        #Calculate the mass limit from all profiles
        mass_lim_2 = np.max([masses_data_2[0][0], masses_data_2[-1][0]])
        
        #Remove elements below the threshold to prevent them from being displayed in the legend
        elements_2, abundances_data_2 = element_remover(dir, abundances_data_2, elements_2, threshold)

        n_frames = np.min([len(profiles_1), len(profiles_2)])
    
    #Colorblind friendly colors
    colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
              '#0072B2', '#D55E00', '#CC79A7', '#000000']
              
    linestyles = ['-', '--']
    style_cycle = itertools.cycle(
        (color, ls) for color in colors for ls in linestyles
    )

    element_styles = {}

    os.makedirs(f"{dir}_Abundances_frames", exist_ok=True)
    frames = []

    fontsize = 20

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'sans-serif'
        
    for i in range(n_frames):

        if PM:
            fig = plt.figure(figsize=(12, 10))
            fig.subplots_adjust(right=0.8) 
            gs = gridspec.GridSpec(1, 1, figure=fig)
            ax1 = plt.subplot(gs[0, 0])
        if not PM:
            fig = plt.figure(figsize=(22, 10))
            gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.25)
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[0, 1])
        
        ax1.clear()
        
        plot_abuns(ax1, abundances_data_1, masses_data_1,
            ages_1, i, mass_lim_1, "Donor Star", threshold,
            element_styles, style_cycle, fontsize, left=True)
        
        handles_dict = {}
        for ax in [ax1]:
            for handle, label in zip(*ax.get_legend_handles_labels()):
                handles_dict[label] = handle
        
        handles = list(handles_dict.values())
        labels = list(handles_dict.keys())

        ax1.legend(
            handles, labels,
            loc="center",
            bbox_to_anchor=(0.9, 0.5), 
            borderaxespad=0,
            bbox_transform=fig.transFigure,
            fontsize=fontsize
        )

        if not PM:
            ax2.clear()
                
            plot_abuns(ax2, abundances_data_2, masses_data_2,
                ages_2, i, mass_lim_2, "Companion Star", threshold,
                element_styles, style_cycle, fontsize, left=False)
    
            handles_dict = {}
            for ax in [ax1, ax2]:
                for handle, label in zip(*ax.get_legend_handles_labels()):
                    handles_dict[label] = handle
            
            handles = list(handles_dict.values())
            labels = list(handles_dict.keys())

            pos1 = ax1.get_position()
            pos2 = ax2.get_position()
            mid_x = (pos1.x1 + pos2.x0) / 2

            ax1.legend(
                    handles, labels,
                    loc="center",
                    bbox_to_anchor=(mid_x, 0.5), 
                    borderaxespad=0,
                    bbox_transform=fig.transFigure,
                    fontsize=fontsize
                )

        filename = f"{dir}_Abundances_frames/frame_{i:03d}.png"
        plt.savefig(filename, dpi=100)
        plt.close(fig)
        frames.append(filename)

    #Create animation
    clip = ImageSequenceClip(frames, fps=fps)
    clip.write_videofile(f"{dir}_Abundances.mp4", codec="libx264")
    
    return
    
#Plotting function
def plot_abuns(ax, abundances, masses, ages, i, mass_lim, which_star, threshold, element_styles, style_cycle, fontsize, left):
    
    for element, profiles in abundances.items():
        if np.max(profiles[i]) > threshold:
            color, ls = get_style(element, element_styles, style_cycle)
            ax.plot(masses[i], profiles[i], label=element,
                    color=color, linestyle=ls, linewidth=5)
            
    ax.set_title(which_star, fontsize=fontsize)

    ax.set_xlabel("Mass [$M_\\odot$]", fontsize=fontsize)
    ax.set_xlim(0, mass_lim)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(top=True, bottom=True,
        axis='x', direction='in',
        which='both', labelsize=fontsize)
    
    ax.set_yscale("log")
    ax.set_ylim(threshold, 2)
        
    ax.text(0.5, 0.97, ages[i],
            transform=ax.transAxes, ha='center', fontsize=fontsize-5)
    ax.text(0.5, 0.96, f"Total Mass: {masses[i][0]:.2f} $M_\\odot$",
            ha="center", va="top", transform=ax.transAxes, fontsize=fontsize-5)
        
    if left == True:
        ax.set_ylabel("Mass Fraction", fontsize=fontsize)
        ax.tick_params(right=False, left=True,
            axis='y', direction='in', which='both', labelsize=fontsize)
               
    else:
        ax.set_ylabel("Mass Fraction", loc='center', fontsize=fontsize)
        ax.yaxis.set_label_position("right")
        ax.tick_params(right=True, left=False,
            axis='y', direction='in', which='both', labelsize=fontsize)
        ax.tick_params(labelleft=False, labelright=True)

##########################################################################################################

#Kippenhahn diagram
#Unlike the rest of the functions, this function does not generate a movie.
#However, because of the amount of data being processed, this generation can take a while.

#For examples, see documentation.

#For a given profile, pull the pp and cno energy generation data in each zone. 
#Say a zone is burning if pp + cno => threshold. 
#Other nuclear reactions are not included but could easily be added. See get_energy_gen for more details. 
#Produces a list of True (meaning a burning zone) or false (meaning not enough energy produced to burn).

def get_profile_and_ages(dir, logs):
    try:
        profiles_index = MesaProfileIndex(file_name=f'{dir}/{logs}/profiles.index')
        profiles = profiles_index.profile_numbers

        #Ages
        m = MesaLogDir(log_path=f'{dir}/{logs}')
        h = MesaData(file_name=f'{dir}/{logs}/history.data')
        models = m.model_numbers

        ages = []
        for i in models:
            ages.append(h.data_at_model_number("star_age", i))

        return ages, profiles
    
    except (OSError, IOError, FileNotFoundError) as e: #Handle file reading errors, such as missing files or incorrect paths.
        print(f"Error: Could not read the data files. Check your paths.")
        print(f"Details: {e}")
        return None
    except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
        print(f"An unexpected error occurred: {e}")
        return None

#Use pp and cno by default
def get_energy_gen(dir, logs, index):

    index += 1
    profile = MesaData(f"{dir}/{logs}/profile{index}.data")
        
    #If you wish to add more reaction networks, follow the methods below. Commented code left as an example. 
    pp = profile.data('pp')
    cno = profile.data('cno')
    #tri_alpha = profile.data('tri_alpha')
    #x = profile.data('x')
    energy = pp + cno
    #energy = pp + cno + tri_alpha + x
    energy_gen = []
    
    for i in range(len(energy)):
        if energy[i] > 1: #Threshold 1 produces reasonable results & fits with convention (not convection!).
            energy_gen.append(True)
        else:
            energy_gen.append(False)
            
    return energy_gen

#For a given profile, pull the mass array for each zone. 
#Produces a list of mass values descending from the stellar surface.
def get_mass(dir, logs, index):

    index += 1
    profile = MesaData(f"{dir}/{logs}/profile{index}.data")
    masses = profile.data('mass')
    
    return masses

#For a given profile, pull adiabatic and radiative gradients in each zone. 
#In zone[i], if gradr > grada, that region is convective. 
#Produces a list of convective or radiative based on the type of zone.
def get_conv_or_rad(dir, logs, index):

    index += 1
    
    profile = MesaData(f"{dir}/{logs}/profile{index}.data")
        
    gradr = profile.data('gradr')
    grada = profile.data('grada')
    conv_or_rad = []
    
    for i in range(len(gradr)):
        if gradr[i] > grada[i]:
            conv_or_rad.append("convective")
        else:
            conv_or_rad.append("radiative")
            
    return conv_or_rad

#For a given profile (taken from get_energy_gen), make a list of every zone where the zone type changes. Starts with False, meaning the first index recorded is where burning stars. List looks like [True, False, True, False, True, ...] but the indexes corresponding to the zones are recorded.
def get_energy_gen_indexes(dir, logs, profiles):

    energy_indexes = []

    for i in range(len(profiles)): #
        energy_gen =  get_energy_gen(dir, logs, i) #Energy generation (T/F) for each zone in a profile
        indexes = [] #Empty list of indexes
        stored = False #Store False, meaning first (and third, fifth, etc) index is True (burning)
        for j in range(1, len(energy_gen)): #Loop over all data in a given profile
            if energy_gen[j] != stored:
                stored = energy_gen[j] #Change the stored value
                indexes.append(j) #Append the index
        energy_indexes.append(indexes) #The full list of lists of indexes
    
    return energy_indexes

#For a given profile (taken from get_conv_or_rad), make a list of every zone where the zone type changes. 
#Starts with radiative, meaning the first index recorded is where convection starts. 
#List looks like [convective, radiative, convective, radiative, convective, ...] but the indexes corresponding to the zones are recorded.
def get_conv_or_rad_indexes(dir, logs, profiles):

    conv_or_rad_indexes = []
        
    # This loop works the same as get_energy_gen_indexes
    for i in range(len(profiles)):
        conv_or_rad = get_conv_or_rad(dir, logs, i)
        indexes = []
        stored = "radiative"
        for j in range(1, len(conv_or_rad)):
            if conv_or_rad[j] != stored:
                stored = conv_or_rad[j]
                indexes.append(j)
            conv_or_rad_indexes.append(indexes)
    
    return conv_or_rad_indexes

#Pull all the energy gen indexes (a list of lists) for a given star. 
#For each list, only pull the mass values corresponding to the zone given by energy_gen_indexes. 
#Produces a list of lists.
def get_mass_indexes_by_energy_gen(dir, logs, profiles):
        
    energy_gen_indexes = get_energy_gen_indexes(dir, logs, profiles) #Get the list of lists based on which star
    energy_indexed_masses = [] #Empty list
    for i in range(len(profiles)):
        profile_indexes = energy_gen_indexes[i] #Pull the energy gen indexes for a given profile
        masses = get_mass(dir, logs, i) #Pull the masses for a given profile
        inner_list = [] #Empty list
        for index in profile_indexes:
            inner_list.append(masses[index]) #To the inner list, append only the mass values corresponding to the energy gen indexes. Remember, these are only the indexes where energy gen changes.
        energy_indexed_masses.append(inner_list) #Append the inner list to the outer list, creating a list of lists

    return energy_indexed_masses

#Pull all the convecive or radiative indexes (a list of lists) for a given star. 
#For each list, only pull the mass values corresponding to the zone given by conv_or_rad_indexes. 
#Produces a list of lists.
def get_mass_indexes_by_conv_or_rad(dir, logs, profiles):
    
    # This loop works the same as get_mass_indexes_by_energy_gen
    conv_or_rad_indexes = get_conv_or_rad_indexes(dir, logs, profiles)
    conv_indexed_masses = []
    for i in range(len(profiles)):
        profile_indexes = conv_or_rad_indexes[i]
        masses = get_mass(dir, logs, i)
        inner_list = []
        for index in profile_indexes:
            # index = index - 1 #Convective or radiative indexes are pulled from the conv_or_rad list, which
            inner_list.append(masses[index])
        conv_indexed_masses.append(inner_list)

    return conv_indexed_masses

#Actual plotting function
def plot_kip(dir, logs, which_star, ax, PM):

    alpha = 0.25

    ages, profiles = get_profile_and_ages(dir, logs)
    
    ages = np.array(ages)
    max_age = max(ages)
    scale_power = int(np.floor(np.log10(max_age)))
    scale_factor = 10**scale_power
    ages = ages/scale_factor

    conv_indexed_masses = get_mass_indexes_by_conv_or_rad(dir, logs, profiles)
    energy_indexed_masses = get_mass_indexes_by_energy_gen(dir, logs, profiles)
    
    #Convective zone ploting function
    for i in range(len(conv_indexed_masses)):
        for j in range(0, len(conv_indexed_masses[i]), 2):
    
            # Edge case where j is the second to or last entry in an inner list
            is_last_entry = (j + 1 >= len(conv_indexed_masses[i]))

            y_pos = conv_indexed_masses[i][j] # Pull the y position (really mass value) for plotting
        
            if is_last_entry:
                # If entry J is the last or second last entry, extend rectangle to bottom.
                # !!! This may need to be revisited.
                height = y_pos
            else:
                # Otherwise, height from convective entry to radiative entry
                height = (conv_indexed_masses[i][j] - conv_indexed_masses[i][j+1])
        
            if i == len(conv_indexed_masses) - 1:
                width = 0.1 # For the last entry i, set the width to 0.1. This is an arbitray value but should not matter with the y (time) scales here.
            else:
                width = ages[i+1] - ages[i] # Width set to cover the distance between the current x coordinate and the next x coordinate.
        
            ax.add_patch(Rectangle(
                (ages[i], y_pos),
                width=width,
                height=-height,
                linewidth=0,
                color='dodgerblue',
                fill=True,
                alpha=alpha + 0.1,
                edgecolor='none'
            ))
    
    #Radiative zone ploting function
    for i in range(len(conv_indexed_masses)):
        for j in range(1, len(conv_indexed_masses[i]), 2):
    
            # Edge case where j is the second to or last entry in an inner list
            is_last_entry = (j + 1 >= len(conv_indexed_masses[i]))

            y_pos = conv_indexed_masses[i][j] # Pull the y position (really mass value) for plotting
        
            if is_last_entry:
                # If entry J is the last or second last entry, extend rectangle to bottom.
                # !!! This may need to be revisited.
                height = y_pos
            else:
                # Otherwise, height from convective entry to radiative entry
                height = (conv_indexed_masses[i][j] - conv_indexed_masses[i][j+1])
        
            if i == len(conv_indexed_masses) - 1:
                width = 0.1 # For the last entry i, set the width to 0.1. This is an arbitray value but should not matter with the y (time) scales here.
            else:
                width = ages[i+1] - ages[i] # Width set to cover the distance between the current x coordinate and the next x coordinate.
        
            ax.add_patch(Rectangle(
                (ages[i], y_pos),
                width=width,
                height=-height,
                linewidth=0,
                color='palegoldenrod',
                fill=True,
                alpha=alpha,
                edgecolor='none'
            ))
    for i in range(len(energy_indexed_masses)):
        for j in range(0, len(energy_indexed_masses[i]), 2):
    
            is_last_entry = (j + 1 >= len(energy_indexed_masses[i]))
        
            y_pos = energy_indexed_masses[i][j]
        
            if is_last_entry:
                height = y_pos
            else:
                height = (energy_indexed_masses[i][j] - energy_indexed_masses[i][j+1])
        
            if i == len(energy_indexed_masses) - 1:
                width = 0.1
            else:
                width = ages[i+1] - ages[i]
            
            ax.add_patch(Rectangle(
                (ages[i], y_pos),
                width=width,
                height=-height,
                linewidth=0,
                color='orangered',
                fill=True,
                alpha=alpha-0.1,
                edgecolor='none'
            ))
            
    # This loop works the same way as above. Plotting third to overlay any burning regions over convective and radiative regions.
    #Burning region ploting function
    for i in range(len(energy_indexed_masses)):
        for j in range(0, len(energy_indexed_masses[i]), 2):
    
            is_last_entry = (j + 1 >= len(energy_indexed_masses[i]))
        
            y_pos = energy_indexed_masses[i][j]
        
            if is_last_entry:
                height = y_pos
            else:
                height = (energy_indexed_masses[i][j] - energy_indexed_masses[i][j+1])
        
            if i == len(energy_indexed_masses) - 1:
                width = 0.1
            else:
                width = ages[i+1] - ages[i]
            
            ax.add_patch(Rectangle(
                (ages[i], y_pos),
                width=width,
                height=-height,
                linewidth=0,
                color='orangered',
                fill=True,
                alpha=alpha,
                edgecolor='none'
            ))
    
    ax.set_title(f"{which_star} Star", fontsize = 16)
    
    fig = ax.get_figure()

    # 2. Define the legend handles
    custom_legend = [
        Patch(facecolor='dodgerblue', edgecolor='black', label='Convective', alpha=alpha),
        Patch(facecolor='palegoldenrod', edgecolor='black', label='Radiative', alpha=alpha),
        Patch(facecolor='orangered', edgecolor='black', label='Burning', alpha=alpha),
        Line2D([0], [0], color='black', label='Star Mass')
    ]
        
    #To get the y axis limit. Probably a more efficient way to do this.
    masses = []
    for i in range(len(ages)):
        masses.append(np.max(get_mass(dir, logs, i)))
    max_mass = np.max(masses)
    
    #Plotting the stellar masses
    ax.plot(ages, masses, color='black', linewidth=1)
            
    ax.set_xlim(ages[0], ages[-1])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel(f"Stellar Age [$10^{{{scale_power}}}$ Years]", fontsize=14)
    ax.ticklabel_format(style='plain', axis='x')
    
    ax.set_ylim(0, max_mass)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel(r"Mass [M$_\odot$]", fontsize=14)
    if which_star == "Companion":
    
        ax.yaxis.set_label_position("right")
        ax.tick_params(left=False, right=True, axis='y', labelright=True, labelleft=False)
        ax.tick_params(axis='y', which='minor', right=True, left=False)

    return custom_legend
 
def plot_Kippenhahn(dir, PM):

    if not PM_check(PM):
        return

    if PM:
        fig = plt.figure(figsize=(12, 10))
        fig.subplots_adjust(right=0.8) 
        gs = gridspec.GridSpec(1, 1, figure=fig)
        ax1 = plt.subplot(gs[0, 0])
        handles = plot_kip(dir, 'binary_logs1', "Donor", ax1, PM)

        ax1.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(1.025, 0.5), 
            borderaxespad=0,
            fontsize=12
        )

    if not PM:
        fig = plt.figure(figsize=(14, 6))
        gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)
        ax1 = plt.subplot(gs[0, 0])
        ax2 = plt.subplot(gs[0, 1]) 

        pos1 = ax1.get_position()
        pos2 = ax2.get_position()
        mid_x = (pos1.x1 + pos2.x0) / 2

        handles = plot_kip(dir, 'binary_logs1', "Donor", ax1, PM)
        plot_kip(dir, 'binary_logs2', "Companion", ax2, PM)

        ax1.legend(
                    handles=handles,
                    loc="center",
                    bbox_to_anchor=(mid_x, 0.5), 
                    borderaxespad=0,
                    bbox_transform=fig.transFigure,
                    fontsize=12
        )

        
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig(f"{dir}_Kippenhahn.png")
    plt.show()

##########################################################################################################

#Equation of state diagram
#For examples, see documentation.

#Many of the functions used in this are found in the Kippenhahn diagram section.

#First, get all:
#Temperature, density, and energy generation for each zone in each profile.
#Plot temperature vs density, with the color of the point corresponding to the energy generation.

#Already have get_energy_gen function from the Kippenhahn section.

def plot_EOS(dir, fps, PM):
    
    if not PM_check(PM):
        return
    
    if not fps_check(fps):
        return
    
        
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'sans-serif'
    fontsize = 20
    os.makedirs(f"{dir}_EOS_frames", exist_ok=True)
    frames = []
    
    if PM:
        try:
            profiles_index_1 = MesaProfileIndex(file_name = f"{dir}/binary_logs1/profiles.index")

            profiles_1 = profiles_index_1.profile_numbers

            fig = plt.figure(figsize=(10, 10))
            gs = gridspec.GridSpec(1, 1, figure=fig)
            ax1 = plt.subplot(gs[0, 0])
            fig.patch.set_facecolor('black')

            for i in range(len(profiles_1)):
                
                energy_1 = get_energy_gen(dir, 'binary_logs1', i) #Energy generation for the first profile, used for color coding in the plot.
                temp_1 = get_profile_temp(dir, 'binary_logs1', i) #Temperature
                density_1 = get_profile_density(dir, 'binary_logs1', i) #Density
                mass_1 = get_profile_mass(dir, 'binary_logs1', i) #Mass
                age_1 = get_profile_age(dir, 'binary_logs1', i) #Age

                ax1.clear()

                plot_equation_of_state(temp_1, density_1, energy_1, age_1, mass_1, "Donor", ax1)

                ax1.legend(
                loc="center",
                bbox_to_anchor=(0.28, 0.7925), 
                borderaxespad=0,
                bbox_transform=fig.transFigure,
                fontsize=fontsize
                )   

                filename = f"{dir}_EOS_frames/frame_{i:03d}.png"
                plt.savefig(filename, dpi=100)
                frames.append(filename) 


        except (OSError, IOError, FileNotFoundError) as e:
            print(f"Error: Could not read the data files. Check your paths.")
            print(f"Details: {e}")
            return None
        except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
            print(f"An unexpected error occurred: {e}")
            return None

    if not PM:
        try:
            profiles_index_1 = MesaProfileIndex(file_name = f"{dir}/binary_logs1/profiles.index")
            profiles_index_2 = MesaProfileIndex(file_name = f"{dir}/binary_logs2/profiles.index")

            profiles_1 = profiles_index_1.profile_numbers
            profiles_2 = profiles_index_2.profile_numbers

            n_frames = np.min([len(profiles_1), len(profiles_2)])

            fig = plt.figure(figsize=(22, 10))
            gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.25)
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[0, 1])
            fig.patch.set_facecolor('black')

            for i in range(n_frames):
                energy_1 = get_energy_gen(dir, 'binary_logs1', i) #Energy generation for the first profile, used for color coding in the plot.
                energy_2 = get_energy_gen(dir, 'binary_logs2', i)
                temp_1 = get_profile_temp(dir, 'binary_logs1', i) #Temperature
                temp_2 = get_profile_temp(dir, 'binary_logs2', i)
                density_1 = get_profile_density(dir, 'binary_logs1', i) #Density
                density_2 = get_profile_density(dir, 'binary_logs2', i)
                mass_1 = get_profile_mass(dir, 'binary_logs1', i) #Mass
                mass_2 = get_profile_mass(dir, 'binary_logs2', i)
                age_1 = get_profile_age(dir, 'binary_logs1', i) #Age
                age_2 = get_profile_age(dir, 'binary_logs2', i)

                ax1.clear()
                ax2.clear()

                plot_equation_of_state(temp_1, density_1, energy_1, age_1, mass_1, "Donor", ax1)
                plot_equation_of_state(temp_2, density_2, energy_2, age_2, mass_2, "Companion", ax2)

                ax1.legend(
                    loc="center",
                    bbox_to_anchor=(0.21, 0.7925), 
                    borderaxespad=0,
                    bbox_transform=fig.transFigure,
                    fontsize=fontsize
                )
            
                filename = f"{dir}_EOS_frames/frame_{i:03d}.png"
                plt.savefig(filename, dpi=100)
                frames.append(filename)

        except (OSError, IOError, FileNotFoundError) as e:
            print(f"Error: Could not read the data files. Check your paths.")
            print(f"Details: {e}")
            return None
        except Exception as e: #Catch-all for any other unexpected errors that may arise during file reading.
            print(f"An unexpected error occurred: {e}")
            return None

    #Create animation
    clip = ImageSequenceClip(frames, fps=fps)
    clip.write_videofile(f"{dir}_EOS.mp4", codec="libx264")

#Plotting function. 
def plot_equation_of_state(temp, density, energy_gen, age, mass, which_star, ax):

    ax.set_facecolor('black')

    #Masks for energy generation criteria
    energy_gen = np.array(energy_gen)
    energy_mask_1 = (1000 >= (energy_gen)) & ((energy_gen) >= 1)
    energy_mask_2 = (10000000 >= (energy_gen)) & ((energy_gen) > 1000)
    energy_mask_3 = ((energy_gen) > 10000000)

    fontsize = 20

    #Plot formatting
    for spine in ax.spines.values():
        spine.set_color('white')
        ax.tick_params(which='major', length=7, color='white', direction='in', width=1, labelsize=fontsize)
        ax.tick_params(which='minor', length=4, color='white', direction='in', width=0.75, labelsize=fontsize)

    ax.set_title(f"{which_star} Star EOS", fontsize=fontsize, color='white')

    ax.set_xlabel('log Density [$g/cm^{3}$]', fontsize=fontsize)
    ax.xaxis.label.set_color('white')
    ax.set_xlim(-10.5, 10.5)
    ax.tick_params(left=True, right=True, axis='y')

    ax.set_ylabel('log Temperature [K]', fontsize=fontsize)
    ax.yaxis.label.set_color('white')
    ax.set_ylim(2.5, 10.5)
    ax.set_xticks([i for i in range(-10, 11)])
    ax.set_xticklabels(['-10', '', '-8', '', '-6', '', '-4', '', '-2', '', '0', '', '2', '', '4', '', '6', '', '8', '', '10'])

    ax.tick_params(colors='white', labelsize=fontsize)
    for spine in ax.spines.values():
        spine.set_color('white')

    ax.plot(density[energy_mask_1], temp[energy_mask_1], linewidth=15, color="#EE7733", label=f"$>$ 1 erg $g^{-1}$ $s^{-1}$")
    ax.plot(density[energy_mask_2], temp[energy_mask_2], linewidth=15, color="#CC3311", label=f"$>$ 1000 erg $g^{-1}$ $s^{-1}$")
    ax.plot(density[energy_mask_3], temp[energy_mask_3], linewidth=15, color="#EE3377", label=f"$>$ 10$^7$ erg $g^{-1}$ $s^{-1}$")
    ax.plot(density, temp, linewidth=10, color="#009988")

    #Plotting nuclear burning thresholds, taken from https://sites.uni.edu/morgans/astro/course/Notes/section2/fusion.html
    #Not perfectly accurate
    h_temp = np.log10(13 * 10 ** 6) #Range from -1 to 4
    h_density = np.log10(100) 
    he_temp = np.log10(100 * 10 ** 6) #Range from 1 to 7
    he_density = np.log10(100000)
    c_temp = np.log10(500 * 10 ** 6) #Range from 2 to 10
    c_density = np.log10(200000)
    o_temp = np.log10(1.5 * 10 ** 9) #Range from 3 to 10
    o_density = np.log10(10 * 10 ** 6)

    h_range = np.linspace(-1, 4, 100)
    h_eq = -1/5 * (h_range - h_density) + h_temp
    he_range = np.linspace(1, 7, 120)
    he_eq = -1/8 * (he_range - he_density) + he_temp
    c_range = np.linspace(2, 10, 160)
    c_eq = -1/8 * (c_range -  c_density) + c_temp
    o_range = np.linspace(3, 10, 140)
    o_eq = -1/8 * (o_range - o_density) + o_temp

    ax.plot(h_range, h_eq, linestyle="dashed", color='white', alpha=0.5)
    ax.plot(he_range, he_eq, linestyle="dashed", color='white', alpha=0.5)
    ax.plot(c_range, c_eq, linestyle="dashed", color='white', alpha=0.5)
    ax.plot(o_range, o_eq, linestyle="dashed", color='white', alpha=0.5)

    ax.text(h_range[0], h_eq[0], 'H burn', ha='center', color='white', fontsize=15)
    ax.text(he_range[0], he_eq[0], 'He burn', ha='center', color='white', fontsize=15)   
    ax.text(c_range[0], c_eq[0], 'C burn', ha='center', color='white', fontsize=15)   
    ax.text(o_range[0], o_eq[0], 'O burn',  ha='center', color='white', fontsize=15)      

    return