# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a set of functions to simulate an energy
## management system controlled by a genetic algorithm.
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################

import matplotlib.pyplot as plt
import csvreader
# from numpy import savetxt
from math import ceil
from numpy import arange, array, append
from GA_EMS import *
import time


###############################################################################
################################   Simulations ################################


start_day = 0#200 #200 - 32 - 115 - 200
end_day = 1#207  #207 - 39 - 275 - 215
horizon = 0#4*2#4*1
consecutive_generations = 5
individuals = 25
#
labels = ['19/07', '20/07', '21/07', '22/07', '23/07', '24/07', '25/07', '26/07']
#labels = ['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02']

method = False # [False, 'linear_random', 'RandomForests']
list_format = 'multiple' # ['single', 'multiple']

file_name = 'PF_H' + str(int(horizon/4)) + '_G' + str(int(consecutive_generations)) + '_P' + str(int(individuals)) + '.csv'


time_registry = []

# [SC, TESS, HP, PV, BESS, Grid]
energy_costs = [0, 0, 0, 0, 0.13, 0.483]
CO2_costs = [0, 0, 0, 0, 0, 0.325]

optimization_weights = [1, 1, 1]#[0.483*0.25/2, 10000/1, 1/2.5] # theta_E, theta_T, theta_CO2 
                                                  # RT = 0.483*0.25/2, 1/1, 1/2.5
                                                  # H = 1 h = 0.483*0.25/2, 100/1, 1/2.5                                                  
beta = 6

# PV
n_modules = 10
module_power_ref = 0.315
module_power = 0.400   

CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
CSVDataP_Load = csvreader.read_data(csv='Load_Profile_15min.csv', address='', delim=',')
CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
CSVDataPrices = csvreader.read_data(csv='DA_Prices_15min.csv', address='')
CSVDataPV.data2array()
CSVDataTamb.data2array()
CSVDataP_Load.data2array()
CSVDataRad.data2array()
CSVDataPrices.data2array()
P_PV_av = [i[0]*n_modules*module_power/module_power_ref/1000 for i in CSVDataPV.ar]
T_amb = [i[0]+273 for i in CSVDataTamb.ar[start_day*24*4:(end_day + ceil(horizon/(24*4)))*24*4]]
P_Load = [i for i in CSVDataP_Load.ar[0]]
a = arange(0,len(CSVDataRad.ar),15)
G = array([CSVDataRad.ar[i][0] for i in a[start_day*24*4:(end_day + ceil(horizon/(24*4)))*24*4]])
Energy_price = [i[0]*0.25 for i in CSVDataPrices.ar[start_day*24*4:(end_day + ceil(horizon/(24*4)))*24*4]] # [energy_costs[-1]*0.25 for i in CSVDataPrices.ar]



# Initial conditions
t = array([0])                      # In s
dt = 60*15                          # In s
t_final = int((end_day - start_day)*24*3600/dt)         # In s

# TESS
T_TESS = array([75 + 273])          # In K
Qdot_TESS = array([0, 0])              # In W
TESS_active = True

# HP
HP_active = True
HP_status = HP_Power(HP_active)
P_HP = array([HP_status[0],HP_status[0],HP_status[0]])       # In kW
Qdot_HP = array([HP_status[1],HP_status[1],HP_status[1]])    # In kW


# Thermal demand
T_0 = 20 + 273
T_in = array([T_0, T_0, T_0])          # In K
T_set_day = [T_0 - 3]*int((6-0)*4) + [T_0]*int((22-6)*4)+ [T_0 - 3]*int((24-22)*4)
T_set = array(T_set_day*((end_day + ceil(horizon/(24*4)))-start_day))
Qdot_Losses = array([0, 0, 0])
Qdot_SC_TESS = array([0])

# Solar Collectors
SC_active = True
Qdot_SC = array([0, 0, 0])                # In W

# PV

P_PV = array([0]) 

# BESS
SoC_BESS = array([.50])                     # In %
P_BESS = array([0])                         # In kW
SoCmax = 0.9
SoCmin = 0.2
charge_efficiency = 0.943
discharge_efficiency = 0.943


# Grid
P_Grid = array([0])                         # In kW


ts = 0
day = 0
i_step = 0


start = time.time()    # The timer is initializad.
for step in range(t_final-1):

    i_step += 1

    start_step = time.time()    # The timer is initializad.
    h = step + horizon + 1
    
    horizon_costs = [[],[]]
    for i in range(horizon + 1):
        energy_costs[-1] = Energy_price[step + horizon]
        horizon_costs[0].append(energy_costs)        
        horizon_costs[1].append(CO2_costs)
    
    
    population_registry.append([])
    best_candidate_registry.append([])    
    
    best_candidate = GA_Optimization(
            T_amb = TS_forecast(T_amb[step:h], step, method, forecast_type = 'Temperature'),
            T_0_in = T_in[-1],
            T_set = T_set[step:h],
            T_soil = TS_forecast(T_amb[step:h], step, method, forecast_type = 'Temperature'),
            P_PV_av = TS_forecast(P_PV_av[step:h], step, method, forecast_type = 'PV'),
            P_Load = TS_forecast(P_Load[step:h], step, method, forecast_type = 'Load'),
            G = TS_forecast(G[step:h], step, method, forecast_type = 'Radiation'),
            T_0_TESS = T_TESS[-1],
            SoC_0_BESS = SoC_BESS[-1],
            SoC_BESS = SoC_BESS,
            horizon = horizon,
            costs = horizon_costs,
            consecutive_generations = consecutive_generations,
            individuals = individuals,
            beta = beta,
            theta_E = optimization_weights[0],
            theta_T = optimization_weights[1],
            theta_CO2 = optimization_weights[2])

    SC_active = best_candidate[0]
    TESS_active = best_candidate[1]
    HP_active = best_candidate[2]
    PV_curtail = best_candidate[3]
    P_BESS_active = best_candidate[4]

    model_state = Thermal_Electrical_model(T_amb = T_amb[step], T_0_in = T_in[-1], T_set = T_set[step], T_soil = T_amb[step], P_PV_av = P_PV_av[step], P_Load = P_Load[step], G = G[step], SC_active = SC_active, TESS_active = TESS_active, HP_active = HP_active, PV_curtail = PV_curtail, P_BESS = P_BESS_active, T_0_TESS = T_TESS[-1], SoC_0_BESS = SoC_BESS[-1], dt = dt)

    
    Qdot_SC = append(Qdot_SC, model_state[0])
    Qdot_SC_TESS = append(Qdot_SC_TESS, model_state[1])
    Qdot_TESS = append(Qdot_TESS, model_state[2])
    Qdot_HP = append(Qdot_HP, model_state[3])
    Qdot_Losses = append(Qdot_Losses, model_state[4])
    T_TESS = append(T_TESS, model_state[5])
    T_in = append(T_in, model_state[6])
    P_HP = append(P_HP, model_state[7])
    P_PV = append(P_PV, model_state[8])
    P_BESS = append(P_BESS, model_state[9])
    P_Grid = append(P_Grid, model_state[10])
    SoC_BESS = append(SoC_BESS, model_state[11])

    t = append(t, dt*step/3600)


    end_step = time.time()    # The timer is initializad.    
    time_registry.append(end_step - start_step)    

    ts += 1
    if ts == 96:
        ts = 0
        day += 1
        print('Current day: ', start_day + day, '. Current temperature: ', T_in[-1]-273)
    

end = time.time()    # The timer is initializad.
totalelapsed = end - start  # The total time is calculated.


extras = True

if extras:   
    
    ##################################   Plots  ###################################
    
    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })
    

    plt.figure(constrained_layout=True)
    plt.plot(t, [i-273 for i in T_TESS])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Temperature in the TESS, $T_{TESS}$, [°C]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t, [i/1000 for i in Qdot_TESS[1:]])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Thermal power of the TESS, $\dot{Q}_{TESS}$, [kW]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t, [i-273 for i in T_amb], 'r', label='$T_{amb}$')
    plt.plot(t, [i-273 for i in T_in[2:]], 'b', label='$T_{in}$')
    plt.plot(t, [i-273 for i in T_set[start_day*24*4:end_day*24*4]], 'g', label='$T_{set}$')
    plt.legend(loc='lower right')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Temperature [°C]')
    plt.show()    
    
    plt.figure(constrained_layout=True)
    plt.plot(t, [i/1000 for i in Qdot_Losses[2:]])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylim([0, 1.8])
    plt.ylabel('Thermal losses, $\dot{Q}_{L}$, [kW]')
    plt.show()    
    
    plt.figure(constrained_layout=True)
    plt.plot(t, [i/1000 for i in Qdot_HP[2:]])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Thermal power of the HP, $\dot{Q}_{HP}$, [kW]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t, [i/1000 for i in Qdot_SC[2:]], 'b', label='To the thermal load')
    plt.legend(loc='lower center')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Thermal power of the SC, $\dot{Q}_{SC}$, [kW]')
    plt.show()

    plt.figure(constrained_layout=True)
    plt.plot(t, P_Load[start_day*4*24:end_day*4*24])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Electric load, $P_{L}$, [kW]')
    plt.show()  
    
    plt.figure(constrained_layout=True)
    plt.plot(t, P_PV, 'r', label='PV power curtailed')
    plt.plot(t, P_PV_av[start_day*24*4:end_day*4*24], 'b', label='PV power available')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('PV power load, $P_{PV}$, [kW]')
    plt.show()  
    
    plt.figure(constrained_layout=True)
    plt.plot(t, SoC_BESS)
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('State-of-charge of the BESS, $SoC_{BESS}$ [%]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t, P_BESS, 'b', label='Power delivered by the BESS')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Power of the BESS, $P_{BESS}$, [kW]')
    plt.show()
    
    plt.figure(constrained_layout=True)
    plt.plot(t, P_Grid, 'b', label='Power consumed from the grid')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Power from the grid, $P_{G}$, [kW]')
    plt.show()


    plt.figure(constrained_layout=True)
    plt.plot(t, [i/0.25 for i in Energy_price], 'b', label='Power consumed from the grid')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xlabel('Time [time steps]')
    plt.ylabel('Energy price, $c_{E}^{grid}$, [€/kWh]')
    plt.show()    
    
    
    
    
    
    fig18, (axs18) = plt.subplots(1,2)
    axs18[0].boxplot(time_registry)
    axs18[0].set_ylabel('Computation time per timestep [s]')
    axs18[0].set_xticklabels(['With outliers'])
#    axs18[0].set_title('With outliers')
    
    axs18[1].boxplot(time_registry, showfliers=False)
    axs18[1].set_ylabel('Computation time per timestep [s]')
    axs18[1].set_xticklabels(['Without outliers'])    
#    axs18[1].set_title('Without outliers')
    #plt.title('Computation time')
    plt.show()