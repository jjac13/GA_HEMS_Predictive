# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################
import time
from numpy import append, array, cos, pi, arange

###############################################################################
###########################   Functions definition  ###########################

###############################################################################
##########################   House model definition  ##########################

def heat_loss(U, Area, T_in, T_out):
    return U * Area * (T_out - T_in)
    
def new_house_Temperature(T_0, Qdot_Losses, mc_T, Qdot_HP = 0, Qdot_TESS = 0, Qdot_SC = 0, dt=1):
    T_new = T_0 + dt*(Qdot_Losses + Qdot_HP + Qdot_TESS + Qdot_SC)/(mc_T)
    
    return T_new

def update_TESS(active, T_0, T_soil, Qdot_SC = 0, Qdot_HP = 0, Tmax = 95 + 273, Tmin = 50 + 273, mdot = 0.1, T_network = 40 + 273, Qdot_SD = 100, efficiency = 0.8, m = 4000, c = 4200, dt=1*3600):

    if active and T_0 >= Tmin:
        Qdot_TESS = mdot*c*(T_0*0.97 - T_network)
    else:
        Qdot_TESS = 0

    
    if T_0 <= Tmin:         # TESS discharged, only charge
        if T_0 <= T_soil:   # Code for heat dissipation/absoprtion through the soil needed
            T_new = T_0
            
        else:
            T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS)*dt/(m*c)
            
    elif T_0 <= Tmax:       # TESS available for charge and discharge
        
        T_new = T_0 + (Qdot_SC*efficiency + Qdot_HP*efficiency - Qdot_SD - Qdot_TESS)*dt/(m*c)
        
    else:                   # TESS fully charged, only discharge
        
        T_new = T_0 + (- Qdot_SD - Qdot_TESS)*dt/(m*c)        
        
    return [T_new, Qdot_TESS*efficiency]          

def Qdot_SolarCollector(active, G, A = 6, SC_eff = 0.45, dt = 1*3600): # A = 6
    
    if active:
        return 0*A*SC_eff*G/dt
    else:
        return 0

# Enphase IQ3 https://enphase.com/download/iq-battery-3-data-sheet

def update_BESS(SoC_0, P_BESS, P_Load, P_PV = 0, SoCmax = 0.9, SoCmin = 0.2, P_BESS_max = 1.28, P_Grid_max = 0, Capacity_BESS = 3.36, charge_efficiency = 0.943, discharge_efficiency = 0.943, P_SD = 0, dt = 0.25):

    
    E_BESS_0 = Capacity_BESS*SoC_0
    
    if P_BESS > 0:
        E_BESS = E_BESS_0 - P_BESS*dt/discharge_efficiency
    elif P_BESS <= 0:
        E_BESS = E_BESS_0 - P_BESS*dt*charge_efficiency
    
    E_BESS = E_BESS*(1-P_SD)
    SoC_BESS = E_BESS/Capacity_BESS
    
    return SoC_BESS
    

def BESS_perm_min(SoC, Capacity_BESS = 3.36, SoCmax = 0.9, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmax)/dt, -P_BESS_max, P_BESS_max)
    
    
def BESS_perm_max(SoC, Capacity_BESS = 3.36, SoCmin = 0.2, P_BESS_max = 1.28, dt = 0.25):
    from numpy import clip
    
    return clip(Capacity_BESS*(SoC - SoCmin)/dt, -P_BESS_max, P_BESS_max)
    
        
def HP_Power(active, P_in = 2.7*1000, COP = 4.1):
    if active:
        return [P_in/1000, P_in*COP]
    else:
        return [0,0]
    
def House_Thermal_Losses(T_0_in, T_amb):

    #############################   Paremeters ################################
    
    # Convective heat transfer coefficients [W/m^2K]
    h_air_wall   = 0.9#24       # Indoor air -> walls, scaled to R-value of a C-label house
    h_wall_atm   = 0.9#34       # Walls -> atmosphere, scaled to R-value of a C-label house
    h_air_window = 25            # Indoor air -> windows
    h_window_atm = 32            # Windows -> atmosphere
    h_air_roof   = 12            # Indoor air -> roof
    h_roof_atm   = 38            # Roof -> atmosphere
    
    ## House
    
    # Air
    c_air        = 1005.4        # Specific heat of air at 273 K [J/kgK]
    airDensity   = 1.025         # Densiity of air at 293 K [kg/m^3]
    kAir         = 0.0257        # Thermal conductivity of air at 293 K [W/mK]
    
    
    # Windows (glass)
    n1_window     = 3            # Number of windows in room 1
    n2_window     = 2            # Number of windows in room 2
    n3_window     = 2            # Number of windows in room 3
    n4_window     = 1            # Number of windows in room 4
    htWindows     = 1            # Height of windows [m]
    widWindows    = 1            # Width of windows [m]
    windows_area  = (n1_window + n2_window + n3_window + n4_window) * htWindows * widWindows
    LWindow       = 0.004        # Thickness of a single window pane [m]
    LCavity       = 0.014        # Thickness of the cavity between the double glass window [m]  
    windowDensity = 2500         # Density of glass [kg/m^3]
    c_window      = 840          # Specific heat of glass [J/kgK]
    kWindow       = 0.8          # Thermal conductivity of glass [W/mK]
    U_windows = ((1/h_air_window) + (LWindow/kWindow) + (LCavity/kAir) + (LWindow/kWindow) + (1/h_window_atm))**-1
    m_windows = windowDensity * windows_area * LWindow
    
    # Walls (concrete)
    lenHouse    = 15             # House length [m]
    widHouse    = 8              # House width [m]
    htHouse     = 2.6            # House height [m]
    LWall       = 0.25           # Wall thickness [m]
    wallDensity = 2400           # Density [kg/m^3]
    c_wall      = 750            # Specific heat [J/kgK]
    kWall       = 0.14           # Thermal conductivity [W/mK]
    walls_area = 2*(lenHouse + widHouse) * htHouse - windows_area
    U_wall = ((1/h_air_wall) + (LWall /kWall) + (1/h_wall_atm))**-1
    m_walls = wallDensity * walls_area * LWall
    
    # Roof (glass fiber)
    pitRoof     = 40/180/pi      # Roof pitch (40 deg)
    LRoof       = 0.2            # Roof thickness [m]
    roofDensity = 2440           # Density of glass fiber [kg/m^3]
    c_roof      = 835            # Specific heat of glass fiber [J/kgK]
    kRoof       = 0.04           # Thermal conductivity of glass fiber [W/mK]
    roof_Area = 2 * (widHouse/(2*cos(pitRoof))*lenHouse)
    U_roof = ((1/h_air_roof) + (LRoof/kRoof) + (1/h_roof_atm))**-1
    m_roof = roofDensity * roof_Area * LRoof
    
    
    m_air = airDensity * lenHouse * widHouse * htHouse
    
    mc_T = m_air*c_air + m_roof*c_roof + m_windows*c_window + m_walls*c_wall
     
    
    ############################   Calculations ###############################    
    ################### Thermal carrier ######################
    
    ## Thermal losses
    # Roof
    Qdot_roof = heat_loss(U_roof, roof_Area, T_0_in, T_amb)
    
    # Windows
    Qdot_windows = heat_loss(U_windows, windows_area, T_0_in, T_amb)
    
    # Walls
    Qdot_wall = heat_loss(U_wall, walls_area, T_0_in, T_amb)
    
    return [Qdot_roof + Qdot_windows + Qdot_wall, mc_T]
    
def Thermal_Electrical_model(T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G, SC_active, TESS_active, HP_active, PV_curtail, P_BESS, T_0_TESS, SoC_0_BESS, dt = 0.25):

#    #############################   Paremeters ################################
#    
#    # Convective heat transfer coefficients [W/m^2K]
#    h_air_wall   = 0.9#24       # Indoor air -> walls, scaled to R-value of a C-label house
#    h_wall_atm   = 0.9#34       # Walls -> atmosphere, scaled to R-value of a C-label house
#    h_air_window = 25            # Indoor air -> windows
#    h_window_atm = 32            # Windows -> atmosphere
#    h_air_roof   = 12            # Indoor air -> roof
#    h_roof_atm   = 38            # Roof -> atmosphere
#    
#    ## House
#    
#    # Air
#    c_air        = 1005.4        # Specific heat of air at 273 K [J/kgK]
#    airDensity   = 1.025         # Densiity of air at 293 K [kg/m^3]
#    kAir         = 0.0257        # Thermal conductivity of air at 293 K [W/mK]
#    
#    
#    # Windows (glass)
#    n1_window     = 3            # Number of windows in room 1
#    n2_window     = 2            # Number of windows in room 2
#    n3_window     = 2            # Number of windows in room 3
#    n4_window     = 1            # Number of windows in room 4
#    htWindows     = 1            # Height of windows [m]
#    widWindows    = 1            # Width of windows [m]
#    windows_area  = (n1_window + n2_window + n3_window + n4_window) * htWindows * widWindows
#    LWindow       = 0.004        # Thickness of a single window pane [m]
#    LCavity       = 0.014        # Thickness of the cavity between the double glass window [m]  
#    windowDensity = 2500         # Density of glass [kg/m^3]
#    c_window      = 840          # Specific heat of glass [J/kgK]
#    kWindow       = 0.8          # Thermal conductivity of glass [W/mK]
#    U_windows = ((1/h_air_window) + (LWindow/kWindow) + (LCavity/kAir) + (LWindow/kWindow) + (1/h_window_atm))**-1
#    m_windows = windowDensity * windows_area * LWindow
#    
#    # Walls (concrete)
#    lenHouse    = 15             # House length [m]
#    widHouse    = 8              # House width [m]
#    htHouse     = 2.6            # House height [m]
#    LWall       = 0.25           # Wall thickness [m]
#    wallDensity = 2400           # Density [kg/m^3]
#    c_wall      = 750            # Specific heat [J/kgK]
#    kWall       = 0.14           # Thermal conductivity [W/mK]
#    walls_area = 2*(lenHouse + widHouse) * htHouse - windows_area
#    U_wall = ((1/h_air_wall) + (LWall /kWall) + (1/h_wall_atm))**-1
#    m_walls = wallDensity * walls_area * LWall
#    
#    # Roof (glass fiber)
#    pitRoof     = 40/180/pi      # Roof pitch (40 deg)
#    LRoof       = 0.2            # Roof thickness [m]
#    roofDensity = 2440           # Density of glass fiber [kg/m^3]
#    c_roof      = 835            # Specific heat of glass fiber [J/kgK]
#    kRoof       = 0.04           # Thermal conductivity of glass fiber [W/mK]
#    roof_Area = 2 * (widHouse/(2*cos(pitRoof))*lenHouse)
#    U_roof = ((1/h_air_roof) + (LRoof/kRoof) + (1/h_roof_atm))**-1
#    m_roof = roofDensity * roof_Area * LRoof
#    
#    
#    m_air = airDensity * lenHouse * widHouse * htHouse
#    
#    mc_T = m_air*c_air + m_roof*c_roof + m_windows*c_window + m_walls*c_wall
#     
#    
#    ############################   Calculations ###############################    
#    ################### Thermal carrier ######################
#    
#    ## Thermal losses
#    # Roof
#    Qdot_roof = heat_loss(U_roof, roof_Area, T_0_in, T_amb)
#    
#    # Windows
#    Qdot_windows = heat_loss(U_windows, windows_area, T_0_in, T_amb)
#    
#    # Walls
#    Qdot_wall = heat_loss(U_wall, walls_area, T_0_in, T_amb)
# 
#    losses_registry[0].append(Qdot_roof)
#    losses_registry[1].append(Qdot_windows)
#    losses_registry[2].append(Qdot_wall)    
    
    [Qdot_Losses, mc_T] = House_Thermal_Losses(T_0_in, T_amb)
    
#    Qdot_losses = Qdot_roof + Qdot_windows + Qdot_wall


    # Solar Collector
    Qdot_SC = Qdot_SolarCollector(SC_active, G)
    
    # TESS     
    Qdot_SC_TESS = Qdot_SolarCollector(not SC_active, G)
    TESS_state = update_TESS(TESS_active, T_0_TESS, T_soil = T_soil, Qdot_SC = Qdot_SC_TESS, dt = dt)
    T_TESS = TESS_state[0]
    Qdot_TESS =  TESS_state[1]

    
    # HP        
    HP_state = HP_Power(HP_active)
    P_HP = HP_state[0]
    Qdot_HP = HP_state[1]        
    
    # Thermal demand
    T_in = new_house_Temperature(T_0_in, Qdot_Losses, mc_T, Qdot_HP = Qdot_HP, Qdot_TESS = Qdot_TESS, Qdot_SC = Qdot_SC, dt = dt)  
#    Qdot_Losses = -(Qdot_roof + Qdot_windows + Qdot_wall)


    ################### Electric carrier ######################

    # PV
    P_PV = P_PV_av*PV_curtail

    # BESS
#   [SoC_BESS_state, P_Grid_state, P_BESS_state] = update_BESS(SoC_BESS, P_Load + P_HP/1000, P_PV)
    SoC_BESS = update_BESS(SoC_0_BESS, P_BESS, P_Load + P_HP/1000)
    
    
    # Grid balance
    
    P_Grid = P_Load + P_HP - P_PV - P_BESS
    
    return [Qdot_SC, Qdot_SC_TESS, Qdot_TESS, Qdot_HP, -Qdot_Losses, T_TESS, T_in, P_HP, P_PV, P_BESS, P_Grid, SoC_BESS]

###############################################################################
#############################   GA EMS definition  ############################

## The function ---------------------------------------------------------------

def date_time_list(horizon, start_date_time, list_format = 'multiple', dt = 0.25):
    
    start_date_time = start_date_time*dt/24
    
    if list_format == 'single':
        # Format is in days alone. The decimal parts are associated to the hours and minutes.
        
        return [start_date_time+i*(dt/24) for i in range(horizon+1)]
        
    elif list_format == 'multiple':
        from math import floor
        # Format is in [day, hour, minute]. The decimal parts are associated to the hours and minutes.
    
    return [[floor(start_date_time+i*(dt/24)), floor((start_date_time+i*(dt/24))%1*24), floor((start_date_time+i*(dt/24))*24%1*60)] for i in range(horizon+1)]

    
## The function ---------------------------------------------------------------
    
def TS_forecast(TS, start_date_time = 0, method = 'linear_random', forecast_type='PV', list_format = 'multiple', dt = 0.25):
    
    
    if not method:
        return TS
    
    elif method == 'linear_random':
        from numpy import random        
        
        random_uncertanty_limit = 0.05/96
        return [TS[i]*random.uniform(1-i*random_uncertanty_limit, 1+i*random_uncertanty_limit) for i in range(len(TS))]
        
    elif method == 'RandomForests':
        
        import pandas as pd
        from joblib import load
        
        start_time=pd.Timestamp(year=2023, month=1, day=1) + pd.Timedelta(minutes=start_date_time*dt*60)
        timestamps=[start_time + pd.Timedelta(minutes=i*dt*60) for i in range(len(TS))]
        
        # Create the DataFrame
        X_test = pd.DataFrame({
            'Timestamp': timestamps,
            forecast_type+'_lag1': TS
        })

        # Extract 'Hour' and 'Month' from the Timestamp
        X_test['Hour'] = X_test['Timestamp'].dt.hour
        X_test['Month'] = X_test['Timestamp'].dt.month
        # Drop the 'Timestamp' column as it's no longer needed
        X_test = X_test.drop('Timestamp', axis=1)
        rf = load(forecast_type+'_forecast.joblib')
        y_pred = rf.predict(X_test)
        #date_list = date_time_list(len(TS)-1, start_date_time, list_format, dt)
        date_list=y_pred.tolist()
#        print(date_list)
        return date_list

## The function ---------------------------------------------------------------

def chromosome(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5):
# Thermal_Electrical_model(T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G, SC_active, TESS_active, HP_active, PV_curtail, P_BESS, T_0_TESS, SoC_0_BESS, dt = 0.25):    
    from numpy import append, random

# For discrete PV curtailment.    
#    chromosome = array([random.randint(0,2), random.randint(0,2), random.randint(0,2), random.randint(0,2), round(random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1])),1)])
    
# For continous PV curtailment.
    
    chromosome = array([])
    
    for i in range(horizon + 1):
    
        
        # For discrete PV curtailment.    
#        chromosome = append(chromosome, array([random.randint(0,2), random.randint(0,2), random.randint(0,2), random.randint(0,2), round(random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1])),1)]))
    
        # For continous PV curtailment.
        chromosome = append(chromosome, array([random.randint(0,2), random.randint(0,2), random.randint(0,2), round(random.uniform(0, 1),1), round(random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1])),1)]))
    
       
    # After appending the h elements in the horizon
#    population_registry[-1].append(list(chromosome))
    chromosome = append(chromosome, cost_function(chromosome, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2))
    
    return chromosome

## The function ---------------------------------------------------------------
    
def initial_population(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, individuals, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5):
    
    population = array([chromosome(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon) for i in range(individuals)])
    
    from collections import Counter
    l_strin = [str(i) for i in population]
    diversity_registry.append(len(Counter(l_strin).keys()))
    
    return population

## The function ---------------------------------------------------------------

## Considerations:
##  - Consider the size of the population in case the elite_rate and 
##    non_elite_rate variables are modified, as their sum should result into an
##    even number when multiplied by the size of the population.


def create_new_population(population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, elite_rate = 0.05, non_elite_rate = 0.15, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5):    
    from numpy import append, random, array
    
    # Elite selection
#    print(population)
#    l = array([i[-1] for i in population])    
#    print(l.argsort())
#    population = population[l.argsort()]  
#    print(population)
    population = population[population[:,-1].argsort()]
    elite = population[int((1-elite_rate)*len(population)): len(population), :]
    
    # Non-elite selection
    non_elite = population[0:int((1-elite_rate)*len(population)), :]
    non_elite = non_elite[random.choice(non_elite.shape[0], int(len(population)*non_elite_rate), replace=False), :]
    
    
    # Mating
    mating_population = append(elite, non_elite, axis=0)    
    mating_population = mating_population[random.choice(list(range(0,len(mating_population))),len(mating_population), replace=False)]    

#    print(mating_population)    

    new_generation = create_new_generation(mating_population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2)
    
#    print(new_generation)     
    
    # Inclusion of new generation
    
    new_population = append(population, new_generation, axis=0)
    
    # Selecting survivors
    
    new_population = new_population[random.choice(new_population.shape[0], len(population), replace=False), :]
    
    return new_population

## The function ---------------------------------------------------------------
    
## Considerations:
##  - The number of genes to be crossed has to be at least 1, with a maximum of
##    1 less than the number of genes per chromosome.

def create_new_generation(mating_population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, rand_gene_number = False, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5):  # , rand_gene_number = False
    from numpy import append, ceil, random, reshape
    
    new_generation = array([])   
    number_of_genes = len(mating_population[0,:])
    
    
    if rand_gene_number:        
        cross_genes = random.choice(list(range(0,number_of_genes)),random.randint(1,number_of_genes), replace=False)
    else:
#        cross_genes = random.choice(list(range(0,number_of_genes)),int(ceil(number_of_genes)), replace=False)
        cross_genes = random.choice(list(range(0,number_of_genes)),int(ceil(number_of_genes)/2), replace=False)


    for i in range(int(len(mating_population)/2)):         

        
        chromosome_1 = mating_population[i]
        chromosome_2 = mating_population[int(len(mating_population)/2)+i]
        
        new_chromosome_1 = array([i for i in chromosome_1])
        new_chromosome_2 = array([i for i in chromosome_2])

        
        for i in cross_genes:            
            new_chromosome_1[i] = chromosome_2[i]
            new_chromosome_2[i] = chromosome_1[i]
            
        # Here the position for the cost is replaced with the actual cost
        
#        population_registry[-1].append(list(new_chromosome_1))        
#        new_chromosome_1[-1] = cost_function(new_chromosome_1, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2)
#        population_registry[-1].append(list(new_chromosome_1)) 
        
#        print(new_chromosome_1)
#        print(cost_function(new_chromosome_1, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2))
        
        new_chromosome_1[-4:] = cost_function(new_chromosome_1, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2)
        
#        print(new_chromosome_1)
        
#        population_registry[-1].append(list(new_chromosome_2))        
#        new_chromosome_2[-1] = cost_function(new_chromosome_2, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2)
#        population_registry[-1].append(list(new_chromosome_2))  
        new_chromosome_2[-4:] = cost_function(new_chromosome_2, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2)
        
        new_chromosomes = append(new_chromosome_1, new_chromosome_2, axis=0)
        new_generation = append(new_generation, new_chromosomes, axis=0)
    

    
    new_generation = reshape(new_generation, (-1,number_of_genes))
    return new_generation
    
## The function ---------------------------------------------------------------
    

def best_individual(population):    
    from numpy import argmin, array
#    from statistics import mean
#    l = array([i[-1] for i in population])
#    print(population)
#    print(l)
#    
#    print('Min costs', min([i[-1] for i in population]))
#    a = population[argmin(population[:,-1])]
#    
#    if min([i[-1] for i in population]) == a[-1]:
#        print('True')

    return population[argmin(population[:,-1])]
#    return population[argmin(l[:,-1])]    
    
    
## The function ---------------------------------------------------------------

def cost_function(chromosome, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5):# theta_E = 1/26, theta_CO2 = 1/4
#    from numpy import random
    from numpy import clip, sqrt
#    CO2_costs = [0.071, 0.00155, 0.1, 0.041, 0.1916, 0.325]
    
    energy_costs = costs[0]
    CO2_costs = costs[1]

# New code    
    electric_cost = 0
    thermal_cost = 0
    CO2_cost = 0
    cost = 0
# New code    
    
    
    for i in range(horizon + 1):
        model_state = Thermal_Electrical_model(T_amb = T_amb[i], T_0_in = T_0_in, T_set = T_set[i], T_soil = T_amb[i], P_PV_av = P_PV_av[i], P_Load = P_Load[i], G = G[i], SC_active = chromosome[5*i + 0], TESS_active = chromosome[5*i + 1], HP_active = chromosome[5*i + 2], PV_curtail = chromosome[5*i + 3], P_BESS = chromosome[5*i + 4], T_0_TESS = T_0_TESS, SoC_0_BESS = SoC_0_BESS, dt = dt)
        
        if model_state[9] < 0:
        
            Variable_powers = [(model_state[0] + model_state[1])/1000, model_state[2]/1000, model_state[3]/1000, model_state[8], model_state[9], model_state[10]]
        
        else:
            Variable_powers = [(model_state[0] + model_state[1])/1000, model_state[2]/1000, model_state[3]/1000, model_state[8], 0, model_state[10]]
        
# Original code        
#        cost = cost + theta_E*sum([Variable_powers[j]*energy_costs[j] for j in range(len(energy_costs))]) + theta_T*beta*abs(T_set[i] - model_state[6]) + theta_CO2*sum([clip(Variable_powers[j], 0, a_max = None)*CO2_costs[j] for j in range(len(CO2_costs))])
#        cost = cost + theta_E*model_state[10] + theta_T*beta*abs(T_set[i] - model_state[6]) + theta_CO2*sum([Variable_powers[j]*CO2_costs[j] for j in range(len(CO2_costs))])
# Original code    
            
# New code            
        electric_cost += sum([c*P for c, P in zip(energy_costs[i], Variable_powers)])
        thermal_cost += beta*abs(T_set[i] - model_state[6])
        CO2_cost += sum([clip(c*P,0, a_max = None) for c, P in zip(CO2_costs[i], Variable_powers)]) 
# New code
            
        T_0_in = model_state[6]
        T_0_TESS = model_state[5]
        SoC_0_BESS = model_state[11]            
    

    
#    model_state = Thermal_Electrical_model(T_amb = T_amb, T_0_in = T_0_in, T_set = T_set, T_soil = T_amb, P_PV_av = P_PV_av, P_Load = P_Load, G = G, SC_active = chromosome[0], TESS_active = chromosome[1], HP_active = chromosome[2], PV_curtail = chromosome[3], P_BESS = chromosome[4], T_0_TESS = T_0_TESS, SoC_0_BESS = SoC_0_BESS, dt = dt)
#    Qdot_SC = append(Qdot_SC, model_state[0])
#    Qdot_SC_TESS = append(Qdot_SC_TESS, model_state[1])
#    Qdot_TESS = append(Qdot_TESS, model_state[2])
#    Qdot_HP = append(Qdot_HP, model_state[3])
#    Qdot_Losses = append(Qdot_Losses, model_state[4])
#    T_TESS = append(T_TESS, model_state[5])
#    T_in = append(T_in, model_state[6])
#    P_HP = append(P_HP, model_state[7])
#    P_PV = append(P_PV, model_state[8])
#    P_BESS = append(P_BESS, model_state[9])
#    P_Grid = append(P_Grid, model_state[10])
#    SoC_BESS = append(SoC_BESS, model_state[11])
    
#    if model_state[9] < 0:
#    
#        Variable_powers = [(model_state[0] + model_state[1])/1000, model_state[2]/1000, model_state[3]/1000, model_state[8], model_state[9], model_state[10]]
#    
#    else:
#        Variable_powers = [model_state[0]/1000, model_state[2]/1000, model_state[3]/1000, model_state[8], 0, model_state[10]]
        
#    
        
# Original code        
#    electric_cost = sum([Variable_powers[j]*energy_costs[j] for j in range(len(energy_costs))])
#    thermal_cost = beta*abs(T_set[i] - model_state[6])
#    CO2_cost = sum([clip(Variable_powers[j],0, a_max = None)*CO2_costs[j] for j in range(len(CO2_costs))])
# Original code            
    
#    
#    cost = theta_E*sum([Variable_powers[i]*energy_costs[i] for i in range(len(energy_costs))]) + theta_T*beta*abs(T_set - model_state[6]) + theta_CO2*sum([Variable_powers[i]*CO2_costs[i] for i in range(len(CO2_costs))])
#    print('Electric cost: ', electric_cost, ', Thermal cost: ', thermal_cost, ', Total cost: ', cost)
    
    
#    cost = theta_E*electric_cost + theta_T*thermal_cost + theta_CO2*CO2_cost
    cost = sqrt(electric_cost**2 + thermal_cost**2 + CO2_cost**2)
#    population_registry[-1][-1].append(electric_cost)
#    population_registry[-1][-1].append(thermal_cost)
#    population_registry[-1][-1].append(CO2_cost)
#    population_registry[-1][-1].append(cost)    


#    return cost 
    return [electric_cost, thermal_cost, CO2_cost, cost]
    
## The function ---------------------------------------------------------------

def GA_Optimization(T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, horizon, costs = [0, 0, 0, 0, 0.13, 0.483], consecutive_generations = 5, individuals = 200, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5): # costs = [0.15, 1.4, 0.11, 0.18, 0.13, 0.483]

    if population_registry[-1]:
#        print(1)
#        print(population_registry[-1])
        population = population_registry[-1]
    else:
#        print(2)
        population = initial_population(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G, T_0_TESS, SoC_0_BESS, horizon, individuals)
    
#    print(population[-1])    
    
#    population = initial_population(costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G, T_0_TESS, SoC_0_BESS, horizon, individuals)
    best_candidate = best_individual(population)
    
#    population_registry[-1].append(population)
#    best_candidate_registry[-1].append(best_candidate)
    
    generation = 0
    number_of_candidates=1
    iteration = 0

    population_registry.append(population)    
    
    while iteration < consecutive_generations:
#        print('  Generation: ', generation, ' Iteration: ', iteration)
        new_population = create_new_population(population, costs, T_amb, T_0_in, T_set, T_soil, P_PV_av, P_Load, G,  T_0_TESS, SoC_0_BESS, beta = beta, theta_E = theta_E, theta_T = theta_T, theta_CO2 = theta_CO2)
#        population_registry[-1].append(new_population)
        new_candidate = best_individual(new_population)
#        print('  -----New candidate-----')
#        print(' ', new_candidate[-1])        
        
        if new_candidate[-1] < best_candidate[-1]:
            population_registry.pop()
            population_registry.append(new_population)
            
#            best_candidate_registry[-1].append(new_candidate)
            
            best_candidate = [i for i in new_candidate]
            iteration = -1
            number_of_candidates+=1
#            print('  -----New best candidate found-----')
#            print(' ', best_candidate[-1])
##        
#        print('Current generation: ', generation)
#        print('Current iteration: ', iteration)
        generation+=1  
        iteration+=1
        
    
#    print('The best candidate from the initial population is: ', best_individual(population))
#    print('The best candidate from the last generation is: ', best_candidate)   
#    print('The number of generations was: ', generation)
    
    n_generations.append(generation)
    
    return best_candidate    

###############################################################################
##################################   Plotting #################################

def week_graph_TQP(t, T_in, T_amb, T_set, Qdot_TESS, Qdot_SC, Qdot_HP, Qdot_Load, P_Grid, P_Load, P_HP, P_PV, P_BESS, labels, Energy_price, dt = 0.25, days = 7):
    
    import matplotlib.pyplot as plt
    
    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman"
    })
    
    fig, [T_ax1, T_ax, Qdot_ax] = plt.subplots(3, 1, sharex = True)
    #    fig, [T_ax, Qdot_ax] = plt.subplots(2, 1, sharex = True)
#    fig, [T_ax, Qdot_ax, P_ax] = plt.subplots(3, 1)
 
    T_ax1.plot([i-273 for i in T_amb], 'r', label='$T_{amb}$', linewidth=1) # , label='Ambient temperature, $T_{amb}$'
    T_ax1.set_ylabel('Temperature [°C]')
    T_ax1.legend(loc='upper right', ncol=3)
    T_ax1.xaxis.set_tick_params(labelbottom=False)
    T_ax1.set_xticks([i*24/dt for i in range(days +1)]) # 
#    T_ax.set_xticklabels([str(i) for i in range(days +1)], fontsize = 12)
    T_ax1.set_xlim([0,days*24/dt])
    T_ax1.grid()
    
#    T_ax.plot([i-273 for i in T_amb], 'r', label='$T_{amb}$', linewidth=1) # , label='Ambient temperature, $T_{amb}$'
    T_ax.plot([i-273 for i in T_set], 'g', label='$T_{set}$', linewidth=1) # , label='Setpoint temperature inside the house, $T_{set}$'    
    T_ax.plot([i-273 for i in T_in], 'b', label='$T_{in}$', linewidth=1) # , label='Temperature inside the house, $T_{in}$'
    T_ax.set_ylabel('Temperature [°C]')
    T_ax.legend(loc='lower right', ncol=3)
    T_ax.xaxis.set_tick_params(labelbottom=False)
    T_ax.set_xticks([i*24/dt for i in range(days +1)]) # 
#    T_ax.set_xticklabels([str(i) for i in range(days +1)], fontsize = 12)
    T_ax.set_xlim([0,days*24/dt])
    T_ax.grid()
    
    Qdot_ax.plot([-i/1000 for i in Qdot_Load], 'r', label='$\dot{Q}_{L}$', linewidth=1)    
    Qdot_ax.plot([i/1000 for i in Qdot_SC], 'b', label='$\dot{Q}_{SC}$', linewidth=1)
    Qdot_ax.plot([i/1000 for i in Qdot_TESS], 'g', label='$\dot{Q}_{TESS}$', linewidth=1)    
    Qdot_ax.plot([i/1000 for i in Qdot_HP], 'c', label='$\dot{Q}_{HP}$', linewidth=1)
    Qdot_ax.set_ylabel('Thermal power [kW]')
    Qdot_ax.legend(loc='center right', ncol=4)
#    Qdot_ax.xaxis.set_tick_params(labelbottom=False)
    Qdot_ax.set_xticks([i*24/dt for i in range(days +1)]) # 
    Qdot_ax.set_xticklabels(labels, fontsize = 12)
    Qdot_ax.set_xlim([0,days*24/dt])
    Qdot_ax.grid()    
   
#    P_ax.plot([-i for i in P_Load], 'r', label='$P_{L}$', linewidth=1)    
#    P_ax.plot([-i for i in P_HP], 'c', label='$P_{HP}$', linewidth=1)        
#    P_ax.plot([-i for i in P_PV], 'b', label='$P_{PV}$', linewidth=1)
#    P_ax.plot([i for i in P_Grid], 'y', label='$P_{Grid}$', linewidth=1)    
#    P_ax.plot([i for i in P_BESS], 'g', label='$P_{BESS}$', linewidth=1)
#    P_ax.set_ylabel('Thermal power [kW]')
#    P_ax.legend(loc='upper left', ncol=5)
#    P_ax.set_xticks([i*24/dt for i in range(days +1)]) # 
##    P_ax.set_xticklabels([str(i) for i in range(days +1)], fontsize = 12)
#    P_ax.set_xticklabels(labels, fontsize = 12)    
#    P_ax.set_xlim([0,days*24/dt])
#    P_ax.grid()   
    
    plt.figure()
    plt.plot(t, [i/0.25 for i in Energy_price], 'b', label='Energy price from grid')
    #plt.legend(loc='center right')    
    plt.grid()
    plt.xlim([0,days*24])
    plt.xticks([i*24 for i in range(days +1)], labels)
    #plt.xlabel('Time [h]')
    plt.ylabel('Energy price, $c_{E}^{grid}$, [€/kWh]')
    #plt.title('Power from the grid')
    plt.show()      

def get_Pareto(population_registry, identify_season = True, indentify_source = False):

    import matplotlib.pyplot as plt    
    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman"
    })    
    
    l_e = []
    l_t = []
    l_CO2 = []
    l_c = []
    
    for instant in population_registry:
        for individual in instant: 
            l_e.append(individual[-4])
            l_t.append(individual[-3])
            l_CO2.append(individual[-2])
            l_c.append(individual[-1])  
            
            
#    plt.figure()
#    plt.scatter(range(len(l_c)), l_t)
#    plt.xlim([0,len(l_c)])
#    plt.yticks(arange(0, 2.51, step=0.25), [round(i,2) for i in arange(0, 2.51, step=0.25)])
#    plt.grid()
#    plt.show()
            
    plt.figure()
    plt.scatter(l_e, l_t)
    plt.xlabel('Energy cost')
    plt.ylabel('Thermal comfort cost')
    plt.grid()                
    
#    fig14, (axs14) = plt.subplots(1,3)
#    axs14[0].scatter(l_e, l_t, c='b')
#    #axs14[0].set_ylim([-10, 45])
#    #axs14[0].set_title('Spring')
#    axs14[0].set_xlabel('Energy cost')
#    axs14[0].set_ylabel('Thermal comfort cost')
#    #axs14[0].legend(loc='lower center', ncol=6)
#    axs14[0].grid()
#    
#    axs14[1].scatter(l_e, l_CO2, c='b')
#    #axs14[1.set_ylim([-10, 45])
#    #axs14[1].set_title('Spring')
#    axs14[1].set_xlabel('Energy cost')
#    axs14[1].set_ylabel('CO$_{2}$ emissions cost')
#    #axs14[1].legend(loc='lower center', ncol=6)
#    axs14[1].grid()
#    
#    axs14[2].scatter(l_t, l_CO2, c='b')
#    #axs14[2.set_ylim([-10, 45])
#    #axs14[2].set_title('Spring')
#    axs14[2].set_xlabel('Thermal comfort cost')
#    axs14[2].set_ylabel('CO$_{2}$ emissions cost')
#    #axs14[2].legend(loc='lower center', ncol=6)
#    axs14[2].grid()
#    plt.show()     
    
    if indentify_source:
        c_e_0 = []
        c_t_0 = []
        c_CO2_0 = []
        c_0 = []
        i_0 = []
        
        c_e_1 = []
        c_t_1 = []
        c_CO2_1 = []
        c_1 = []
        i_1 = []
        
        c_e_2 = []
        c_t_2 = []
        c_CO2_2 = []
        c_2 = []
        i_2 = []
        
        c_e_3 = []
        c_t_3 = []
        c_CO2_3 = []
        c_3 = []
        i_3 = []
        
        c_e_4 = []
        c_t_4 = []
        c_CO2_4 = []
        c_4 = []
        i_4 = []
        
        c_e_5 = []
        c_t_5 = []
        c_CO2_5 = []
        c_5 = []
        i_5 = []
        
        c_e_6 = []
        c_t_6 = []
        c_CO2_6 = []
        c_6 = []
        i_6 = []
        
        c_e_7 = []
        c_t_7 = []
        c_CO2_7 = []
        c_7 = []
        i_7 = []
        
        i = 0
        for instant in population_registry:
            for individual in instant:
        
                if individual[0:3] == [0, 0, 0]:     # 0 - No heating
                    c_e_0.append(individual[-4])
                    c_t_0.append(individual[-3])
                    c_CO2_0.append(individual[-2])
                    c_0.append(individual[-1])
                    i_0.append(i)
                elif individual[0:3] == [0, 0, 1]:   # 1- Only HP
                    c_e_1.append(individual[-4])
                    c_t_1.append(individual[-3])
                    c_CO2_1.append(individual[-2])
                    c_1.append(individual[-1])  
                    i_1.append(i)
                elif individual[0:3] == [0, 1, 0]:   # 2- Only TESS
                    c_e_2.append(individual[-4])
                    c_t_2.append(individual[-3])
                    c_CO2_2.append(individual[-2])
                    c_2.append(individual[-1])   
                    i_2.append(i)
                elif individual[0:3] == [0, 1, 1]:   # 3 - TESS + HP
                    c_e_3.append(individual[-4])
                    c_t_3.append(individual[-3])
                    c_CO2_3.append(individual[-2])
                    c_3.append(individual[-1])    
                    i_3.append(i)
                elif individual[0:3] == [1, 0, 0]:   # 4- Only SC
                    c_e_4.append(individual[-4])
                    c_t_4.append(individual[-3])
                    c_CO2_4.append(individual[-2])
                    c_4.append(individual[-1])  
                    i_4.append(i)
                elif individual[0:3] == [1, 0, 1]:   # 5 - SC + HP
                    c_e_5.append(individual[-4])
                    c_t_5.append(individual[-3])
                    c_CO2_5.append(individual[-2])
                    c_5.append(individual[-1]) 
                    i_5.append(i)
                elif individual[0:3] == [1, 1, 0]:   # 6 - SC + TESS
                    c_e_6.append(individual[-4])
                    c_t_6.append(individual[-3])
                    c_CO2_6.append(individual[-2])
                    c_6.append(individual[-1])  
                    i_6.append(i)
                else: # elif i[0:3] == [1, 1, 1]:   # All heating activated
                    c_e_7.append(individual[-4])
                    c_t_7.append(individual[-3])
                    c_CO2_7.append(individual[-2])
                    c_7.append(individual[-1])
                    i_7.append(i)
            i +=1

#        plt.figure()
#        plt.scatter(range(0, len(c_0)), c_0, c='b', label='-')
#        plt.scatter(range(len(c_0), len(c_1)+len(c_0)), c_1, c='r', label='HP')
#        plt.scatter(range(len(c_1)+len(c_0), len(c_2)+len(c_1)+len(c_0)), c_2, c='g', label='TESS')
#        plt.scatter(range(len(c_2)+len(c_1)+len(c_0), len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_3, c='k', label='HP + TESS')
#        plt.scatter(range(len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_4, c='c', label='SC')
#        plt.scatter(range(len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_5, c='y', label='SC + HP')
#        plt.scatter(range(len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_6)+len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_6, c='b', label='SC + TESS')
#        plt.scatter(range(len(c_6)+len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_7)+len(c_6)+len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_7, c='r', label='SC + TESS + HP')
##        plt.xlim([0,len(l_c)])
#        plt.xlabel('Individuals')
##        plt.ylim([0,2])
##        plt.yticks(arange(0, 2.01, step=0.25), [round(i,2) for i in arange(0, 2.01, step=0.25)])
#        plt.ylabel('Overall cost')
#        plt.legend(loc='lower center', ncol=4)
#        plt.grid()
#        plt.show()
#        
#        
#        fig16, (axs16) = plt.subplots(1,3)
#        axs16[0].scatter(c_e_0, c_t_0, c='k', label = '-')
#        axs16[0].scatter(c_e_1, c_t_1, c='tab:orange', label = 'HP')
#        axs16[0].scatter(c_e_2, c_t_2, c='tab:red', label = 'TESS')
#        axs16[0].scatter(c_e_3, c_t_3, c='tab:pink', label = 'HP + TESS')
#        axs16[0].scatter(c_e_4, c_t_4, c='tab:gray', label = 'SC')
#        axs16[0].scatter(c_e_5, c_t_5, c='tab:olive', label = 'SC + HP')
#        axs16[0].scatter(c_e_6, c_t_6, c='tab:cyan', label = 'SC + TESS')
#        axs16[0].scatter(c_e_7, c_t_7, c='b', label = 'SC + TESS + HP')
##        axs16[0].set_xlim([-0.10, 2])
##        axs16[0].set_ylim([0, 1.2])
#        #axs16[0].set_title('Spring')
#        axs16[0].set_xlabel('Energy cost')
#        axs16[0].set_ylabel('Thermal comfort cost')
#        #axs16[0].legend(loc='lower center', ncol=6)
#        axs16[0].grid()
#        
#        axs16[1].scatter(c_e_0, c_CO2_0, c='k', label = '-')
#        axs16[1].scatter(c_e_1, c_CO2_1, c='tab:orange', label = 'HP')
#        axs16[1].scatter(c_e_2, c_CO2_2, c='tab:red', label = 'TESS')
#        axs16[1].scatter(c_e_3, c_CO2_3, c='tab:pink', label = 'HP + TESS')
#        axs16[1].scatter(c_e_4, c_CO2_4, c='tab:gray', label = 'SC')
#        axs16[1].scatter(c_e_5, c_CO2_5, c='tab:olive', label = 'SC + HP')
#        axs16[1].scatter(c_e_6, c_CO2_6, c='tab:cyan', label = 'SC + TESS')
#        axs16[1].scatter(c_e_7, c_CO2_7, c='b', label = 'SC + TESS + HP')
##        axs16[1].set_xlim([-0.10, 2])
##        axs16[1].set_ylim([0, 1.4])
#        #axs16[1].set_title('Spring')
#        axs16[1].set_xlabel('Energy cost')
#        axs16[1].set_ylabel('CO$_{2}$ emissions cost')
#        #axs16[1].legend(loc='lower center', ncol=6)
#        axs16[1].grid()
#        
#        axs16[2].scatter(c_t_0, c_CO2_0, c='k', label = '-')
#        axs16[2].scatter(c_t_1, c_CO2_1, c='tab:orange', label = 'HP')
#        axs16[2].scatter(c_t_2, c_CO2_2, c='tab:red', label = 'TESS')
#        axs16[2].scatter(c_t_3, c_CO2_3, c='tab:pink', label = 'HP + TESS')
#        axs16[2].scatter(c_t_4, c_CO2_4, c='tab:gray', label = 'SC')
#        axs16[2].scatter(c_t_5, c_CO2_5, c='tab:olive', label = 'SC + HP')
#        axs16[2].scatter(c_t_6, c_CO2_6, c='tab:cyan', label = 'SC + TESS')
#        axs16[2].scatter(c_t_7, c_CO2_7, c='b', label = 'SC + TESS + HP')
##        axs16[2].set_xlim([0, 1.2])
##        axs16[2].set_ylim([0, 1.4])
#        #axs16[2].set_title('Spring')
#        axs16[2].set_xlabel('Thermal comfort cost')
#        axs16[2].set_ylabel('CO$_{2}$ emissions cost')
#        
#        #axs16[2].legend(loc='lower center', ncol=6)
#        axs16[2].grid()
#        
#        
#        #plt.xticks(arange(0, 24*4+1, step=4*2), [i for i in range(0,25,2)])  # Set label locations. 
#        #plt.xlim([0, 24*4])
#        labels = ['-','HP','TESS','HP + TESS','SC','SC + HP','SC + TESS','SC + TESS + HP']
#        fig16.legend(labels, loc='upper center', ncol=8)
#        plt.show()
                   
    
        fig2, (axs2) = plt.subplots(2,1)   
        axs2[0].scatter(i_0, c_t_0, c='k', label = '-')
        axs2[0].scatter(i_1, c_t_1, c='tab:orange', label = 'HP')
        axs2[0].scatter(i_2, c_t_2, c='tab:red', label = 'TESS')
        axs2[0].scatter(i_3, c_t_3, c='tab:pink', label = 'HP + TESS')
        axs2[0].scatter(i_4, c_t_4, c='tab:gray', label = 'SC')
        axs2[0].scatter(i_5, c_t_5, c='tab:olive', label = 'SC + HP')
        axs2[0].scatter(i_6, c_t_6, c='tab:cyan', label = 'SC + TESS')
        axs2[0].scatter(i_7, c_t_7, c='b', label = 'SC + TESS + HP')
#        axs2[0].set_xlim([-0.10, 2])
#        axs2[0].set_ylim([0, 1.2])
        #axs2[0].set_title('Spring')
        axs2[0].set_xlabel('Individual')
        axs2[0].set_ylabel('Thermal comfort cost')
        #axs2[0].legend(loc='lower center', ncol=6)
        axs2[0].grid()  
        
  
        axs2[1].scatter(i_0, c_e_0, c='k', label = '-')
        axs2[1].scatter(i_1, c_e_1, c='tab:orange', label = 'HP')
        axs2[1].scatter(i_2, c_e_2, c='tab:red', label = 'TESS')
        axs2[1].scatter(i_3, c_e_3, c='tab:pink', label = 'HP + TESS')
        axs2[1].scatter(i_4, c_e_4, c='tab:gray', label = 'SC')
        axs2[1].scatter(i_5, c_e_5, c='tab:olive', label = 'SC + HP')
        axs2[1].scatter(i_6, c_e_6, c='tab:cyan', label = 'SC + TESS')
        axs2[1].scatter(i_7, c_e_7, c='b', label = 'SC + TESS + HP')
#        axs2[1].set_xlim([-0.10, 2])
#        axs2[1].set_ylim([0, 1.2])
        #axs2[1].set_title('Spring')
        axs2[1].set_xlabel('Individual')
        axs2[1].set_ylabel('Energy cost')
        #axs2[1].legend(loc='lower center', ncol=6)
        axs2[1].grid()           
    

        labels = ['-','HP','TESS','HP + TESS','SC','SC + HP','SC + TESS','SC + TESS + HP']
        fig2.legend(labels, loc='upper center', ncol=8)
        plt.show()

def Save_CSV(variable, file_name = ''):
    import csv

    with open(file_name, 'w') as file:
        writer = csv.writer(file)
        writer.writerows(variable)   

###############################################################################
################################   Simulations ################################

import matplotlib.pyplot as plt
import csvreader
from numpy import savetxt
from math import ceil


start_day = 200 #200 - 32 - 115 - 200
end_day = 207  #207 - 39 - 275 - 215
horizon = 8#4*2#4*1
consecutive_generations = 5
individuals = 2000
#
labels = ['19/07', '21/07', '22/07', '23/07', '24/07', '25/07', '26/07', '27/07']
#labels = ['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02']

method = 'RandomForests' # [False, 'linear_random', 'RandomForests']
list_format = 'multiple' # ['single', 'multiple']

#file_name = 'PF_H' + str(int(horizon/4)) + '_G' + str(int(consecutive_generations)) + '_P' + str(int(individuals)) + '.csv'

global n_generations
n_generations = [0]
global population_registry
population_registry = []
global best_candidate_registry
best_candidate_registry = []
global diversity_registry
diversity_registry = []
global losses_registry
losses_registry = [[],[],[]]

time_registry = []

# [SC, TESS, HP, PV, BESS, Grid]
#energy_costs = [0.15, 1.4, 0.11, 0.18, 0.13, 0.483]
energy_costs = [0, 0, 0, 0, 0.13, 0.483]
#CO2_costs = [0.071, 0.00155, 0.1, 0.041, 0.1916, 0.325]
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
G = array([CSVDataRad.ar[i][0] for i in a])
Energy_price = [i[0]*0.25 for i in CSVDataPrices.ar] # [energy_costs[-1]*0.25 for i in CSVDataPrices.ar]



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
#    print('Current step: ', i_step)
    i_step += 1
    
    start_step = time.time()    # The timer is initializad.
    h = step + horizon + 1
    
    
#    energy_costs[-1] = Energy_price[step]
    
    
    horizon_costs = [[],[]]
    for i in range(horizon + 1):
        energy_costs[-1] = Energy_price[step + horizon]
        horizon_costs[0].append(energy_costs)        
        horizon_costs[1].append(CO2_costs)
    
    
    population_registry.append([])
    best_candidate_registry.append([])    
#    from numpy import random
#    SC_active = random.randint(0,2)
#    TESS_active = random.randint(0,2)
#    HP_active = random.randint(0,2)
#    PV_curtail = random.uniform(0,1)
#    P_BESS_active = random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1]))
    
    
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


#    P_BESS_active = random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1]))
#    SoC_BESS_state = update_BESS(SoC_BESS[-1], P_BESS[-1], P_Load[-1] + P_HP[-1]/1000)
#    SoC_BESS = append(SoC_BESS[-1], SoC_BESS_state)    
    
#    return [Qdot_SC, Qdot_TESS, Qdot_HP, Qdot_losses, T_TESS, T_in, P_HP, P_PV, P_BESS, P_Grid, SoC_BESS]
    
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


#    SC_active = random.randint(0,2)
#    TESS_active = random.randint(0,2)
#    HP_active = random.randint(0,2)
#    PV_curtail = random.uniform(0,1)
#    P_BESS_active = random.uniform(BESS_perm_min(SoC_BESS[-1]), BESS_perm_max(SoC_BESS[-1]))

    end_step = time.time()    # The timer is initializad.    
    time_registry.append(end_step - start_step)    

    ts += 1
    if ts == 96:
        ts = 0
        day += 1
        print('Current day: ', start_day + day, '. Current temperature: ', T_in[-1]-273)
    



end = time.time()    # The timer is initializad.
totalelapsed = end - start  # The total time is calculated.

from statistics import mean
print('Time per day: ', totalelapsed/(end_day - start_day), ', T_av: ', mean(T_in)-273, 'T_max', max(T_in)-273)


week_graph_TQP(t[0:7*24*4], T_in[2:2+7*24*4], T_amb[0:7*24*4], T_set[0:7*24*4], Qdot_TESS[1:1+7*24*4], Qdot_SC[2:2+7*24*4], Qdot_HP[2:2+7*24*4], Qdot_Losses[2:2+7*24*4], P_Grid[0:7*24*4], P_Load[0:7*24*4], P_HP[2:2+7*24*4], P_PV[0:7*24*4], P_BESS[0:7*24*4], labels, Energy_price[start_day*24*4+1:end_day*4*24+1])

from statistics import median
print('Grid purchase: ', sum(P_Grid)*dt/3600, ' kWh')
print('Grid costs', ' €', (sum([P*c for P,c in zip(P_Grid, Energy_price[start_day*24*4+1:end_day*4*24+1])]))*dt/3600)
print('Energy cost: ', ' €', (sum([P*c for P,c in zip(P_Grid, Energy_price[start_day*24*4+1:end_day*4*24+1])]) + sum([i for i in P_BESS if i>0])*energy_costs[-2])*dt/3600)
print('Thermal cost: ', beta*sum([abs(T_set[i] - T_in[i+2]) for i in range(len(T_in)-2)]))
print('CO2 cost: ', sum(P_Grid)*CO2_costs[-1]*dt/3600, ' gCO_{2}')
print('Median time: ', median(time_registry), ' s')


extras = True

if extras:

    ###################################   CSVs  ###################################
    
    
#    T_in_file_name = 'T_in_H' + str(int(horizon/4)) + '_G' + str(int(consecutive_generations)) + '_P' + str(int(individuals)) + '.csv'
#    savetxt(T_in_file_name, T_in, delimiter =", ", fmt ='% s')
#    #savetxt('Thermal_load_15min.csv', Qdot_Losses, delimiter =", ", fmt ='% s')
#    
#    results = zip(t.tolist(), Qdot_SC[2:].tolist(), Qdot_SC_TESS.tolist(), Qdot_TESS[1:].tolist(), Qdot_HP[2:].tolist(), Qdot_Losses[2:].tolist(), T_TESS.tolist(), T_in[2:].tolist(), P_HP[2:].tolist(), P_PV.tolist(), P_BESS.tolist(), P_Grid.tolist(), SoC_BESS.tolist(), n_generations)
#    
#    if 'file_name' in globals():
#        import csv
#        with open(file_name, "w") as f:
#            writer = csv.writer(f)
#            for row in results:
#                writer.writerow(row)
#            f.close()
    
    #savetxt(file_name, results, delimiter =", ", fmt ='% s')
    
    
    ##################################   Plots  ###################################
    
    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman"
    })
    
    
    plt.figure()
    plt.plot(t, [i-273 for i in T_TESS])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Temperature in the TESS, $T_{TESS}$, [°C]')
    #plt.title('Temperature in the TESS')
    plt.show()
    
    plt.figure()
    plt.plot(t, [i/1000 for i in Qdot_TESS[1:]])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Thermal power of the TESS, $\dot{Q}_{TESS}$, [kW]')
    #plt.title('Thermal power provided by the TESS')
    plt.show()
    
    plt.figure()
    plt.plot(t, [i-273 for i in T_amb[0:-24*4]], 'r', label='Ambient temperature, $T_{amb}$')
    plt.plot(t, [i-273 for i in T_in[2:]], 'b', label='Temperature inside the house, $T_{in}$')
    #plt.plot(t, [i-273 for i in T_set[start_day*24*4:end_day*24*4]], 'g', label='Setpoint temperature inside the house, $T_{set}$')
    plt.legend(loc='lower center')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Temperature [°C]')
    #plt.title('Temperature inside the house')
    plt.show()    
    
    plt.figure()
    plt.plot(t, [i/1000 for i in Qdot_Losses[2:]])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylim([0, 1.8])
    plt.ylabel('Thermal losses, $\dot{Q}_{L}$, [kW]')
    #plt.title('Thermal losses in the house')
    plt.show()    
    
    plt.figure()
    plt.plot(t, [i/1000 for i in Qdot_HP[2:]])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Thermal power of the HP, $\dot{Q}_{HP}$, [kW]')
    #plt.title('Heat provided by the HP')
    plt.show()
    
    plt.figure()
    plt.plot(t, [i/1000 for i in Qdot_SC[2:]], 'b', label='To the thermal load')
    #plt.plot(t, [i/1000 for i in Qdot_SC_TESS], 'r', label='To the TESS')
    plt.legend(loc='lower center')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Thermal power of the SC, $\dot{Q}_{SC}$, [kW]')
    #plt.title('Thermal power provided by the SC')
    plt.show()
    #
    plt.figure()
    plt.plot(t, P_Load[start_day*4*24:end_day*4*24])
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Electric load, $P_{L}$, [kW]')
    #plt.title('Electric load')
    plt.show()  
    #
    #plt.figure()
    #plt.plot(t, [i/1000 for i in P_HP[2:]])
    #plt.grid()
    #plt.xlim([0, end_day*24])
    #plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    ##plt.xlabel('Time [h]')
    #plt.ylabel('Power of the HP, $P_{HP}$, [kW]')
    ##plt.title('Power consumed by the HP')
    #plt.show()
    
    plt.figure()
    plt.plot(t, P_PV, 'r', label='PV power curtailed')
    plt.plot(t, P_PV_av[start_day*24*4:end_day*4*24], 'b', label='PV power available')
    plt.legend(loc='lower center')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('PV power load, $P_{PV}$, [kW]')
    #plt.title('Electric load')
    plt.show()  
    
    plt.figure()
    plt.plot(t, SoC_BESS)
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('State-of-charge of the BESS, $SoC_{BESS}$ [%]')
    #plt.title('SoC of the BESS')
    plt.show()
    
    plt.figure()
    plt.plot(t, P_BESS, 'b', label='Power delivered by the BESS')
    #plt.legend(loc='center right')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Power of the BESS, $P_{BESS}$, [kW]')
    #plt.title('Power of the BESS')
    plt.show()
    
    plt.figure()
    plt.plot(t, P_Grid, 'b', label='Power consumed from the grid')
    #plt.legend(loc='center right')
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Power from the grid, $P_{G}$, [kW]')
    #plt.title('Power from the grid')
    plt.show()


    plt.figure()
    plt.plot(t, [i/0.25 for i in Energy_price[start_day*24*4+1:end_day*4*24+1]], 'b', label='Power consumed from the grid')
    #plt.legend(loc='center right')    
    plt.grid()
    plt.xlim([0, end_day*24])
    plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    #plt.xlabel('Time [h]')
    plt.ylabel('Energy price, $c_{E}^{grid}$, [€/kWh]')
    #plt.title('Power from the grid')
    plt.show()    
    
    
    
    ###############################################################################
    # Pareto fronts
    day= 1
    d_hour = 12
    
    instant = (day-1)*24*4 + d_hour*4
    
    c_e_0 = []
    c_t_0 = []
    c_CO2_0 = []
    c_0 = []
    
    c_e_1 = []
    c_t_1 = []
    c_CO2_1 = []
    c_1 = []
    
    c_e_2 = []
    c_t_2 = []
    c_CO2_2 = []
    c_2 = []
    
    c_e_3 = []
    c_t_3 = []
    c_CO2_3 = []
    c_3 = []
    
    c_e_4 = []
    c_t_4 = []
    c_CO2_4 = []
    c_4 = []
    
    c_e_5 = []
    c_t_5 = []
    c_CO2_5 = []
    c_5 = []
    
    c_e_6 = []
    c_t_6 = []
    c_CO2_6 = []
    c_6 = []
    
    c_e_7 = []
    c_t_7 = []
    c_CO2_7 = []
    c_7 = []
    
    for i in population_registry[instant]:
        if i[0:3] == [0, 0, 0]:     # 0 - No heating
            c_e_0.append(i[-4])
            c_t_0.append(i[-3])
            c_CO2_0.append(i[-2])
            c_0.append(i[-1])  
        elif i[0:3] == [0, 0, 1]:   # 1- Only HP
            c_e_1.append(i[-4])
            c_t_1.append(i[-3])
            c_CO2_1.append(i[-2])
            c_1.append(i[-1])  
        elif i[0:3] == [0, 1, 0]:   # 2- Only TESS
            c_e_2.append(i[-4])
            c_t_2.append(i[-3])
            c_CO2_2.append(i[-2])
            c_2.append(i[-1])          
        elif i[0:3] == [0, 1, 1]:   # 3 - TESS + HP
            c_e_3.append(i[-4])
            c_t_3.append(i[-3])
            c_CO2_3.append(i[-2])
            c_3.append(i[-1])    
        elif i[0:3] == [1, 0, 0]:   # 4- Only SC
            c_e_4.append(i[-4])
            c_t_4.append(i[-3])
            c_CO2_4.append(i[-2])
            c_4.append(i[-1])   
        elif i[0:3] == [1, 0, 1]:   # 5 - SC + HP
            c_e_5.append(i[-4])
            c_t_5.append(i[-3])
            c_CO2_5.append(i[-2])
            c_5.append(i[-1])  
        elif i[0:3] == [1, 1, 0]:   # 6 - SC + TESS
            c_e_6.append(i[-4])
            c_t_6.append(i[-3])
            c_CO2_6.append(i[-2])
            c_6.append(i[-1])  
        else: # elif i[0:3] == [1, 1, 1]:   # All heating activated
            c_e_7.append(i[-4])
            c_t_7.append(i[-3])
            c_CO2_7.append(i[-2])
            c_7.append(i[-1])          
            
    
    l = population_registry[instant]
    l_e = [i[-4] for i in l]
    l_t = [i[-3] for i in l]
    l_CO2 = [i[-2] for i in l]
    l_c = [i[-1] for i in l]
    
    plt.figure()
    plt.scatter(range(len(l_c)), l_c)
    plt.xlim([0,len(l_c)])
#    plt.ylim([0,2])
#    plt.yticks(arange(0, 2.51, step=0.25), [round(i,2) for i in arange(0, 2.51, step=0.25)])
    plt.grid()
    plt.show()
    
    fig14, (axs14) = plt.subplots(1,3)
    axs14[0].scatter(l_e, l_t, c='b')
    #axs14[0].set_ylim([-10, 45])
    #axs14[0].set_title('Spring')
    axs14[0].set_xlabel('Energy cost')
    axs14[0].set_ylabel('Thermal comfort cost')
    #axs14[0].legend(loc='lower center', ncol=6)
    axs14[0].grid()
    
    axs14[1].scatter(l_e, l_CO2, c='b')
    #axs14[1.set_ylim([-10, 45])
    #axs14[1].set_title('Spring')
    axs14[1].set_xlabel('Energy cost')
    axs14[1].set_ylabel('CO$_{2}$ emissions cost')
    #axs14[1].legend(loc='lower center', ncol=6)
    axs14[1].grid()
    
    axs14[2].scatter(l_t, l_CO2, c='b')
    #axs14[2.set_ylim([-10, 45])
    #axs14[2].set_title('Spring')
    axs14[2].set_xlabel('Thermal comfort cost')
    axs14[2].set_ylabel('CO$_{2}$ emissions cost')
    #axs14[2].legend(loc='lower center', ncol=6)
    axs14[2].grid()
    
    
    #plt.xticks(arange(0, 24*4+1, step=4*2), [i for i in range(0,25,2)])  # Set label locations. 
    #plt.xlim([0, 24*4])
    #labels = ['Surface', '$y = 10$ cm', '$y = 20$ cm', '$y = 30$ cm', '$y = 40$ cm', '$y = 50$ cm',]
    #fig3.legend(labels, loc='lower center', ncol=6)
    plt.show()
    
    ###############################################################################
    plt.figure()
    plt.scatter(range(0, len(c_0)), c_0, c='b', label='-')
    plt.scatter(range(len(c_0), len(c_1)+len(c_0)), c_1, c='r', label='HP')
    plt.scatter(range(len(c_1)+len(c_0), len(c_2)+len(c_1)+len(c_0)), c_2, c='g', label='TESS')
    plt.scatter(range(len(c_2)+len(c_1)+len(c_0), len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_3, c='k', label='HP + TESS')
    plt.scatter(range(len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_4, c='c', label='SC')
    plt.scatter(range(len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_5, c='y', label='SC + HP')
    plt.scatter(range(len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_6)+len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_6, c='b', label='SC + TESS')
    plt.scatter(range(len(c_6)+len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0), len(c_7)+len(c_6)+len(c_5)+len(c_4)+len(c_3)+len(c_2)+len(c_1)+len(c_0)), c_7, c='r', label='SC + TESS + HP')
    plt.xlim([0,len(l_c)])
    plt.xlabel('Individuals')
#    plt.ylim([0,2])
#    plt.yticks(arange(0, 2.01, step=0.25), [round(i,2) for i in arange(0, 2.01, step=0.25)])
    plt.ylabel('Overall cost')
    plt.legend(loc='lower center', ncol=4)
    plt.grid()
    plt.show()
    
    
    fig16, (axs16) = plt.subplots(1,3)
    axs16[0].scatter(c_e_0, c_t_0, c='k', label = '-')
    axs16[0].scatter(c_e_1, c_t_1, c='tab:orange', label = 'HP')
    axs16[0].scatter(c_e_2, c_t_2, c='tab:red', label = 'TESS')
    axs16[0].scatter(c_e_3, c_t_3, c='tab:pink', label = 'HP + TESS')
    axs16[0].scatter(c_e_4, c_t_4, c='tab:gray', label = 'SC')
    axs16[0].scatter(c_e_5, c_t_5, c='tab:olive', label = 'SC + HP')
    axs16[0].scatter(c_e_6, c_t_6, c='tab:cyan', label = 'SC + TESS')
    axs16[0].scatter(c_e_7, c_t_7, c='b', label = 'SC + TESS + HP')
#    axs16[0].set_xlim([-0.10, 2])
#    axs16[0].set_ylim([0, 1.2])
    #axs16[0].set_title('Spring')
    axs16[0].set_xlabel('Energy cost')
    axs16[0].set_ylabel('Thermal comfort cost')
    #axs16[0].legend(loc='lower center', ncol=6)
    axs16[0].grid()
    
    axs16[1].scatter(c_e_0, c_CO2_0, c='k', label = '-')
    axs16[1].scatter(c_e_1, c_CO2_1, c='tab:orange', label = 'HP')
    axs16[1].scatter(c_e_2, c_CO2_2, c='tab:red', label = 'TESS')
    axs16[1].scatter(c_e_3, c_CO2_3, c='tab:pink', label = 'HP + TESS')
    axs16[1].scatter(c_e_4, c_CO2_4, c='tab:gray', label = 'SC')
    axs16[1].scatter(c_e_5, c_CO2_5, c='tab:olive', label = 'SC + HP')
    axs16[1].scatter(c_e_6, c_CO2_6, c='tab:cyan', label = 'SC + TESS')
    axs16[1].scatter(c_e_7, c_CO2_7, c='b', label = 'SC + TESS + HP')
#    axs16[1].set_xlim([-0.10, 2])
#    axs16[1].set_ylim([0, 1.4])
    #axs16[1].set_title('Spring')
    axs16[1].set_xlabel('Energy cost')
    axs16[1].set_ylabel('CO$_{2}$ emissions cost')
    #axs16[1].legend(loc='lower center', ncol=6)
    axs16[1].grid()
    
    axs16[2].scatter(c_t_0, c_CO2_0, c='k', label = '-')
    axs16[2].scatter(c_t_1, c_CO2_1, c='tab:orange', label = 'HP')
    axs16[2].scatter(c_t_2, c_CO2_2, c='tab:red', label = 'TESS')
    axs16[2].scatter(c_t_3, c_CO2_3, c='tab:pink', label = 'HP + TESS')
    axs16[2].scatter(c_t_4, c_CO2_4, c='tab:gray', label = 'SC')
    axs16[2].scatter(c_t_5, c_CO2_5, c='tab:olive', label = 'SC + HP')
    axs16[2].scatter(c_t_6, c_CO2_6, c='tab:cyan', label = 'SC + TESS')
    axs16[2].scatter(c_t_7, c_CO2_7, c='b', label = 'SC + TESS + HP')
#    axs16[2].set_xlim([0, 1.2])
#    axs16[2].set_ylim([0, 1.4])
    #axs16[2].set_title('Spring')
    axs16[2].set_xlabel('Thermal comfort cost')
    axs16[2].set_ylabel('CO$_{2}$ emissions cost')
    
    #axs16[2].legend(loc='lower center', ncol=6)
    axs16[2].grid()
    
    
    #plt.xticks(arange(0, 24*4+1, step=4*2), [i for i in range(0,25,2)])  # Set label locations. 
    #plt.xlim([0, 24*4])
    labels = ['-','HP','TESS','HP + TESS','SC','SC + HP','SC + TESS','SC + TESS + HP']
    fig16.legend(labels, loc='upper center', ncol=8)
    plt.show()
    
    
    #plt.figure()
    #plt.plot(t[0:-1], time_registry, 'b', label='Cumputation time')
    ##plt.legend(loc='center right')
    #plt.grid()
    #plt.xlim([0, end_day*24])
    #plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
    ##plt.xlabel('Date')
    #plt.ylabel('Cumputation time')
    ##plt.title('Cumputation time')
    
    
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
    
#    from collections import Counter
#    c2 = Counter(diversity_registry)
#    plt.figure(19)
#    #plt.hist(diversity_registry, density=True)
#    plt.bar(c2.keys(), c2.values())
#    plt.xlabel('Population diversity')
#    plt.ylabel('Frequency')
#    #plt.title('Population diversity')
#    plt.show()
#    
#    diveristy_file_name = 'Diversity_H' + str(int(horizon/4)) + '_G' + str(int(consecutive_generations)) + '_P' + str(int(individuals)) + '.csv'
#    savetxt(diveristy_file_name, diversity_registry, delimiter =", ", fmt ='% s')
    
    