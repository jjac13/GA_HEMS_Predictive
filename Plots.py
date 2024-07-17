# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0

###############################################################################
###########################   Imported modules   ##############################
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from numpy import arange, array, zeros, linspace, ones
import csvreader


def get_file_names(folder_path = 'Files\\'):
    from os import listdir
    names = listdir(folder_path)
    
    return names


def plot_boxplot_from_files(folder_path = 'Files\\', labels = [], normalized = True, fliers = False):
    boxplot_data = []
    names = get_file_names(folder_path)
    
    for file in names:
        file_data = csvreader.read_data(csv=file, address=folder_path)
        file_data.data2cols()
        if normalized:       
            boxplot_data.append([100*i/max(file_data.col[0]) for i in file_data.col[0]])
        else:
            boxplot_data.append(file_data.col[0])
 

    boxplot_dictionary = {labels[i]: boxplot_data[i] for i in range(len(labels))}
        


    fig, ax = plt.subplots()
    ax.boxplot(boxplot_dictionary.values(), showfliers=fliers)
    ax.set_xticklabels(boxplot_dictionary.keys())
    ax.set_xlabel('Population size')
    ax.set_ylabel('Diversity, [%]')

#    return boxplot_dictionary

labels=['10', '25', '40', '75', '100', '200']
folder_path = 'Diversity\\'

##############

start_day = 0
end_day = 365

t = arange(0, end_day*24, step=0.25)

CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
CSVDataTamb.data2array()
T_amb = [i[0] for i in CSVDataTamb.ar[start_day*24*4:end_day*24*4]]

CSVDataTin_GA = csvreader.read_data(csv='T_in_dev_GA.csv', address='')
CSVDataTin_GA.data2array()
Tin_GA = [i[0] for i in CSVDataTin_GA.ar[start_day*24*4:end_day*24*4]]

CSVDataTin_RB = csvreader.read_data(csv='T_in_dev_RB.csv', address='')
CSVDataTin_RB.data2array()
Tin_RB = [i[0] for i in CSVDataTin_RB.ar[start_day*24*4:end_day*24*4]]

plt.rcParams.update({
#    "text.usetex": True,
    "font.family": "Times New Roman"
})


#plt.figure(1)
##plt.plot(t, [i for i in T_amb], 'r', label='Ambient temperature, $T_{amb}$')
#plt.plot(t, [i-273 for i in Tin_GA], 'b', label='With GA, $T_{in}$')
#plt.plot(t, [i-273 for i in Tin_RB], 'g', label='With RB, $T_{in}$')
#plt.legend(loc='upper left')
#plt.grid()
#plt.xlim([0, end_day*24])
#plt.xticks(arange(1, 8760, step=24*31), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])  # Set label locations. 
##plt.xlabel('Time [h]')
#plt.ylabel('Temperature inside the house, $T_{in}$, [°C]')
##plt.title('Temperature inside the house')
#plt.show()  



data = [Tin_GA, Tin_RB]

bins =linspace(-1, 4, 50)
plt.style.use('seaborn-deep')
plt.figure()
plt.hist([Tin_GA, Tin_RB], bins, label=['GA', 'RB'], weights=[ones(len(i))/ len(i) for i in data])
plt.legend(loc='upper right')
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlim([-1, 4])
plt.xlabel('Temperature deviation, $\Delta$T, [°C]')
plt.ylabel('Frequency')
plt.show()

fig7, ax7 = plt.subplots()
bplot1 = ax7.boxplot(data)
plt.xticks([1, 2], ['GA', 'RB'])
ax7.set_ylabel('Temperature deviation [°C]')    
plt.show()

#plt.xlabel('Time [h]')
#plt.ylabel('Temperature inside the house, $T_{in}$, [°C]')
#plt.title('Temperature inside the house')


labels = ['PV', 'BESS', 'Grid', 'SC', 'TESS', 'HP']
#GA_util = [1108, 300, 2221, 1271, 448, 945]
GA_util = [1316, 104, 1988, 1196, 591, 921]
RB_util = [1511, 464, 2045, 69, 1020, 1124]
x_axis = arange(len(labels))

plt.figure()
#plt.grid()
plt.bar(x_axis -0.2, GA_util, width=0.4, label = 'GA', color = 'blue')
plt.bar(x_axis +0.2, RB_util, width=0.4, label = 'RB', color = 'green')
plt.xticks(x_axis, labels)
plt.ylabel('Energy used [kWh]')
plt.legend()
plt.show()  



method = (
    "GA",
    "RB",
)
#weight_counts_energy = {
#    "PV": array([199.44, 271.98]),
#    "BESS": array([39, 60.32]),
#    "Grid": array([1072.74, 987.74]),
#    "SC": array([190.65, 10.35]),
#    "TESS": array([627.2, 1428]),    
#    "HP": array([103.95, 123.64]),
#}

weight_counts_energy = {
    "PV": array([0, 0]),
    "BESS": array([13.58, 60.32]),
    "Grid": array([960.01, 987.74]),
    "SC": array([0, 0]),
    "TESS": array([0, 0]),    
    "HP": array([0, 0]),
}
width = 0.5

plt.figure()

fig, ax = plt.subplots()
bottom = zeros(2)

for boolean, weight_counts_energy in weight_counts_energy.items():
    p = ax.barh(method, weight_counts_energy, width, label=boolean, left=bottom)
    bottom += weight_counts_energy

#ax.set_title("")
plt.xlabel('Accumulated cost [€]')    
ax.legend(loc="center", ncol=6)
#plt.grid()
plt.show()


weight_counts_CO2 = {
    "PV": array([45.43, 61.95]),
    "BESS": array([57.48, 88.9]),
    "Grid": array([721.83, 664.63]),
    "SC": array([91, 4.94]),
    "TESS": array([0.69, 1.58]),    
    "HP": array([94.50, 112.4]),
}

plt.figure()

fig, ax = plt.subplots()
bottom = zeros(2)

for boolean, weight_counts_CO2 in weight_counts_CO2.items():
    p = ax.barh(method, weight_counts_CO2, width, label=boolean, left=bottom)
    bottom += weight_counts_CO2

#ax.set_title("")
plt.xlabel('Accumulated CO$_{2}$ emmissions [kg$_{CO_{2}}$]')    
ax.legend(loc="center", ncol=6)
#plt.grid()
plt.show()


Cost = csvreader.read_data(csv='cost_population.csv', address='Diversity\\')
Cost.data2cols()

plt.figure()
plt.scatter(Cost.col[0], Cost.col[1], c='tab:red', label='Normalized time')
plt.scatter(Cost.col[0], Cost.col[2], c='tab:olive', label='Normalized overall cost')
plt.scatter(Cost.col[0], Cost.col[3], c='tab:blue', label='Normalized shoulder')
plt.xlim([0,max(Cost.col[0])])
plt.xlabel('Population size')
plt.ylim([0,1.1])
#plt.yticks(arange(0, 2.01, step=0.25), [round(i,2) for i in arange(0, 2.01, step=0.25)])
#plt.ylabel('Overall cost')
plt.legend(loc='lower center', ncol=1)
plt.grid()
plt.show()

