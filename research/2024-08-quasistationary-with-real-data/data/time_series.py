import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob


plot_names = glob.glob('*.csv')
if 'coord_heights.csv' in plot_names:
    plot_names.remove('coord_heights.csv')

rawData = list()
for plot in plot_names:
    df = pd.read_csv(f'{plot}', header = 0, names = ['time', 'value'], encoding='windows-1251', sep=';')
    df['time'] = pd.to_datetime(df['time'].values, format = "%d.%m.%Y %H:%M:%S")
    df['value'] = df['value'].str.replace(',', '.').astype(float)

    rawData.append(df)

plotsCount = len(plot_names) 

fig = plt.figure(figsize=(8, 5))

axes = [plt.subplot(plotsCount, 1, _ + 1) for _ in range(plotsCount)]

def init_func():
    for i in range(len(axes)):
        axes[i].clear()
        axes[i].grid(visible=True)
        axes[i].set_xlabel("Время, ч")
        axes[i].set_ylabel(plot_names[i])
        data = rawData[i]
        coordData = data['time'].tolist()
        paramData = data['value'].tolist()
        #xLim = [(min(coordData) - 0.1 * (max(coordData) - min(coordData)))/3600, (max(coordData) + 0.1 * (max(coordData) - min(coordData)))/3600]
        #xLim = [coordData., coordData.]
        yLim = [min(paramData) - 0.1 * (max(paramData) - min(paramData)), max(paramData) + 0.1 * (max(paramData) - min(paramData))]
        #axes[i].set_xlim(xLim)
        axes[i].set_ylim(yLim)
    ini_draw()
        
plots = list()
def ini_draw(step=0):
    for i in range(plotsCount):
        data = rawData[i]
        coordData = data['time'].tolist()
        paramData = data['value'].tolist()
        plots.append(axes[i].plot(coordData, paramData, 'ro', markersize=1))

init_func()

plt.subplots_adjust(left=0.06, bottom=0.06, right=0.971, top=0.973, 
         wspace=0.2, hspace=0.24)

plt.show()


dummy = 1