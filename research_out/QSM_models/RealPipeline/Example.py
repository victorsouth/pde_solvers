import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox
import os

path = 'QuasiStationaryFull'
while True:
    plots = os.listdir(path)
    while True:
        for i in range(len(plots)):
            print(f'{i+1}. {plots[i]}')
        try:
            choose = input("Выберите файл: ")
            if len(choose) != 1:
                print("неверный ввод")
                continue
            choice = int(choose) - 1
            file = plots[choice]
            break
        except:
            print("неверный ввод")
                

    rawData = pd.read_csv(f'{path}/{file}', encoding='windows-1251')
    parametersNames = rawData.columns.tolist()[2:]
    [timeLabel, coordLabel] = rawData.columns.tolist()[:2]
    plotsCount = len(parametersNames) 

    fig = plt.figure(figsize=(8, 5))

    axes = [plt.subplot(plotsCount, 1, _ + 1) for _ in range(plotsCount)]

    data_skip = len(set(rawData[coordLabel]))
    time_moments = list(set(rawData[timeLabel]))
    time_moments.sort()

    def init_func():
        for i in range(len(axes)):
            axes[i].clear()
            axes[i].grid(visible=True)
            axes[i].set_xlabel(coordLabel)
            axes[i].set_ylabel(parametersNames[i])
            coordData = rawData[coordLabel]
            paramData = rawData[parametersNames[i]]
            xLim = [min(coordData) - 0.1 * (max(coordData) - min(coordData)), max(coordData) + 0.1 * (max(coordData) - min(coordData))]
            yLim = [min(paramData) - 0.1 * (max(paramData) - min(paramData)), max(paramData) + 0.1 * (max(paramData) - min(paramData))]
            axes[i].set_xlim(xLim)
            axes[i].set_ylim(yLim)
        ini_draw()
            
    plots = list()
    def ini_draw(step=0):
        for i in range(len(parametersNames)):
            coordData = rawData[coordLabel]
            paramData = rawData[parametersNames[i]]
            plots.append(axes[i].plot(coordData[step * data_skip: (step + 1) * data_skip], paramData[step * data_skip: (step + 1) * data_skip], color='b'))

    def change(step):
        for i in range(len(axes)):
            paramData = rawData[parametersNames[i]]
            plots[i][0].set_ydata(paramData[step * data_skip: (step + 1) * data_skip])
            text_box.set_val(f"{time_moments[step]}")

    init_func()

    if len(time_moments) != 1:
        ax_time = plt.axes([0.15, 0.0001, 0.5, 0.04])
        time_slider = Slider(ax_time, 'Steps', 0, len(time_moments)-1, valstep=1)
        ax_textbox = fig.add_axes([0.75, 0.001, 0.2, 0.04])
        text_box = TextBox(ax_textbox, "t =", textalignment="center")
        text_box.set_val("0")

        time_slider.on_changed(change)

    plt.subplots_adjust(bottom=0.12)
    #plt.tight_layout()

    plt.show()


