import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import os

plots = [file for file in os.listdir() if '.py' not in file]
plots.sort()

titles = {
    'rare' : 'Разреженный исходный профиль',
    'norm' : 'Нормальный исходный профиль',
    'tight' : 'Плотный исходный профиль'
    }

axes = [plt.subplot(int(len(plots) / 2), 1, _ + 1) for _ in range(int(len(plots) / 2))]

for i in range(0, len(plots), 2):
    for n, t in titles.items():
        if n in plots[i]:
            axes_index = int(i / 2)
            df_prof = pd.read_csv(plots[i], encoding='windows-1251')
            axes[axes_index].plot(df_prof['km'], df_prof['height'], 'bo-', markersize=6, label = 'Исходный профиль')    
            df_prof = pd.read_csv(plots[i + 1], encoding='windows-1251')
            axes[axes_index].plot(df_prof['km'], df_prof['height'], 'ro-', markersize=3, label = 'Новый профиль')
            axes[axes_index].set_title(t)
            axes[axes_index].set_xlabel('Координаты, м')
            axes[axes_index].set_ylabel('Высотки, м')
            axes[axes_index].legend(loc=4)
            axes[axes_index].grid()
        
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.971, top=0.954, 
         wspace=0.465, hspace=0.24)
        
plt.show()
            
    
    
    