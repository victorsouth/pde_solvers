import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
from matplotlib.widgets import Slider

# folders = ['/' + folder + '/' for folder in os.listdir()]
folders = [folder for folder in os.listdir()]
filename = '/diff_press.csv'

dfs = []
df_r = pd.read_csv(os.getcwd() + '/PipeIdentificationWithPrinter' + filename, encoding='windows-1251')
df_r['diff_press_before'] = df_r['diff_press_before'] / 1000.0
df_r['diff_press_after'] = df_r['diff_press_after'] / 1000.0
df_before = df_r[df_r.columns[:2].to_list() + ['diff_press_before']]
df_before.columns.name = 'Погрешность до идентификации'
df_after = df_r[df_r.columns[:2].to_list() + ['diff_press_after']]
df_after.columns.name = 'Погрешность после идентификации'

dfs = [df_before, df_after]
print(dfs)

parameters_names = [df.columns.tolist()[2] for df in dfs]

interval = 30

def find_prob_x(quant, prob, p, q):
    min_flag = True
    for i in range(len(prob)):
        if prob[i] >= q and min_flag:
            i_min = i
            min_flag = False
        if prob[i] > p + q:
            i_max = i - 1
            return quant[i_min], quant[i_max]

axes = [plt.subplot(1, len(parameters_names), _ + 1) for _ in range(len(parameters_names))]

    
def init_func():
    for i in range(len(axes)):
        axes[i].clear()
        axes[i].grid(visible=True)
        axes[i].set_xlabel('Погрешность давления в конце ЛУ, кПа')
        print(dfs[i].columns.name)
        axes[i].set_ylabel(dfs[i].columns.name)
        axes[i].set_xlim(-200, 200)

def draw_fun(p):
    for i in range(len(parameters_names)):
        res = stats.ecdf(dfs[i][parameters_names[i]])
        axes[i].ecdf(dfs[i][parameters_names[i]])

        q = (1 - p) / 2

        [imin, imax] = find_prob_x(res.cdf.quantiles, res.cdf.probabilities, p, q)

        fsize = 11
        msize = 5
        axes[i].text(imin - 50, 0.08, f'{imin:.1f}', fontsize=fsize)
        axes[i].text(imax , 0.90, f'{imax:.1f}', fontsize=fsize)
        axes[i].text(-180, 0.5, f'delta = {(imax - imin):.1f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
        axes[i].text(-190, p, f'p = {p:.2f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
        
        axes[i].axhline(y = q, color = 'r', linestyle = '--') 
        axes[i].axhline(y = p + q, color = 'r', linestyle = '--') 
        axes[i].plot(imin, q, 'bo', markersize=msize)
        axes[i].plot(imax, p + q, 'bo', markersize=msize)
        

def change(pe):
    init_func()
    draw_fun(pe)

init_func()

draw_fun(0.90)

ax_time = plt.axes([0.15, 0.0001, 0.5, 0.04])
time_slider = Slider(ax_time, 'p', 0, 1, valstep=0.05, valinit=0.90)

time_slider.on_changed(change)    

plt.subplots_adjust(left=0.05, bottom=0.15, right=0.976, top=0.967, 
        wspace=0.2, hspace=0.136)

plt.show()
