import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import tabulate
import math

experiments_type = {
    'Выбор задания реологии в стационарной модели' : ['StationaryInitialReology', 'StationaryCurrentReology', 'StationaryMeanReology'],
    'Сравнение стационарной и квазистационарной модели' : ['QuasiStationaryFullReology', 'StationaryCurrentReology'],
    'Исследование влияния плотности и вязкости на квазистац' : ['QuasiStationaryDensityOnly', 'QuasiStationaryFullReology', 'QuasiStationaryViscosityOnly']
}

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

x_right = dfs[0][parameters_names[0]].max()
x_left = dfs[0][parameters_names[0]].min()
y_bot = 0
y_top = 0
interval = 30


axes = [plt.subplot(len(parameters_names), 1, _ + 1) for _ in range(len(parameters_names))]

def find_y_top(df, name, n):
    df['bins'] = pd.cut(df[name], bins=n)
    freqs = df.groupby(['bins'])[name].count().reset_index()
    return freqs[name].max()

for name, df in zip(parameters_names, dfs):
    if x_right < df[name].max():
        x_right = df[name].max()
    if x_left > df[name].min():
        x_left = df[name].min()

    bins_count = int((df[name].max() - df[name].min()) / interval)
    freq = find_y_top(df, name, bins_count)
    if freq > y_top:
        y_top = freq
    
def init_func():
    global dfs
    for i in range(len(axes)):
        axes[i].clear()
        axes[i].grid(visible=True)
        axes[i].set_xlabel('Погрешность давления в конце ЛУ, кПа')
        print(dfs[i].columns.name)
        axes[i].set_ylabel(dfs[i].columns.name)
        axes[i].set_xlim(x_left, x_right)
        axes[i].set_ylim(y_bot, y_top + 100)

def draw_fun():
    global dfs
    for i in range(len(parameters_names)):
        bins_count = int((dfs[i][parameters_names[i]].max() - dfs[i][parameters_names[i]].min()) / interval)
        top_pos = y_top - 400
        axes[i].hist(dfs[i][parameters_names[i]], bins=bins_count, color='skyblue', edgecolor='black')
        axes[i].text(x_right - 200, top_pos, f'СКО: {dfs[i][parameters_names[i]].std():.4f}', fontsize=12, bbox={'facecolor': 'white', 'alpha': 1})
        axes[i].text(x_right - 200, top_pos - 500, f'Среднее: {dfs[i][parameters_names[i]].mean():.4f}', fontsize=12, bbox={'facecolor': 'white', 'alpha': 1})
        dfs[i]['sum_square'] = dfs[i][parameters_names[i]].apply(lambda x: x ** 2)
        axes[i].text(x_right - 200, top_pos - 900, f'ЦФ: {math.sqrt((dfs[i]["sum_square"].sum()) / len(dfs[i]["sum_square"])):.4f}', fontsize=12, bbox={'facecolor': 'white', 'alpha': 1})
init_func()

draw_fun()

plt.subplots_adjust(left=0.05, bottom=0.06, right=0.976, top=0.967, 
        wspace=0.2, hspace=0.136)

plt.show()


