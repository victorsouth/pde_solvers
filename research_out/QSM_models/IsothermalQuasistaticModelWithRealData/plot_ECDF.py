import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
from matplotlib.widgets import Slider

# folders = ['/' + folder + '/' for folder in os.listdir()]
folders = [folder for folder in os.listdir()]
filename = '/diff_press.csv'

experiments_type = {
    'Выбор задания реологии в стационарной модели Тек+Начальная' : ['StationaryCurrentReology','StationaryInitialReology', ],
    'Выбор задания реологии в стационарной модели Тек+Ср' : ['StationaryCurrentReology', 'StationaryMeanReology'],
    'Выбор задания реологии Full' : ['StationaryCurrentReology','StationaryInitialReology', 'StationaryMeanReology'],
    'Сравнение стационарной и квазистационарной модели' : ['QuasiStationaryFullReology', 'StationaryCurrentReology'],
    'Исследование влияния плотности и вязкости на квазистац FULL' : ['QuasiStationaryDensityOnly', 'QuasiStationaryFullReology', 'QuasiStationaryViscosityOnly'],
    'Исследование влияния плотности и вязкости на квазистац RHOvsVISC' : ['QuasiStationaryDensityOnly','QuasiStationaryViscosityOnly'],
    'Исследование влияния плотности и вязкости на квазистац KVvsVSIC' : ['QuasiStationaryFullReology', 'QuasiStationaryViscosityOnly']
}
 
while True:
    ch_dict = {}
    for i, research_name in enumerate(experiments_type.keys()):
        ch_dict[i] = research_name
        print(f'{i+1}. {research_name}')
    try:
        choice = int(input('Выберите эксперимент: ')) - 1


        dfs = []
        for folder in folders:
            if folder in experiments_type[ch_dict[choice]]:
                df_r = pd.read_csv(os.getcwd() + '/' + folder + filename, encoding='windows-1251')
                df_r['diff_press'] = df_r['diff_press'] / 1000.0
                df_r.columns.name = folder
                dfs.append(df_r)
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
                axes[i].set_xlabel('Погрешность давления в конце ЛУ, кПа', fontsize=20)
                print(dfs[i].columns.name)
                axes[i].set_ylabel(dfs[i].columns.name, fontsize=20)
                axes[i].set_xlim(-200, 200)

        def draw_fun(p):
            for i in range(len(parameters_names)):
                res = stats.ecdf(dfs[i][parameters_names[i]])
                axes[i].ecdf(dfs[i][parameters_names[i]])

                q = (1 - p) / 2

                [imin, imax] = find_prob_x(res.cdf.quantiles, res.cdf.probabilities, p, q)

                fsize = 20
                msize = 5
                axes[i].text(imin - 50, 0.08, f'{imin:.1f}', fontsize=fsize)
                axes[i].text(imax , 0.90, f'{imax:.1f}', fontsize=fsize)
                axes[i].text(-180, 0.5, f'delta = {(imax - imin):.1f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
                axes[i].text(-190, p, f'p = {p:.2f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
                axes[i].tick_params(axis='both', labelsize=20)
                axes[i].axhline(y = q, color = 'r', linestyle = '--') 
                axes[i].axhline(y = p + q, color = 'r', linestyle = '--') 
                axes[i].plot(imin, q, 'bo', markersize=msize)
                axes[i].plot(imax, p + q, 'bo', markersize=msize)
                

        def change(pe):
            init_func()
            draw_fun(pe)

        init_func()

        draw_fun(0.95)

        ax_time = plt.axes([0.15, 0.0001, 0.5, 0.04])
        time_slider = Slider(ax_time, 'p', 0, 1, valstep=0.05, valinit=0.95)

        time_slider.on_changed(change)    

        plt.subplots_adjust(left=0.05, bottom=0.15, right=0.976, top=0.967, 
                wspace=0.2, hspace=0.136)

        plt.show()
    except:
        print('Ошибка выбора!')
        continue
