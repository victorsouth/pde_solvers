import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
from matplotlib.widgets import Slider

# folders = ['/' + folder + '/' for folder in os.listdir()]
project_path = r'C:\Users\Nikita Andrianov\Documents\GitHub\pde_solvers\research_out\QSM_models\NonisothermalQuasistaticModelWithRealData'  # для Windows

os.chdir(project_path)
print("Рабочая директория установлена:", os.getcwd())
folders = [folder for folder in os.listdir()]
print("Найденные элементы в текущей папке:")
print(folders)
filename = '/diff_temp.csv'

experiments_type = {
    'Выбор задания реологии в стационарной модели' : ['StationaryInitialReology', 'StationaryCurrentReology', 'StationaryMeanReology'],
    'Сравнение стационарной и квазистационарной модели' : ['QuasiStationaryFullReology', 'StationaryCurrentReology'],
    'Сравнение модели без кт и идентифицированным кт' : ['QuasiStationaryFullReology', 'QuasiStationaryFullReologyIdeal'],
    'Исследование влияния плотности и вязкости на квазистац' : ['QuasiStationaryDensityOnly', 'QuasiStationaryFullReology', 'QuasiStationaryViscosityOnly']
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
                #df_r['diff_temp'] = df_r['diff_temp'] / 1000.0
                df_r.columns.name = folder
                dfs.append(df_r)
        print(dfs)

        parameters_names = [df.columns.tolist()[4] for df in dfs]

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
                axes[i].set_xlabel('Погрешность температуры в конце ЛУ, К')
                print(dfs[i].columns.name)
                axes[i].set_ylabel(dfs[i].columns.name)
                axes[i].set_xlim(-10, 10)

        #def ecdf(data):
        #    x = np.sort(data)
        #    y = np.arange(1, len(x)+1) / len(x)
        #    return x, y

        def draw_fun(p):
            for i in range(len(parameters_names)):
                res = stats.ecdf(dfs[i][parameters_names[i]])
                axes[i].ecdf(dfs[i][parameters_names[i]])
                #x, y = ecdf(dfs[i][parameters_names[i]])

                q = (1 - p) / 2

                [imin, imax] = find_prob_x(res.cdf.quantiles, res.cdf.probabilities, p, q)
                print(f'Квантильные значения: imin = {imin}, imax = {imax}')
                #imin, imax = find_prob_x(x, y, p, q)

                fsize = 20
                msize = 5

                x_min, x_max = axes[i].get_xlim()  # получить текущие границы X

                axes[i].text(imin - (x_max - x_min) * 0.05, 0.08, f'{imin:.1f}', fontsize=fsize)
                axes[i].text(imax, 0.90, f'{imax:.1f}', fontsize=fsize)
                axes[i].text(x_min + (x_max - x_min) * 0.05, 0.5, f'delta = {(imax - imin):.3f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
                axes[i].text(x_min + (x_max - x_min) * 0.05, p, f'p = {p:.2f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
                
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
    except:
        print('Ошибка выбора!')
        continue
