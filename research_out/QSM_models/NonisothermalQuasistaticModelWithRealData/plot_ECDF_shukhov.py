import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
from matplotlib.widgets import Slider

# Установка рабочей директории
project_path = r'C:\Users\Nikita Andrianov\Documents\GitHub\pde_solvers\research_out\QSM_models\NonisothermalQuasistaticModelWithRealData'
os.chdir(project_path)
print("Рабочая директория установлена:", os.getcwd())

# Получаем список всех папок
folders = [folder for folder in os.listdir()]
print("Найденные элементы в текущей папке:")
print(folders)

# Название файла, который мы ищем в папках
filename = 'diff_temp.csv'

# Список возможных экспериментов
experiments_type = {
    'Выбор задания реологии в стационарной модели' : ['StationaryInitialReology', 'StationaryCurrentReology', 'StationaryMeanReology'],
    'Сравнение стационарной и квазистационарной модели' : ['QuasiStationaryFullReology', 'StationaryCurrentReology'],
    'Сравнение модели без кт и идентифицированным кт' : ['QuasiStationaryFullReology', 'QuasiStationaryFullReologyIdeal'],
    'Исследование влияния плотности и вязкости на квазистац' : ['QuasiStationaryDensityOnly', 'QuasiStationaryFullReology', 'QuasiStationaryViscosityOnly'],
    'Сравнение методов расчета температуры' : ['QuasiStationaryFullReology']
}

# Главный цикл
while True:
    ch_dict = {}
    for i, research_name in enumerate(experiments_type.keys()):
        ch_dict[i] = research_name
        print(f'{i+1}. {research_name}')
    try:
        choice = int(input('Выберите эксперимент: ')) - 1

        dfs = []

        if ch_dict[choice] == 'Сравнение методов расчета температуры':
            folder = 'QuasiStationaryFullReology'
            df_r = pd.read_csv(os.path.join(os.getcwd(), folder, filename), encoding='windows-1251')

            df_r.columns = ['time1', 'time2', 'etalon', 'calculated', 'diff_temp', 'calculated_shukhov', 'temp_delta_shukhov']

            df1 = df_r['diff_temp']
            df2 = df_r['temp_delta_shukhov']

            df1.name = 'Осн. метод'
            df2.name = 'Шухов'
            labels = ['Погрешность осн. метод', 'Погрешность Шухова']

            dfs = [df1.to_frame(), df2.to_frame()]
        else:
            for folder in folders:
                if folder in experiments_type[ch_dict[choice]]:
                    df_r = pd.read_csv(os.path.join(os.getcwd(), folder, filename), encoding='windows-1251')
                    df_r.columns.name = folder
                    dfs.append(df_r)

        parameters_names = [df.columns.tolist()[0 if ch_dict[choice] == 'Сравнение методов расчета температуры' else 4] for df in dfs]

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
                axes[i].set_xlabel(f'Погрешность температуры ({labels[i]}), К', fontsize=12)
                axes[i].set_ylabel(dfs[i].columns.name)
                axes[i].set_xlim(-10, 10)

        def draw_fun(p):
            for i in range(len(parameters_names)):
                res = stats.ecdf(dfs[i][parameters_names[i]])
                quant = res.cdf.quantiles
                prob = res.cdf.probabilities

                [imin, imax] = find_prob_x(quant, prob, p, (1 - p) / 2)
                print(f'Квантильные значения: imin = {imin}, imax = {imax}')

                fsize = 20
                msize = 5

                x_min, x_max = axes[i].get_xlim()

                axes[i].plot(quant, prob, 'k-', linewidth=2)
                axes[i].text(imin - (x_max - x_min) * 0.05, 0.08, f'{imin:.1f}', fontsize=fsize)
                axes[i].text(imax, 0.90, f'{imax:.1f}', fontsize=fsize)
                axes[i].text(x_min + (x_max - x_min) * 0.05, 0.5, f'delta = {(imax - imin):.3f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
                axes[i].text(x_min + (x_max - x_min) * 0.05, p, f'p = {p:.2f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})

                axes[i].axhline(y=(1 - p) / 2, color='r', linestyle='--')
                axes[i].axhline(y=(1 + p) / 2, color='r', linestyle='--')
                axes[i].plot(imin, (1 - p) / 2, 'bo', markersize=msize)
                axes[i].plot(imax, (1 + p) / 2, 'bo', markersize=msize)

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

    except Exception as e:
        print('Ошибка выбора!')
        print('Подробности ошибки:', e)
        import traceback
        traceback.print_exc()
        continue
