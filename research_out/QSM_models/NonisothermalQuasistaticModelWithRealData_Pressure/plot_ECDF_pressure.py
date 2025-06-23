import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
from matplotlib.widgets import Slider

# Установка рабочей директории
project_path = r'C:\Users\Nikita Andrianov\Documents\GitHub\pde_solvers\research_out\QSM_models\NonisothermalQuasistaticModelWithRealData_Pressure'
os.chdir(project_path)
print("Рабочая директория установлена:", os.getcwd())

# Получаем список всех папок
folders = [folder for folder in os.listdir()]
print("Найденные элементы в текущей папке:")
print(folders)

# Название CSV-файла
filename = 'diff_temp.csv'

# Возможные эксперименты
experiments_type = {
    'Выбор задания реологии в стационарной модели' : ['StationaryInitialReology', 'StationaryCurrentReology', 'StationaryMeanReology'],
    'Сравнение стационарной и квазистационарной модели' : ['QuasiStationaryFullReology', 'StationaryCurrentReology'],
    'Сравнение модели без кт и идентифицированным кт' : ['QuasiStationaryFullReology', 'QuasiStationaryFullReologyIdeal'],
    'Исследование влияния плотности и вязкости на квазистац' : ['QuasiStationaryDensityOnly', 'QuasiStationaryFullReology', 'QuasiStationaryViscosityOnly'],
    'Сравнение методов расчета температуры' : ['QuasiStationaryFullReology_Pressure']
}

# Функция поиска границ интервала
def find_prob_x(quant, prob, p, q):
    min_flag = True
    for i in range(len(prob)):
        if prob[i] >= q and min_flag:
            i_min = i
            min_flag = False
        if prob[i] > p + q:
            i_max = i - 1
            return quant[i_min], quant[i_max]

# Главный цикл
while True:
    ch_dict = {}
    for i, research_name in enumerate(experiments_type.keys()):
        ch_dict[i] = research_name
        print(f'{i+1}. {research_name}')
    try:
        choice = int(input('Выберите эксперимент: ')) - 1
        selected_folders = experiments_type[ch_dict[choice]]

        # Найти первый подходящий файл
        for folder in folders:
            if folder in selected_folders:
                df = pd.read_csv(os.path.join(os.getcwd(), folder, filename), encoding='windows-1251')
                data = df.iloc[:, 6] / 1000.0 # 6-й столбец
                title = folder
                break
        else:
            raise Exception("Не найден подходящий CSV-файл")

        fig, ax = plt.subplots()

        def init_func():
            ax.clear()
            ax.grid(visible=True)
            ax.set_xlabel('Погрешность давления в конце ЛУ, кПа', fontsize=12)
            ax.set_ylabel(title)
            ax.set_ylim(0, 1)
            ax.set_xlim(-200, 200)

        def draw_fun(p):
            sorted_data = np.sort(data)
            n = len(sorted_data)
            prob = np.linspace(1/n, 1, n)

            q = (1 - p) / 2
            imin, imax = find_prob_x(sorted_data, prob, p, q)
            print(f'Квантильные значения: imin = {imin}, imax = {imax}')

            init_func()

            ax.plot(sorted_data, prob, 'k-', linewidth=2)

            fsize = 14
            msize = 6
            x_min, x_max = ax.get_xlim()

            ax.text(imin - 50, 0.08, f'{imin:.1f}', fontsize=fsize)
            ax.text(imax , 0.90, f'{imax:.1f}', fontsize=fsize)
            ax.text(-180, 0.5, f'delta = {(imax - imin):.1f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
            ax.text(-190, p, f'p = {p:.2f}', fontsize=fsize, bbox={'facecolor': 'white', 'alpha': 1})
                
            ax.axhline(y = q, color = 'r', linestyle = '--') 
            ax.axhline(y = p + q, color = 'r', linestyle = '--') 
            ax.plot(imin, q, 'bo', markersize=msize)
            ax.plot(imax, p + q, 'bo', markersize=msize)


        def change(pe):
            draw_fun(pe)

        draw_fun(0.90)

        ax_time = plt.axes([0.15, 0.0001, 0.5, 0.04])
        time_slider = Slider(ax_time, 'p', 0, 1, valstep=0.05, valinit=0.90)
        time_slider.on_changed(change)

        plt.subplots_adjust(bottom=0.15)
        plt.show()

    except Exception as e:
        print('Ошибка выбора!')
        print('Подробности ошибки:', e)
        import traceback
        traceback.print_exc()
        continue
