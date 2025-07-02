import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox
import os

path = './ColdLU_DynamicTemperature/'

# Ожидаем, что пользователь выберет файл из списка
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
    
    # Загрузка данных
    rawData = pd.read_csv(f'{path}/{file}', encoding='windows-1251')

    # Названия столбцов
    columns = rawData.columns.tolist()
    [_, timeLabel] = columns[:2]  # Используем второй time
    etalon_col = 'T_out_etalon'
    calc_col = 'T_out_calc'
    diff_col = 'T_out_error'
    
    # Один временной вектор (второй столбец)
    time = rawData[timeLabel]

    # Значения
    etalon = rawData[etalon_col]
    calculated = rawData[calc_col]
    diff_temp = rawData[diff_col]
    
    # Построение графиков
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Верхний график: Эталон и обе модели
    axs[0].plot(time, etalon, color='red', label='Эталонные значения')
    axs[0].plot(time, calculated, color='blue', label='Вычисленные значения (основная модель)')
    axs[0].set_ylabel('Температура')
    axs[0].legend()
    axs[0].grid(True)
    axs[0].set_title('Сравнение эталонных и вычисленных температур')

    # Нижний график: Отклонения обеих моделей
    axs[1].plot(time, diff_temp, color='green', label='Отклонения (основная модель)')
    axs[1].set_xlabel('Время')
    axs[1].set_ylabel('Δ Температура')
    axs[1].legend()
    axs[1].grid(True)
    axs[1].set_title('Отклонения моделей')

    plt.tight_layout()
    plt.show()