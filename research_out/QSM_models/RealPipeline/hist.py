import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


path = 'QuasiStationaryFull'

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

data = pd.read_csv(path + '/' + file, encoding='windows-1251')
name = data.columns.tolist()[2]

bins_count = int((data[name].max() - data[name].min()) / 30000)

plt.hist(data[name], bins=bins_count, color='skyblue', edgecolor='black')
plt.grid(visible=True)
#Смещение: {data[name].skew()}\n
print(f'Медиана: {data[name].median()}\nСреднее: {data[name].mean()}')
print(f'СКО: {data[name].std()}')

plt.show()

