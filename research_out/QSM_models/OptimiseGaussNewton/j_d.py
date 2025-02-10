import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import tabulate

df_r = pd.read_csv(os.getcwd() + '/PipeIdentificationWithPrinter/j.csv', encoding='windows-1251')

plt.plot(df_r['d'], df_r['J'], 'bo')

plt.grid()
plt.xlabel('d')
plt.ylabel('J')
plt.show()