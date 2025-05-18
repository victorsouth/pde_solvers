import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator

file_name = 'ident_history.csv'

# folder = 'DiameterWithPrint' # указать нужную папку
folder = 'FrictionWithPrinter' # указать нужную папку

# Чтение CSV в DataFrame
df = pd.read_csv(f'{folder}/{file_name}')  

# Создание фигуры с тремя подграфиками (3 строки, 1 столбец)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 12))

# График 1: target_function от step
axes[0].plot(df['step'], df['target_function'], label='Целевая функция', color='blue')
# axes[0].set_xlabel('Steps')
axes[0].set_xticks(df['step'])
axes[0].set_ylabel('L(Δθ)', fontsize=20)
axes[0].tick_params(axis='both', labelsize=20)
axes[0].grid(True)
axes[0].legend()

# График 2: argument_history от step
axes[1].set_yticks([tick/10 for tick in range(int(df['argument_history'].min() - 1) * 10, int(df['argument_history'].max() + 1) * 10)])
axes[1].plot(df['step'], df['argument_history'], label=f"Поправочный коэффициент на {'λ' if folder == 'FrictionWithPrinter' else 'd'}", color='green')
axes[1].set_xlabel('Steps', fontsize=20)
axes[1].set_xticks(df['step'] if folder == 'FrictionWithPrinter' else df['step'])
# axes[1].set_yticks([i/100 for i in range(9999990, 110, 10)])
axes[1].yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
axes[1].yaxis.get_major_formatter().set_scientific(False)
axes[1].yaxis.get_major_formatter().set_useOffset(False)
axes[1].yaxis.set_major_locator(MultipleLocator(0.00005))  # не более 5 отметок
axes[1].set_ylabel('θfric' if folder == 'FrictionWithPrinter' else 'θdiam', fontsize=20)
axes[1].grid(True)
axes[1].tick_params(axis='both', labelsize=20)
axes[1].legend()

# # График 3: steps от step 
# axes[2].plot(df['step'], df['steps'], label='Steps', color='red')
# axes[2].set_title('Steps')
# axes[2].set_xlabel('Step')
# axes[2].set_ylabel('Steps')
# axes[2].grid(True)
# axes[2].legend()

# Автоматическая регулировка расстояний между графиками
plt.tight_layout()

# Показать графики
plt.show()