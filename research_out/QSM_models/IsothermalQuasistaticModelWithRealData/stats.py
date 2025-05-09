import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import levene, kruskal

def compare_samples(compare_data, file_name='diff_press.csv', alpha=0.05):
    """
    Сравнивает несколько выборок и добавляет вывод о значимости различий.
    
    Параметры:
        compare_data (dict): Словарь с данными.
        file_name (str): Название файла с данными.
        alpha (float): Уровень значимости (по умолчанию 0.05).
    
    Возвращает:
        pd.DataFrame: Таблица с результатами и выводами.
    """
    results = []
    
    for exp_name, groups in compare_data.items():
        # Загрузка данных
        data = {}
        for group_name, folder in groups.items():
            df = pd.read_csv(f'{folder}/{file_name}')
            data[group_name] = df['diff_press'].dropna().values
        
        if len(data) < 2:
            raise ValueError(f"В эксперименте '{exp_name}' должно быть минимум 2 выборки!")
        
        # Анализ средних и СКО
        mean_abs = {name: np.abs(np.mean(values)) for name, values in data.items()}
        best_mean_group = min(mean_abs, key=mean_abs.get)
        mean_diff = {name: val - min(mean_abs.values()) for name, val in mean_abs.items()}
        
        std_values = {name: np.std(values, ddof=1) for name, values in data.items()}
        best_std_group = min(std_values, key=std_values.get)
        std_diff = {name: val - min(std_values.values()) for name, val in std_values.items()}
        
        # Проверка значимости
        if len(data) == 2:
            stat_mean, p_mean = stats.mannwhitneyu(*data.values(), alternative='two-sided')
            test_name_mean = 'Манна-Уитни'
            mean_significant = p_mean < alpha
        else:
            stat_mean, p_mean = kruskal(*data.values())
            test_name_mean = 'Краскел-Уоллис'
            mean_significant = p_mean < alpha
        
        stat_std, p_std = levene(*data.values(), center='median')
        test_name_std = 'Браун-Форсайт'
        std_significant = p_std < alpha
        
        # Формирование выводов
        mean_conclusion = "Различия значимы" if mean_significant else "Различия незначимы"
        std_conclusion = "Различия значимы" if std_significant else "Различия незначимы"
        
        # Запись результатов
        for group in data.keys():
            results.append({
                'Эксперимент': exp_name,
                'Группа': group,
                'Мат. ожидание': np.mean(data[group]),
                '|Мат. ожидание|': mean_abs[group],
                'Отличие мат. ожидания от лучшего': mean_diff[group],
                'СКО': std_values[group],
                'Отличие СКО от лучшего': std_diff[group],
                'Лучшая группа по среднему': best_mean_group,
                'Лучшая группа по СКО': best_std_group,
                'Тест для средних': test_name_mean,
                'p-value (средние)': p_mean,
                'Вывод о средних': mean_conclusion,
                'Тест для СКО': test_name_std,
                'p-value (СКО)': p_std,
                'Вывод о СКО': std_conclusion
            })
    
    return pd.DataFrame(results)

# Пример использования
compare_data = {
    'Текщая реология VS Начальная реология': {
        'Начальная реология': 'StationaryInitialReology',
        'Текущая реология': 'StationaryCurrentReology'
    }, 
    'Текщая реология VS Средняя реология': {
        'Текщая реология': 'StationaryCurrentReology',
        'Средняя реология': 'StationaryMeanReology'
    }, 
    'Начальная реология VS Средняя реология': {
        'Начальная реология': 'StationaryInitialReology',
        'Средняя реология': 'StationaryMeanReology'
    }, 
    'Стац VS Квазистац': {
        'Стац': 'StationaryCurrentReology',
        'Квазистац': 'QuasiStationaryFullReology'
    },
    'Квазистац на плотности VS Квазистац на вязкости': {
        'Квазистац на плотности': 'QuasiStationaryDensityOnly',
        'Квазистац на вязкости': 'QuasiStationaryViscosityOnly'
    },
    'Квазистац VS Квазистац на плотности': {
        'Квазистац на плотности': 'QuasiStationaryDensityOnly',
        'Квазистац': 'QuasiStationaryFullReology'
    },
    'Квазистац VS Квазистац на вязкости': {
        'Квазистац на вязкости': 'QuasiStationaryViscosityOnly',
        'Квазистац': 'QuasiStationaryFullReology'
    },
}

results_df = compare_samples(compare_data)
print(results_df)

# Сохранение в CSV
results_df.to_excel('stat_comparison_with_conclusions.xlsx', index=False)