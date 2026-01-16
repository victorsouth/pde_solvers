#pragma once

#ifndef M_G
/// @brief Ускорение свободного падения
#define M_G 9.81
#endif
#ifndef M_R
/// @brief Газовая постоянная (более точное значение 8.31446261815324);
#define M_R 8.31441 
#endif
#ifndef KELVIN_OFFSET
/// @brief Разность между градусами Кельвина и Цельсия
#define KELVIN_OFFSET 273.15
#endif
#ifndef ATMOSPHERIC_PRESSURE
/// @brief Стандартное давление 
#define ATMOSPHERIC_PRESSURE 101325
#endif
#ifndef DEFAULT_AMBIENT_TEMPERATURE
/// @brief Внешняя температура (дефолтная)
#define DEFAULT_AMBIENT_TEMPERATURE 300
#endif
#ifndef TECHNICAL_ATHMOSPERIC_PRESSURE
/// @brief Техническая атмосфера
#define TECHNICAL_ATHMOSPERIC_PRESSURE 98067
#endif
#ifndef NORMAL_TEMPERATURE
/// @brief Нормальная температура 
#define NORMAL_TEMPERATURE KELVIN_OFFSET + 20
#endif
#ifndef STANDART_TEMPERATURE
/// @brief Стандартная лабораторная температура KELVIN_OFFSET + 25
#define STANDART_TEMPERATURE 298.15
#endif
#ifndef Kv2CvCoeff
/// @brief Коэффициент пересчёта между европейским и американским коэффициентами
/// пропускной способности
#define Kv2CvCoeff 1.17
#endif
#ifndef PLANCK_CONST
/// @brief Постоянная Планка
#define PLANCK_CONST 6.626176e-34
#endif
#ifndef BOLTZMANN_CONST
/// @brief Постоянная Больцмана
#define BOLTZMANN_CONST 1.380662e-23
#endif
#ifndef AVOGADRO_CONST
/// @brief Постоянная Авогадро
#define AVOGADRO_CONST 6.02214E+23
#endif
#ifndef CALORIE_JOULE
/// @brief  Количество джоулей в одной калории.
/// Почему столько? Всегда же 4.184 было? Ответ: это старое «международное» определение.
/// Для термохимических расчётов следует использовать 4.184.
#define CALORIE_JOULE 4.1868
#endif

