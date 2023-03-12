#pragma once


/// @brief Гидравлическое сопротивление по Шифринсону
/// \param reynolds_number
/// \param relative_roughness
/// \return
inline double hydraulic_resistance_shifrinson(double reynolds_number, double relative_roughness)
{
    return 0.11 * pow(relative_roughness, 0.25);
}

/// @brief Гидравлическое сопротивление по Альтушулю
/// @param reynolds_number 
/// @param relative_roughness 
/// @return 
inline double hydraulic_resistance_altshul(double reynolds_number, double relative_roughness)
{
    return 0.11 * pow(relative_roughness + 68 / reynolds_number, 0.25);
}


/// @brief Расчет гидравлического сопротивления в широком диапазоне чисел Рейнольдса
/// Квадратичное трение по формуле Исаева [Морозова, Коршак, ф-ла (1)]
/// @param reynolds_number 
/// @param relative_roughness 
/// @return 
inline double hydraulic_resistance_isaev(double reynolds_number, double relative_roughness) {
    const double Re = fabs(reynolds_number);
    const double& Ke = relative_roughness;

    double lam;

    // Формулы по РД-75.180.00-КТН-258-10. 
    // В конце нет форулы для больших чисел Рейнольдса

    if (Re < 1) {
        lam = 64; // Стокс при Re = 1
    }
    else if (Re < 2320)
    {
        lam = 64 / Re; // Стокс 
    }
    else if (Re < 4000)
    {
        // Стокс + Блазиус, сглаженный переход
        double gm = 1 - exp(-0.002 * (Re - 2320));
        lam = 64 / Re * (1 - gm) + 0.3164 / pow(Re, 0.25) * gm;
    }
    //else if(Re < 1e4)//min(1e5,27/pow(Ke,1.143))  27/pow(Ke,1.143)   1e4   //10/Ke
    //{
    //    lam=0.3164/pow(Re,0.25);
    //}
    else if (Re < 560 / Ke)
    {
        // Исаев по [Морозова, Коршак], ф-ла (1)
        lam = 1.0 / sqr(-1.8 * log10(6.8 / Re + pow(Ke / 3.7, 1.1)));
    }
    else
    {
        // Шифринсон
        lam = 0.11 * pow(Ke, 0.25);
    }

    return lam;
}


