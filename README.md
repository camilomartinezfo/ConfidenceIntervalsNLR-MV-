# Intervalos de confianza para un modelo de regresión no lineal

El presente repositorio contiene tres métodos para calcular el intervalo de confianza para la respuesta media de un modelo no lineal, siendo las aproximaciones del método delta, Bates & Watts (1988, p.59) y propagación del error por medio de una simulación Monte Carlo (Spiess, 2018)


### Propagación del error por Simulación Monte Carlo

Este consiste en la propagación del error de la variable predictora, en la cual se usa como entrada la matriz de varianzas y covarianzas de los parámetros ajustados del modelo y sus parámetros estimados como media para crear una simulación normal multivariada de n muestras para cada observación i de la variable predictora (Spiess, 2013). Es decir, se generan combinaciones de parámetros simulados siguiendo la matriz de varianzas y covarianzas estimada siguiendo la tendencia central de los parámetros ajustados, por lo tanto, es importante notar que se asume una distribución normal de la variable dependiente por cada i (Spiess, 2013). 

Con cada grupo de parámetros simulados se calcula la función de la variable respuesta para cada observación i. Finalmente, se obtienen estadísticos como la media, desviación estándar y cuantiles asociados al intervalo de confianza por cada observación i (Spiess, 2013). 


 

