# Intervalos de confianza para un modelo de regresión no lineal

El presente repositorio contiene tres métodos para calcular los intervalos de confianza para la respuesta media de un modelo no lineal, se emplearon la aproximación de Bates & Watts (1988, p.59), el método delta y la propagación del error por medio de una simulación Monte Carlo (Spiess, 2018).


### Método delta

El método delta es un enfoque numérico utilizado para estimar la varianza aproximada de una función que depende de una o más variables aleatorias. Es una aproximación de primer orden de Taylor, que permite que la función sea descrita por una recta tangente que la interseca en su media. 

Es un método muy empleado para estimar los intervalos de confianza de la respuesta media cuando se emplean estimadores de máxima verosimilitud. También puede usarse a partir de estimadores de mínimos cuadrados ordinarios lineales y no lineales. El método delta toma una función compleja, continua, y diferenciable, a partir de la cual, estima una aproximación lineal de la función, y calcula su varianza.

En el código se estima el vector gradiente para cada una de las derivadas parciales, en función de los estimadores paramétricos, para luego multiplicarlo por la matriz de varianzas y covarianzas de los estimadores. Finalmente, la adición o sustracción del producto de la raíz cuadrada del gradiente, para cada estimación, por el cuantil t, constituye el intervalo de confianza para la respuesta media.


### Aproximación Bates & Watts (1988, p. 59)

La aproximación de Bates & Watts (1988) , para regresión no lineal emplea la descomposición QR para la matriz V0, que corresponde a la derivada parcial de las observaciones respecto a cada parámetro del modelo, usando los estimados óptimos. Una vez hecha la descomposición QR, se recuperan Q1 y R1, en la cual a R1 se le calcula la inversa para así obtener los intervalos de confianza.



### Propagación del error por Simulación Monte Carlo

Este método consiste en la propagación del error de la variable predictora. Se usa como entrada la matriz de varianzas y covarianzas de los estimadores paramétricos, y su valor estimado como media. A partir de estos valores, se genera una simulación normal multivariada de n muestras, para cada observación i de la variable predictora (Spiess, 2013). Es decir, se reproducen combinaciones de estimadores paramétricos simulados siguiendo la matriz de varianzas y covarianzas estimada, y la tendencia central de los &beta;<sub>i</sub> estimados. Además, se asume una distribución normal de la variable dependiente por cada i (Spiess, 2013). 

Con cada grupo de estimadores paramétricos simulados se calcula la función de la variable respuesta para cada observación i. Finalmente, se obtienen estadísticos como la media, desviación estándar y cuantiles asociados al intervalo de confianza por cada observación i (Spiess, 2013). 

  
  
### Referencias 


Bates, D. M., & Watts, D. G. (1988). Nonlinear regression analysis and its applications (Vol. 2). New York: Wiley.

Spiess, A. (2013). predictNLS (Part 1, Monte Carlo simulation): confidence intervals for ‘nls’ models. R-bloggers. Obtenido de https://www.r-bloggers.com/2013/08/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/

Spiess, A. (2018). propagate: Propagation of Uncertainty. R package version 1.0-6. Obtenido de https://CRAN.R-project.org/package=propagate


