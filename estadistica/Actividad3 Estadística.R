
### Master Universitario en Bioinformática
### Estadistica y R para la Salud 
### Actividad grupal 3: Análisis de un caso práctico en R 


# cargar las librerias necesarias
library(stats)   # librería para el PCA
library(ggplot2) # librería para hacer la representación gráfica
library(gtsummary) # para dibujar las tablas

# Cargar los datos 
Dataset <- read.csv("C:/Users/CESAR HINOJOSA/Downloads/Dataset5.csv")

# Escalar los datos 
Data_scaled <- scale(Dataset5)

# Verificar si las columnas han sido correctamente estandarizadas
summary(Data_scaled)

# Aplicar PCA
Pca_result <- prcomp(Data_scaled, center = TRUE, scale. = TRUE)

# Ver el resumen del PCA
summary(Pca_result)



### Graficos de la varianza acumulada y explicada de cada componente ###

# Obtener la varianza explicada por cada componente
explained_variance <- summary(Pca_result)$importance[2,]

# Convertir la varianza explicada en un data frame
variance_df <- data.frame(
  Componente = 1:length(explained_variance),
  varianza_explicada = explained_variance
)
print(variance_df)

# Graficar la varianza explicada
ggplot(variance_df, aes(x = Componente, y = varianza_explicada)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Scree Plot - Varianza Explicada por Componente Principal",
    x = "Componente Principal",
    y = "Varianza Explicada"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Calcular la varianza acumulada
cumulative_variance <- cumsum(explained_variance)

# Crear un dataframe con la varianza acumulada
cumulative_df <- data.frame(
  Componente = 1:length(cumulative_variance),
  Varianza_Acumulada = cumulative_variance
)

# Graficar la varianza acumulada 
ggplot(cumulative_df, aes(x = Componente, y = Varianza_Acumulada)) +
  geom_line(color = "darkgreen", size = 1.2) +
  geom_point(color = "darkorange", size = 2) +
  theme_minimal() +
  labs(
    title = "Varianza Acumulada de los Componentes Principales",
    x = "Componente Principal",
    y = "Varianza Acumulada"
    
  )


### Gráfico de dispersión de los primeros dos componentes principales  ###

# Mostrar los componentes necesarios para explicar al menos el 70% de la varianza
components_to_keep <- which(cumulative_variance >= 0.70)[1]
components_to_keep

# Crear un dataframe con los componentes principales
Pca_data <- data.frame(Pca_result$x)

# Graficar los primeros dos componentes principales con ggplot2 
ggplot(Pca_data, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, color = "blue") +
  labs(title = "Dispersión de los Primeros Dos Componentes Principales", x = "Componente Principal 1", y = "Componente Principal 2") +
  theme_minimal()



### Tabla resumen con la varianza explicada y la varianza acumulada por cada componente ###

# Crear una tabla con la varianza explicada y acumulada por cada componente
variance_table <- data.frame(
  Componente = 1:length(explained_variance),
  Varianza_Explicada = explained_variance,
  Varianza_Acumulada = cumulative_variance
)
# Ver tabla resumen
print(variance_table)



### tabla descriptiva con la media y sd de la varianza acumulada y la varianza explicada ###

# Calcular media y desviación estándar de la varianza acumulada y explicada
mean_explained_variance <- mean(variance_table$Varianza_Explicada)
sd_explained_variance <- sd(variance_table$Varianza_Explicada)
mean_cumulative_variance <- mean(variance_table$Varianza_Acumulada)
sd_cumulative_variance <- sd(variance_table$Varianza_Acumulada)

# Agregar la media y la desviación estándar calculadas a un resumen
cat("\nTabla Resumen de Varianza Explicada y Acumulada\n")
cat("Varianza Explicada - Media:", mean_explained_variance, "Desviación Estándar:", sd_explained_variance, "\n")
cat("Varianza Acumulada - Media:", mean_cumulative_variance, "Desviación Estándar:", sd_cumulative_variance, "\n")

# Crear tabla descriptiva con gtsummary con la media y sd de la varianza acumulada 
# y la varianza explicada
variance_summary <- variance_table %>%
  tbl_summary(
    by = NULL,  # No dividir por ningún grupo
    missing = "no",  # Excluir filas con valores faltantes
    statistic = all_continuous() ~ "{mean} ({sd})"  # Mostrar la media y la desviación estándar
  )

# Ver tabla descriptiva
View(variance_summary)



