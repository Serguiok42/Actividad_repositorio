
### Master Universitario en Bioinformática
### Estadistica y R para la Salud 
### Actividad grupal 3: Análisis de un caso práctico en R 


# cargar las librerias necesarias
library(stats)   # librería para el PCA
library(ggplot2) # librería para hacer la representación gráfica
library(gtsummary) # para dibujar las tablas

### preparación de los datos ###

# Cargar los datos 
Dataset <- read.csv("C:/Users/CESAR HINOJOSA/Downloads/Dataset expresión genes.csv")

# para ver los nombres de las columnas
colnames(Dataset)

# Sacar una lista de las columnas que contienen la expresión de los genes y la columnna id
columnas_genes <- c("id", "AQ_ADIPOQ", "AQ_ALOX5", "AQ_ARG1", "AQ_BMP2", "AQ_CCL2", "AQ_CCL5", "AQ_CCR5", "AQ_CD274",      
                    "AQ_CD36", "AQ_CHKA", "AQ_CPT1A", "AQ_CSF2", "AQ_CXCR1", "AQ_FASN", "AQ_FOXO3", "AQ_FOXP3", "AQ_G6PD",        
                    "AQ_IL10", "AQ_IL1B", "AQ_IL6", "AQ_IRS1", "AQ_JAK1", "AQ_JAK3", "AQ_LDHA", "AQ_LIF", "AQ_MAPK1", "AQ_NFE2L2",      
                    "AQ_NFKB1", "AQ_NLRP3", "AQ_NOS2", "AQ_NOX5", "AQ_PDCD1", "AQ_PPARG", "AQ_PTAFR",       
                    "AQ_PTGS2", "AQ_SLC2A4", "AQ_SOD1", "AQ_SREBF1", "AQ_STAT3",       
                    "AQ_TGFB1", "AQ_TLR3", "AQ_TLR4", "AQ_TNF", "AQ_GPD2", "AQ_GPX1", "AQ_IFNG")


# Crear un nuevo dataset con solo esas columnas
Dataset_nuevo <- Dataset[, columnas_genes]

# ver dataset nuevo
print(Dataset_nuevo)

# Escalar los datos 
Data_scaled <- scale(Dataset_nuevo)

# Verificar si las columnas han sido correctamente estandarizadas
summary(Data_scaled)

# Aplicar PCA
Pca_result <- prcomp(Data_scaled, center = TRUE, scale. = TRUE)

# Ver el resumen del PCA
summary(Pca_result)


## Tabla de la contribucion de cada cada gen por componentes

# Matriz de cargas para ver la contribución de cada gen a los componentes
Pca_result$rotation

# Convertir la matriz de cargas en un dataframe
tabla_cargas <- as.data.frame(Pca_result$rotation)

# elimina la primera fila "id" a todas las columnas del dataframe tabla_cargas
tabla_cargas_1 <- tabla_cargas[-1, ]
print(tabla_cargas_1)

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



