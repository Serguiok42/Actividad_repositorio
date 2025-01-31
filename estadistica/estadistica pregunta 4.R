rm(list = ls())

getwd()
path <- "C:/Users/Asus/OneDrive/bioinformatica/primer_quatri/estadistica y R/Actividad 3"
setwd(path)
### 4 ###
# cargar las librerias necesarias
library(stats)   # librería para el PCA
library(ggplot2) # librería para hacer la representación gráfica
library(gtsummary) # para dibujar las tablas

### preparación de los datos ###

# Cargar los datos 
Dataset <- read.csv("Dataset expresión genes.csv")

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
Pca_result <- prcomp(Data_scaled, center = TRUE, scale. = FALSE)

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

# Cargar librerías necesarias
library(dplyr)
library(gtsummary)

# Crear terciles para cada componente del PCA
Pca_data <- as.data.frame(Pca_result$x)  # Convertir PCA a DataFrame

for (pc in colnames(Pca_data)) {
  Pca_data[[paste0(pc, "_tercil")]] <- cut(
    Pca_data[[pc]], 
    breaks = quantile(Pca_data[[pc]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("T1", "T2", "T3"))}

# Unir los datos PCA con las variables originales
Dataset_final <- cbind(Dataset_nuevo, Pca_data)

# Función para determinar si la variable es paramétrica o no
is_parametric <- function(x) {
  shapiro.test(x)$p.value > 0.05  # Prueba de normalidad de Shapiro-Wilk
}

# Crear tabla descriptiva con gtsummary
tabla_descriptiva <- Dataset_final %>%
  select(PC1_tercil,starts_with("AQ")) %>%  # Excluir las columnas de PCA
  tbl_summary(
    by = PC1_tercil,  # Agrupar por los terciles de PC1
    statistic = all_continuous() ~ "{mean} ({sd})",  # Estadísticas si paramétrico
    type = all_continuous() ~ "continuous2",  # Permitir cálculo de mediana e IQR
    missing = "no",
    digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE) #Formato cientifico
  ) %>%
  add_p(
    test = all_continuous() ~ "kruskal.test", 
    pvalue_fun = ~style_pvalue(.x, digits = 3)
  ) %>%
  modify_header(label = "**Variable**") %>%
  bold_p(t = 0.05, q = FALSE) %>%    # Pvalues significativos en negrita
  modify_caption("**Componente principal 1**") 

# Mostrar la tabla
tabla_descriptiva

