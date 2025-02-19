rm(list = ls())

getwd()
path <- "C:/Users/Asus/OneDrive/bioinformatica/primer_quatri/estadistica y R/Actividad 3"
setwd(path)

### ACTIVIDAD COMPLETA ###
# cargar las librerias necesarias
library(stats)   # librería para el PCA
library(ggplot2) # librería para hacer la representación gráfica
library(gtsummary) # para dibujar las tablas
library(nnet)
library(broom)
library(gt)
library(dplyr)

### preparación de los datos ###

# Cargar los datos 
Dataset <- read.csv("Dataset expresión genes.csv")

# para ver los nombres de las columnas
colnames(Dataset)

# Sacar una lista de las columnas que contienen la expresión de los genes 
columnas_genes <- c("AQ_ADIPOQ", "AQ_ALOX5", "AQ_ARG1", "AQ_BMP2", "AQ_CCL2", "AQ_CCL5", "AQ_CCR5", "AQ_CD274",      
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
# Convertir a dataframe
Data_scaled_df <- as.data.frame(Data_scaled)

# Verificar si las columnas han sido correctamente estandarizadas
summary(Data_scaled)

# Aplicar PCA

Pca_result <- prcomp(Data_scaled, center = TRUE, scale. = FALSE)


# Ver el resumen del PCA
summary(Pca_result)


## Tabla de carga de la contribucion de cada gen por componentes

# Matriz de cargas (contribución de cada gen a los componentes)
Pca_result$rotation

# Convertir la matriz de cargas en un data frame
tabla_cargas <- as.data.frame(Pca_result$rotation)

# elimina la primera fila a todas las columnas del dataframe tabla_cargas
tabla_cargas_1 <- tabla_cargas[-1, ]
print(tabla_cargas_1)

### PARTE 3 ###
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



###########FUNCIONES EXTRA PARA LA GENERACIÓN DE GRÁFICOS y OBTENCIÓN DE INFORMACIÓN RELEVANTE########

library(factoextra)

#VisualizaCIión de los eigenvalues de los componentes principales
eigenvalues <- fviz_eig(Pca_result, addlabels = TRUE, ylim = c(0, 50))
eigenvalues # los eigenvalues (valores propios) miden el nivel de estiramiento o encogimiento de sus vectores asociados tras la transformació de los datos (PCA). 
            # Por tanto, a mayor magnitud del valor propio (%)  su  dirección asociada (vector propio) mayor  información importante contiene.En nuestro contexto, a mayor eigenvalue mayor información de los datos (varianza explicada) contiene la dimension del PCA  . 



#Extraer los resultados de las variables
var <- get_pca_var(Pca_result)

#Para visualizar como las variables se asocian con las dimensiones del PCA realizamos un gráfico de correlación:

correlación_var_PC1y2 <-  fviz_pca_var(Pca_result, col.var = "black")

correlación_var_PC1y3 <- fviz_pca_var(Pca_result, col.var = "black", axes = c(1, 3))
#```{r}
#Resultado: prácticamente todos los genes se correlación positivamente con la dimensión 1. Para la dimensión 2 más o menos la mitad de las variables se correlación positivamente y la mitad negativamente (algo similar ocurre para la dimensión 3 (PC3)).
# la longitud de las flechas indican el grad de correlaciión de la variable con la dimensión. En este caso a tener tantas variables es difíl distinguirlas, ya que los títulos se solapan, pero esta gráfica nos permite ver que variablles (genes) están más relacionados unos con otros. 
# ```

# También podemos Visualizar el gráfico de correlación de las variables por el coseno (cos2). A mayor valor del coseno mayor contribución con la dimension 
correlación_var_COs2 <- fviz_pca_var(Pca_result, col.var = "cos2", 
                                     gradient.cols = c("blue", "yellow", "red"), 
                                     repel = TRUE) # evitar el overlapping
#```{r}
#Podemos observar en este gráfico como AQ_Nox5 y AQ_ADIPOQ tienen una contribución negativa baja  para las dimensiones 1 y 2.
#Además de ver que tienen poca calidad ( no contribuyen aenas a estas dimensiones), podemos ver que estos genes correlacionan entre sí y no con el resto.
# ```


#Para poder mejorar la visualización del gráfico anterior y poder agrupar las variables genes por clusters:
# Realizar k-means clustering
kmeans <- kmeans(Data_scaled_df, centers=5)
grupo <- as.factor(kmeans$cluster)
# Asignar colores a los clusters
colores <- c("grey", "yellow", "blue", "pink", "red")
colores_clusters <- colores[grupo] 
# Asegurarse de que la longitud de colores_clusters coincida con el número de variables
colores_clusters <- colores_clusters[1:nrow(Pca_result$var$coord)]
# Visualizar PCA con colores según los clusters
fviz_pca_var(Pca_result, col.var = colores_clusters, #ERROR EN COL.VAR. No soy capaz de solucionarlo
             gradient.cols = c("grey", "yellow", "blue", "pink", "red"), 
             legend.title = "Cluster", 
             repel = TRUE) # evitar el overlapping


#Visualizar la importancia de cada variable por dimensiones
fviz_contrib(Pca_result, choice= "var", axes= 1, top = 500) 
fviz_contrib(Pca_result, choice= "var", axes= 2, top = 500)
fviz_contrib(Pca_result, choice= "var", axes= 45, top = 500)

#En el PCA1 podemos observar que hay 27 variables (de AQ_JAK1 a AQ_PPAR, en el gráfico) de mayor contribución, >2%  (por encima de la linea que corta el eje de ordenadas)
 #otros 14 genes  contribuyen aunque en menor medida, entre el 2 y 1% (de AQ_FOXP3 a AQ_NOX2, en el gráfico). El resto de genes tiene una contribución al PCA1 menor, = o < al 1% (para cada uno de ellos)
# Cargar librerías necesarias

#En el PCA2 realmente hay 6 variables son las que producen mayor contribución, a partir de aquí el peso es mucho menor (parece lógico,  si pensamos que el PC2 representa sólo el 6.5% de la varianza, que sean pocas variables que explican la mayor parte de la contribución)
#Cuanto más dimensiones evaluemos esta tendencia será más acusada, ues cada vez explica menos varianza. Por ejemplo, en el pcA45, tan sólo un gen (AQ_PTAFR) explica más del 20% de la contribución.


#vALORAMOS LOS PATRONES DE LOS PACIENTES

fviz_pca_ind(Pca_result, col.ind = "cos2", 
             gradient.cols = c("blue", "yellow", "red"), 
             repel = TRUE)

fviz_pca_ind(Pca_result, pointsize = "cos2", 
             pointshape = 21, 
             fill = "yellow", 
             repel = TRUE)

#Se pueden distinguir distintos patrónes de pacientes en funcíon de la expresión de los genes.
fviz_pca_biplot(Pca_result, repel = TRUE,
                col.var = "#2E9FDF", # Color de las variables
                col.ind = "#696969"  # Color de los individuos
)
#La mayor parte de los genes se relacionan con pacientes que se correlacionan negativamente con la dimension 1. Se puede observar un clúster de pacientes (correlación positiva con Dim1) que no se asocian prácticamente a genes, salvo a los genes AQ_NOX5 y AQ_ADIPOQ para algunos de ellos


### PARTE 4 ###
###### GENERACIÓN DE TABLAS DESCRIPTIVAS######

#Creación de terciles para los componentes del PCA#

library(dplyr)
library(gtsummary)


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
tabla_descriptiva_PC1 <- Dataset_final %>%
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

# Crear tabla descriptiva para los terciles de PC2
tabla_descriptiva_PC2 <- Dataset_final %>%
  select(PC2_tercil, starts_with("AQ")) %>%  # Incluir las columnas que comienzan con "AQ"
  tbl_summary(
    by = PC2_tercil,  # Agrupar por los terciles de PC2
    statistic = all_continuous() ~ "{mean} ({sd})",  # Estadísticas si paramétrico
    type = all_continuous() ~ "continuous2",  # Permitir cálculo de mediana e IQR
    missing = "no",
    digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE) # Formato científico
  ) %>%
  add_p(
    test = all_continuous() ~ "kruskal.test", 
    pvalue_fun = ~style_pvalue(.x, digits = 3)
  ) %>%
  modify_header(label = "**Variable**") %>%
  bold_p(t = 0.05, q = FALSE) %>%  # P-values significativos en negrita
  modify_caption("**Componente principal 2**")

# Combinar las tablas descriptivas
tabla_descriptiva_combinada <- tbl_merge(
  tbls = list(tabla_descriptiva_PC1, tabla_descriptiva_PC2),
  tab_spanner = c("**Componente principal 1**", "**Componente principal 2**")
)

# Modificar el título superior de la tabla combinada
tabla_combinada <- tabla_descriptiva_combinada %>%
  modify_caption("**Componentes principales**")

tabla_combinada_gt <- as_gt(tabla_combinada)

tabla_combinada_gt <- tabla_combinada_gt %>%
  tab_source_note(
    source_note = "Tabla estadistica. Analisis estadisticos de los diferentes genes por terciles para los componentes principales 1 y 2.   Los valores estadisticos se encuentran representados por la media y la desviacion estandard y los rangos intercuartilicos (p25-p75). 
    Uso de Kruskal-Wallis como prueba estadistica. Los p-values < 0.05 se encuentran destacados en negrita")
    
tabla_combinada_gt




###PARTE 5########


# Necesitamos que la variable metástasis sea binomial, (si = clase 1/no = clase 0) asi que creo una nueva columna usando los datos de la variable extensión.
Dataset_Regresion <- Dataset %>%
  mutate(metastasis = ifelse(extension == "metastasico", 1, 0))


# Uno los datos con las variables clínicas y la nueva variable metastasis
data_final_Regresion <- cbind(Dataset_Regresion[, c("metastasis", "edad", "sexo", "tumor")], Pca_data)

# convierto la variable metastasis y terciles a factor
data_final_Regresion$metastasis <- as.factor(data_final_Regresion$metastasis)
data_final_Regresion$PC1_tercil <- as.factor(data_final_Regresion$PC1_tercil)
data_final_Regresion$PC2_tercil <- as.factor(data_final_Regresion$PC2_tercil)
data_final_Regresion$PC3_tercil <- as.factor(data_final_Regresion$PC3_tercil)

# Aplico el modelo de regresión logística
model <- glm(metastasis ~ PC1_tercil + PC2_tercil + PC3_tercil + edad + sexo + tumor, 
             data = data_final_Regresion, 
             family = "binomial")

# Creo la tabla con OR, IC 95% y p-value
table_results_Regresion <- tidy(model, exponentiate = TRUE) %>% 
  mutate(`IC 95%` = paste0(round(exp(confint.default(model)[, 1]), 2), " - ", round(exp(confint.default(model)[, 2]), 2))) %>%
  select(term, estimate, `IC 95%`, p.value)

# Formateo la tabla usando gt
table_results_Regresion %>%
  gt() %>%
  tab_header(title = md("**Tabla de Regresión Logística**"),
             subtitle = "Metástasis (Sí/No) por terciles de PCA") %>%
  fmt_number(columns = c(estimate, p.value), decimals = 3) %>%
  cols_label(term = "Variable", estimate = "OR", `IC 95%` = "IC 95%", p.value = "P-valor") %>%
  opt_table_outline() 
