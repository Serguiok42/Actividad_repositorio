rm(list=ls())

path <- '/Users/Karlos/OneDrive - Universidad Rey Juan Carlos/MÁSTER BIOINFORMÁTICA/TEMARIO/MASTER/asignaturas 1er Cuatri/1_ALGORTIMOS E IA/ARCHIVOS/ACTIVIDAD GRUPAL'
setwd(path)
library(stats)   # librería para el MDS
library(ggplot2)
library(readr)

classes_df <- read_delim("classes.csv", col_names = c("Sample", "Class"), delim = ";", show_col_types = FALSE)


# Leer el archivo de nombres de columnas
column_names <- readLines("column_names.txt")
column_names

# Leer el archivo expression_values.csv
expression_genes <- read_delim("gene_expression.csv", delim = ";", show_col_types = FALSE)

# Asignar los nombres de las columnas del archivo column_names.txt al dataframe de valores de expresión
colnames(expression_genes) <- column_names


#Comprobación de columnas vacías
sumas <- colSums(expression_genes) # sumo los datos por columnas
columnascero <- names(sumas[sumas==0]) # veo cuantas sumas son == 0
print(columnascero) # veo qué datos hay que excluir
Expresión_de_genes_df_final <- expression_genes[, !names(expression_genes) %in% columnascero] # reemplazo el dataset df sin esas columnas
Expresión_de_genes_df_final <- scale(Expresión_de_genes_df_final) #normalización de los datos (z-score)


# Existe discrepancia en el número muestral entre los dataframes classes_df y expression_genes

# Eliminar la fila sample_0 en classes_df 
classes_df_depurada <- classes_df[classes_df$Sample != "sample_0", ]

# Combinar los dataframes

combined_df <- cbind(classes_df_depurada, Expresión_de_genes_df_final)


# Verificar el contenido del dataframe combinado
head(combined_df)



## MÉTODOS NO SUPERVISADOS ######

       ### REDUCCIÓN DE LA DIMENSIONALIDAD #####


############## Algoritmo tSNE #############
library(Rtsne)
# Seteamos la semilla para que sea replicable el algoritmo
set.seed(1234)


# Guardado en un dataframe de los 500 primeros genes 
data3 <- sapply(combined_df[3:499], as.numeric)

# Algoritmo
# funcion Rtsne()
#   X: datos sobre los que reduciremos la dimensionalidad
#   dims: tamaño final del conjunto de datos (mejor <=3) por eficiencia
#   num_threads: hilos a utilizar (procesadores). No hace falta usarlo de momento
# 
#   Variable Y con matriz del t-SNE

tsne <- Rtsne(X=data3)
tsne_result <- data.frame(tsne$Y)

# Graficamos
ggplot(tsne_result, aes(x = X1, y = X2, color = classes_df_depurada$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método t-SNE", x = "PC1", y = "PC2", color = "TIPO") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))




            ### CLUSTERING #####

#### ---- Clustering jerarquico aglomerativo (de abajo a arriba): datos con poca cantidad de observaciones ----
library(ggdendro)
library(cluster)
library(factoextra)
library(gridExtra)


# Calcular la matriz de distancia
dist_matrix <- dist(combined_df)

# Se ejecuta el algoritmo de clusterización jerárquica aglomerativa

hclust_model_ward <- hclust(dist_matrix, method = "ward.D") # agrupa los clusters tratando de que sean lo más COMPACTOS posible minimizando la dispersión interna

colors <- rainbow(5)
k= 5
# Ward: crea clusters redondeados y homogéneos similares a los que genera k-means, PERO no funciona tan bien si los clusters tienen formas raras o tamaños muy diferentes
clust_ward <- fviz_dend(hclust_model_ward, 
                        cex = 0.5,
                        k = k,
                        palette = colors,
                        main = "Ward",
                        xlab = "Índice de Observaciones",
                        ylab = "Distancia",
                        label_row = TRUE) + 
  theme_classic()

grid.arrange(clust_ward, nrow = 2)
# ward.D: datos en los que esperas clusters compactos y homogéneos -> "Los datos deberían agruparse en clusters compactos con baja variabilidad interna."

expression_genes$cluster_ward <- as.factor(cutree(hclust_model_ward, k = 5))  #este código nos permite ver los subgrupos (númerados del 1 al 5) en el dataframe original, que no estaba etiquetado y , por tanto, poder estudiar y comparar las caracteristicas de la muestra de cada cluster en el dataframe)
