rm(list=ls())

path <- "C:/Users/Asus/OneDrive/bioinformatica/primer_quatri/algoritmos/actividad 3"
setwd (path)

#Carga de librerias
library(readr)
library(glmnet) # ElasticNet
library(tidyverse)
library(caret) # ML
library(rpart) # DT
library(rpart.plot) # DT plot
library(rattle) # DT plot
library(pROC) # ROC
library(PRROC) # PR-Curve
library(MASS) # LDA
library(klaR) # RDA
library(gridExtra) # juntar los gráficos

classes_df <- read_delim("classes.csv", col_names = c("Sample", "Class"), delim = ";", show_col_types = FALSE)


# Leer el archivo de nombres de columnas
column_names <- readLines("column_names.txt")
column_names

# Leer el archivo expression_values.csv
expression_genes <- read_delim("gene_expression.csv", delim = ";", show_col_types = FALSE)

# Asignar los nombres de las columnas del archivo column_names.txt al dataframe de valores de expresión
colnames(expression_genes) <- column_names

sumas <- colSums(expression_genes) # sumo los datos por columnas
columnascero <- names(sumas[sumas == 0]) # veo cuantas sumas son == 0
expression_genes <- expression_genes[,!names(expression_genes) %in% columnascero] #Sustituyo las columnas por las que no suman 0

# Existe discrepancia entre los dataframes classes_df y expression_genes

# Eliminar la fila sample_0 en classes_df 
classes_df_depurada <- classes_df[classes_df$Sample != "sample_0", ]

# Combinar los dataframes
df <- cbind(classes_df_depurada, expression_genes)

# Verificar el contenido del dataframe combinado
head(df)

genes <- names(df[3:499])


# Preparar los datos para el modelo LASSO
x <- as.matrix(df[, genes])
y <- factor(df$Class)

# Ajustar el modelo LASSO
set.seed(1995)
lasso_model <- cv.glmnet(x, y, family = "multinomial", alpha = 1)
selected_genes <- coef(lasso_model, s = "lambda.min")
selected_genes <- as.matrix(selected_genes) # Convertir a matriz densa si es necesario (esto convierte el formato disperso a un formato manejable)
selected_genes_df <- as.data.frame(selected_genes) # Convertir la matriz a data frame para facilitar su manipulación
non_zero_indices <- selected_genes_df[selected_genes_df$s1 != 0, , drop = FALSE] # Filtrar los coeficientes no cero
dim(non_zero_indices)
non_zero_indices
