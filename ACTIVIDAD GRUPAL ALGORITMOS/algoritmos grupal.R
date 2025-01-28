# Instalo los paquetes necesarios y cargo las librerias 

library(dplyr)
library(tidyr)
library(ggplot2)
install.packages("data.table")
install.packages("VIM")
library(data.table)
library(VIM)

# cargamos el dateset 

library(readr)
classes <- read.csv("C:/Users/Maria/Desktop/classes.csv", sep = ";", header = FALSE, fill = TRUE)
colnames(classes) <- c("Sample", "Class")
gene_expression <- read_delim("C:/Users/Maria/Desktop/gene_expression.csv", delim = ";", col_names = FALSE)
column_names <- readLines("C:/Users/Maria/Desktop/column_names.txt")

# Asignar nombres de las columnas

colnames(gene_expression) <- column_names

# Convertir el tibble a un data.frame

gene_expression_df <- as.data.frame(gene_expression)

# Asignar nombres de las filas

rownames(gene_expression_df) <- classes$Sample
gene_expression_df$Class <- classes$Class # Añadimos la columna Class
gene_expression_df <- gene_expression_df %>%
  select(last_col(), everything())  # Movemos la columna Class al principio 

# Imputamos los datos NA

missing_values <- colSums(is.na(gene_expression_df))
print(missing_values)

# Guardamos el df depurado 

write.csv(gene_expression_df, "processed_gene_expression.csv", row.names = FALSE)


