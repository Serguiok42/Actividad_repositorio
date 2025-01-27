rm(list=ls())

path <- '/Users/Karlos/OneDrive - Universidad Rey Juan Carlos/MÁSTER BIOINFORMÁTICA/TEMARIO/MASTER/ALGORTIMOS E IA/ARCHIVOS/ACTIVIDAD GRUPAL'
setwd(path)

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

# Existe discrepancia entre los dataframes classes_df y expression_genes

# Eliminar la fila sample_0 en classes_df 
classes_df_depurada <- classes_df[classes_df$Sample != "sample_0", ]

# Combinar los dataframes
combined_df <- cbind(classes_df_depurada, expression_genes)

# Verificar el contenido del dataframe combinado
head(combined_df)
