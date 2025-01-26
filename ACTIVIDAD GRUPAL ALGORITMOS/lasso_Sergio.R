rm(list = ls())

library(dplyr)
library(readr)

classes_df <- read_delim("classes.csv", col_names = c("Sample", "Class"), delim = ";", show_col_types = FALSE)


# Leer el archivo de nombres de columnas
column_names <- readLines("column_names.txt")
column_names

# Leer el archivo expression_values.csv
expression_genes <- read_delim("gene_expression.csv", delim = ";", show_col_types = FALSE)

# Asignar los nombres de las columnas del archivo column_names.txt al dataframe de valores de expresión
colnames(expression_genes) <- column_names

# Existe discrepancia entre los dataframes classes_df y expression_genes

# Eliminar la fila sample_0 en classes_df 
classes_df_depurada <- classes_df[classes_df$Sample != "sample_0", ]

# Combinar los dataframes
combined_df <- cbind(classes_df_depurada, expression_genes)

# Verificar el contenido del dataframe combinado
head(combined_df)

library(glmnet)
library(e1071)

# 1. Remover la columna Sample y convertir Class a factor
combined_df$Class <- as.factor(combined_df$Class)
x <- as.matrix(combined_df[, !names(combined_df) %in% c("Sample", "Class")])
y <- combined_df$Class

# 2. Asegurarse que x sea numérica
x <- apply(x, 2, as.numeric)

set.seed(123)
cv_lasso <- cv.glmnet(x, y, alpha = 1, family = "multinomial")
best_lambda <- cv_lasso$lambda.min

# Seleccionar características usando Lasso
lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda, family = "multinomial")
# Obtener coeficientes para todas las clases
coef_list <- coef(lasso_model)

# Combinar los coeficientes de todas las clases
all_coefs <- do.call(cbind, lapply(coef_list, as.matrix))

# Encontrar características con coeficientes no nulos en cualquier clase
selected_features <- which(rowSums(abs(all_coefs[-1,])) != 0)

# Crear nuevo conjunto de datos con características seleccionadas
data_selected <- combined_df[, c("Class", names(combined_df)[selected_features])]




