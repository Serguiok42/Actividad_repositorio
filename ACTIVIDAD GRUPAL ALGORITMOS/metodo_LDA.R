rm(list=ls())

path <- "C:/Users/Asus/OneDrive/bioinformatica/primer_quatri/algoritmos/actividad 3"
setwd (path)

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


###########------LDA------##########

# Dividir el conjunto de datos en conjuntos de entrenamiento y prueba
set.seed (123) #Determinamos la semilla de aleatoriedad
trainIndex <- createDataPartition(data_selected$Class, p = 0.8, list = FALSE) #Division de los datos 80% para entrenar, 20% para test)

train_data <- data_selected[trainIndex,] #Datos para entrenar el modelo
test_data <- data_selected[-trainIndex,] #Datos para testearlo (los contrarios al de entreno)

table(train_data$Class) #Tabla con el numero de muestras para cada tipo de clase para los datos de entreno
table(test_data$Class)#Tabla con el numero de muestras para cada tipo de clase para los datos de test 

numerical_columns <- train_data[, sapply(train_data, is.numeric)] #Seleccionamos todos los datos numericos
scaled_data <- scale(numerical_columns) #Escalado de los datos seleccionados
scaled_data <- as.data.frame(scaled_data)
training_data <-  cbind(scaled_data, Class = train_data$Class) #Creacion de tabla con los datos escalados + columna de diagnosis

numerical_columns <- test_data[, sapply(test_data, is.numeric)] #Seleccionamos todos los datos numericos
scaled_data <- scale(numerical_columns) #Escalado de los datos seleccionados
scaled_data <- scale(numerical_columns) #Escalado de los datos seleccionados
scaled_data <- as.data.frame(scaled_data)
testing_data <- cbind(scaled_data, Class = test_data$Class) #Creacion de tabla con los datos escalados + columna de diagnosis


training_data <- as.data.frame(training_data)
testing_data <- as.data.frame(testing_data)

training_data$Class <- as.factor (train_data$Class)
testing_data$Class <- as.factor (test_data$Class)

# Crear la formula sumando cada gen
length(training_data) #Para saber cuantas variables se encuentran
names <- colnames(training_data[1:100]) #La variable clases no se incluye
formula <- as.formula(paste("Class ~", paste(names, collapse = "+"))) #Se obtiene la formula del diagnosis junto a la suma de todos los parametros
formula

#### ---- LDA (lineal) ----

library(MASS) # Libreria LDA

# Ajuste del modelo LDA en el entrenamiento
lda_model <- lda(formula, data = training_data) 
lda_model$scaling # contribuciones/coeficientes


lda_pred <- predict(lda_model, newdata = training_data)
lda_pred$x

# Realizar predicciones sobre el conjunto de prueba
lda_predictions <- predict(lda_model, newdata = testing_data)
lda_predictions$x

# Obtener la predicción (predicciones de la clase)
predicted_classes <- lda_predictions$class
predicted_classes

#Comprobacion de que se ha obtenido el mismo numero de predicciones que de los datos de testing
length(predicted_classes) 
length (testing_data$Class)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testing_data$Class)
true_classes
length(true_classes)

# Crear la matriz de confusion 
confusion <- confusionMatrix(predicted_classes, true_classes)
print(confusion)

probabilities_lda <- predict(lda_model, newdata = testing_data, type = "prob") # Obtener probabilidades
print(probabilities_lda)
