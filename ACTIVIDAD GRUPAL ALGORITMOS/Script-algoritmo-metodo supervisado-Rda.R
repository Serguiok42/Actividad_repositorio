## Algoritmo e Inteligencia Artificial
## Actividad grupal 3: Análisis de un conjunto de datos de origen biológico mediante 
## técnicas de machine learning supervisadas
## Freyda Luisa Encarnación Alejandro 

## Cargar los paquetes necesarios
library(dplyr)
library(caret)
library(ggplot2)
library(randomForest)

## lectura de los datos y preparación de los datos 

expression_gene <- read.csv2("C:/Users/CESAR HINOJOSA/Downloads/gene_expression.csv", header = FALSE)
head(expression_gene)
classes <- read.csv2("C:/Users/CESAR HINOJOSA/Downloads/classes.csv", header=FALSE)
print(classes)
                         
## Con colnames le asignas los nombres a las columnas, que ya vienen del .txt
column_names <- "C:/Users/CESAR HINOJOSA/Downloads/column_names.txt"
nombres <- readLines(column_names)
head(nombres, 10)
colnames(expression_gene) <- nombres

#Esto asegurará que rownames tenga la longitud correcta.
rownames(expression_gene) <- classes$V1[1:nrow(expression_gene)]

## Rownames es para asignar nombres a las filas y le puse classes V1
## para que entendiera que son los que aparecen en la primera columna del data set classes
rownames(expression_gene) <- classes$V1

## Luego la columna 2 del data set es la que dice la clase
## Y con cbind se agrega esa columna
columna_clases <- classes$V2
expression_gene <- cbind(expression_gene, columna_clases)
View(expression_gene)

# Conversión de la columna clases a factor 
expression_gene$columna_clases <- as.factor(expression_gene$columna_clases)

## Comprobación de las dimensiones del DataFrame 
dim(expression_gene)

## Convertir la columna de clases a factor, si no está en ese formato
expression_gene$columna_clases <- as.factor(expression_gene$columna_clases)

## Para reproducubilidad
set.seed(42)

## Separar los datos en conjunto de entrenamiento y prueba
trainIndex <- createDataPartition(data$columna_clases, p = 0.8, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

table(trainData$columna_clases)
table(testData$columna_clases)

## para verificar si son caracteres
str(trainData)

dim(trainData)

## Convertir a numerico excepto la ultima columna
trainData[, -ncol(trainData)] <- lapply(trainData[, -ncol(trainData)], as.numeric)
str(trainData)

testData[, -ncol(testData)] <- lapply(testData[, -ncol(testData)], as.numeric)
str(testData)

## Escalar datos excepto la última columna
trainData[, 1:500] <- scale(trainData[, 1:500])
testData[, 1:500] <- scale(testData[, 1:500])

## convertir a dataframe los datos de prueba y de entrenamiento ##
trainDataframe <- as.data.frame(trainData)
testDataframe <- as.data.frame(testData)

View(trainDataframe)

##Elimina los guiones bajos para que no se malinterpreten como operadores
names(trainDataframe) <- gsub("-", "_", names(trainDataframe))
names(testDataframe) <- gsub("-", "_", names(testDataframe))

## Ajustar el modelo Random Forest para elegir las mejores variables
rf_model <- randomForest(columna_clases ~ ., data = trainDataframe, importance = TRUE)
print(rf_model)

# Obtener las variables o las características mas importantes
importance(rf_model)

# Almacenar la importancia de las características en un dataframe
importance_df <- data.frame(
  Feature = rownames(importance(rf_model)),
  Importance = importance(rf_model)[, 1]
)

## Ordenar las características por importancia
importance_df <- importance_df[order(-importance_df$Importance), ]

## Ordenar las 380 variables mas importantes
variables_importantes <- importance_df$Feature[1:380]

## definí una formula reducida  porque no corría con las 500 variables
## Crear la fórmula de manera dinámica
formula_reduce <- as.formula(paste("columna_clases ~", paste(variables_importantes, collapse = " + ")))
print(formula_reduce)


## verificar que todo el dataset esté correcto 
View(trainDataframe)
formula_reduce
is.factor(trainDataframe$columna_clases)
str(trainDataframe)

# verificar si hay valores faltantes
summary(trainDataframe)
dim(trainDataframe)
str(trainDataframe[, 501]) 


--RDA (regularizado)--
  
## instalar y cargar los paquete para el modelo RDA 
  
remove.packages("klaR")
install.packages("klaR")
library(klaR)
library(ggplot2)
library(caret)

## Ajustar el modelo RDA en el entrenamiento

rda_model <- rda(formula_reduce, data = trainDataframe)
rda_model

rda_pred <- predict(rda_model, newdata = trainDataframe)
View(rda_pred)

# Realizar predicciones sobre el conjunto de prueba
rda_predictions <- predict(rda_model, newdata = testDataframe)
View(rda_predictions)

# Obtener la predicción (predicciones de la clase)
predicted_classesA <- rda_predictions$class
predicted_classesA
length(predicted_classesA)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classesA <- as.factor(testDataframe$columna_clases)
true_classesA
length(true_classesA)

# Crear la matriz de confusión (Clase predicho testing vs. clase real testing)
confusionA <- confusionMatrix(predicted_classesA, true_classesA)
print(confusionA)

# Crear la matriz de confusión como un dataframe
confusion_matrix_df <- as.data.frame(as.table(confusionA))
View(confusion_matrix_df)                     

## graficar la matriz de confusión con ggplot2 
ggplot(confusion_matrix_df, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile() +  # Crear los tiles con los valores
  geom_text(aes(label = Freq), color = "white", size = 5) +  # Añadir etiquetas con los conteos
  scale_fill_gradient(low = "white", high = "blue") +  # Colores de la matriz
  theme_minimal() +  # Tema minimalista
  labs(title = "Matriz de Confusión", x = "Clase Real", y = "Clase Predicha")  # Etiquetas de los ejes


