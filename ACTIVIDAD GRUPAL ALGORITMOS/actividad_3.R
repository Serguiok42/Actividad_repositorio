rm (list = ls())
path <- "C:/Users/Asus/OneDrive/bioinformatica/primer_quatri/algoritmos/actividad 3"
setwd (path)

#Carga de las librerias
library (stats)
library (ggplot2)
library(readr) 
library(caret) 
library(plotly) 
library(vegan) 
library(Rtsne) 
library(FNN)
library(ica)
library(uwot)

#Cargo el dataset de la expresion genica con los nombres de las columnas del archivo .txt
df <- read.csv("gene_expression.csv", sep = ";", col.names = read_lines("column_names.txt"))
clases <- read.csv("classes.csv", sep = ";") #Cargo las clases 
sumas <- colSums(df) # sumo los datos por columnas
columnascero <- names(sumas[sumas == 0]) # veo cuantas sumas son == 0
df <- df[,!names(df) %in% columnascero] #Sustituyo las columnas por las que no suman 0
df$clases <- clases$classes

data <- sapply(df[,1:498], as.numeric)
anyNA(data)
