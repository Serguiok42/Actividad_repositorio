rm(list=ls())
#Para poner el directorio base
getwd()
path <- "C:/Users/Asus/OneDrive/bioinformatica/primer_quatri/estadistica y R"
setwd(path)

#Definimos el dataset de los genes
df <-read.csv("Dataset expresión genes.csv")

#Carga de librerias
library(moments) # asmietrias y curtosis
library(ggplot2) # graficos
library(dplyr) # programación
library(gridExtra) # figuras
library(car) # levene test
library(gtsummary) # tabla resumen
library (gt) #Tabla gt
library (tidyverse)

#--Ejercicio 1----
#Dataset con solo los genes
df_genes <- df %>% select(starts_with("AQ_"))

#Hacemos una matriz vacia para los p-values
length(df_genes) #Para saber de cuantas filas es la matriz
pvalues1 <- matrix(NA, nrow = 46, ncol = 1, dimnames = list(c(NULL), c("valor_p")))

#definimos la lista de variables con los genes
variables <- colnames (df_genes)

#Bucle para cada uno de los genes
for (i in 1:length(variables)) {
  shapiro_result <- shapiro.test(df_genes[[variables[i]]])
  pvalues1[i,1] <- shapiro_result$p.value
}

#Definimos las filas de las columnas como el nombre de los genes 
rownames(pvalues1) <- variables


pvalues1 <- as.data.frame (pvalues1)

pvalues1 <- pvalues1 %>%
  tibble::rownames_to_column("Variable")

pvalues1 <- pvalues1 %>%
  mutate(
    Prueba = "Shapiro-Wilk",
    Interpretacion = ifelse(valor_p < 0.05, "Variable no paramétrica", "Variable paramétrica"),
  )
# Crear una tabla con gt
tabla1 <- pvalues1 %>%
  gt() %>%
  tab_header(
    title = "Estudio de normalidad de los genes",
  ) %>%
  cols_label(
    Variable = "Variable",
    valor_p = "Valor p",
    Prueba = "Prueba utilizada",
    Interpretacion = "Interpretacion",
    
  ) %>%
  fmt_scientific(
    columns = valor_p,
    decimals = 3
  ) %>%
  tab_source_note(
    source_note = "Nota: Los valores p se obtuvieron usando la prueba Shapiro-Wilk."
  )

# Mostrar la tabla
tabla1

#Guardar la tabla en formato png para el word
gtsave(tabla1, filename = "tabla1.png")

####-----Variables bioquimicas, sociodemograficas y sintomaticas-----
df_variables <- df[37:55]

#Hacemos una matriz vacia para los p-values
length(df_variables) #Para saber de cuantas filas es la matriz
pvalues_bioquimicos <- matrix(NA, nrow = 19, ncol = 1, dimnames = list(c(NULL), c("valor_p")))

#definimos la lista de variables con los genes
variables2 <- colnames (df_variables)

#Bucle para cada uno de los genes
for (i in 1:length(variables2)) {
  shapiro_result <- shapiro.test(df_variables[[variables2[i]]])
  pvalues_bioquimicos[i,1] <- shapiro_result$p.value
}

#Definimos las filas de las columnas como el nombre de los genes 
rownames(pvalues_bioquimicos) <- variables2


pvalues_bioquimicos <- as.data.frame (pvalues_bioquimicos)

pvalues_bioquimicos <- pvalues_bioquimicos %>%
  tibble::rownames_to_column("Variable")

pvalues_bioquimicos <- pvalues_bioquimicos %>%
  mutate(
    Prueba = "Shapiro-Wilk",
    Interpretacion = ifelse(valor_p < 0.05, "Variable no paramétrica", "Variable paramétrica"),
  )
# Crear una tabla con gt
tabla_bioquimica <- pvalues_bioquimicos %>%
  gt() %>%
  tab_header(
    title = "Estudio de normalidad de las variables bioquimicas",
  ) %>%
  cols_label(
    Variable = "Variable",
    valor_p = "Valor p",
    Prueba = "Prueba utilizada",
    Interpretacion = "Interpretacion",
  ) %>%
  fmt_scientific(
    columns = valor_p,
    decimals = 3
  ) %>%
  tab_source_note(
    source_note = "Nota: Los valores p se obtuvieron usando la prueba Shapiro-Wilk."
  )
tabla_bioquimica

#-----Ejercicio 2: contraste de hipotesis entre tel tratamiento y los tipos de tumor

df$trat <- as.factor (df$trat)
levels (df$trat) #tratA y tratB
df$tumor <- as.factor (df$tumor)
levels(df$tumor) #CCR, CM y CP

table (df$trat, df$tumor)
prop.table (table(df$trat, df$tumor))*100
# H0: CCR[trat] = CM[trat] = CP[trat] -> p≥0.05
# H1: CCR[trat] ≠/= CM[trat] ≠/= CP[trat] -> p<0.05


#---Normalidad en funcion de grupo tratA y tratB----

# H0: Hay normalidad en los dos grupos TratA o tratB -> p≥0.05
# H1: No existe normalidad en los genes en tratA o tratB  -> p<0.05

#Tabla en funcion del tratamiento
table (df$trat)

#Definimos los dataframes para tratamiento A y B
df_tumor_tratA <- df %>% filter (trat=="tratA")  %>% select (tumor, starts_with("AQ_"))
df_tumor_tratB <- df %>% filter (trat == "tratB") %>% select (tumor, starts_with("AQ_"))


#Hacemos una matriz vacia para los dos tratamientos 
pvalues_tumor_trat <- matrix(NA, nrow = 46, ncol = 2, dimnames = list(c(NULL), c("tratA", "tratB")))

for (i in 1:length(variables)) {
  shapiro_result <- shapiro.test(df_tumor_tratA[[variables[i]]])
  pvalues_tumor_trat[i,1] <- shapiro_result$p.value
  
  shapiro_result <- shapiro.test(df_tumor_tratB[[variables[i]]])
  pvalues_tumor_trat[i,2] <- shapiro_result$p.value
}

rownames(pvalues_tumor_trat) <- variables
format(pvalues_tumor_trat, scientific = FALSE)



#p > 0,05 en todos los casos (podria poner una tabla gt para la normalidad)
#No se estudia la normalidad en los subgrupos de cancer porque la n es muy pequeña

# k > 2; Datos no parametricos --> kruskal wallis

#---Tabla descriptiva ejercicio 2----

#Tabla descriptiva en funcion del tratamiento y el tumor para los genes

#Para saber cuantas k tiene la variable tratamiento y tumor



#Filtrar para obtener la tabla con solo los genes, el tumor y el tratamiento
df_tabla <- df %>% select (trat, tumor, starts_with("AQ_"))

#Tabla con los genes segun tratamiento y tumor
df_tabla <- df %>% select (trat, tumor, starts_with("AQ_"))

#Tabla con los genes segun tratamiento y tumor
tabla2 <- df_tabla %>% 
  tbl_strata (strata = trat,
              .tbl_fun = ~.x %>%
                tbl_summary (by = tumor, 
                             statistic = all_continuous () ~ "{p50} ({p25} - {p75})",
                             digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
                add_p ( test = list (all_continuous () ~"kruskal.test",
                                     all_categorical () ~ "chisq.test"),
                        pvalue_fun = ~ style_pvalue (.x, digits = 3)
                ) %>%
                modify_header(label = "**Gen**") %>%
                bold_p(t = 0.05, q = FALSE) %>%    # Pvalues significativos en negrita
                bold_labels() %>%
                modify_caption("**Tabla estadistica de los genes respecto a tratamiento y tumor**")
  )
tabla2gt <- as_gt(tabla2)

tabla2gt <- tabla2gt %>%
  tab_source_note(
    source_note = "Figura 2. Analisis estadisticos de los genes segun tratamiento y tipo de tumor. Los valores estadisticos se encuentran representados por la mediana y los rangos intercuartilicos (p25-p75). 
    Uso de Kruskal-Wallis como prueba estadistica. Los p-values < 0.05 se encuentran destacados en negrita.
    "
  )
#Guardar la tabla en formato png para el word
gtsave(tabla2gt, filename = "tabla2gt.png")

gtsave(tabla2gt, filename = "tabla2gt.html")

#----Ejercicio 3----
#Contraste de hipotesis entre la expresion genica y la edad 
# H0: Expresion[Edad > p50] = Expresion[Edad < p50] -> p≥0.05
# H1: Expresion[Edad > p50] ≠/= Expresion[Edad < p50] -> p<0.0
# k = 2; Variables independientes; n > 30; Homogeneidad de Varianzas

#----Separar en dos grupos segun > o < de la mediana de edad----

#Definimos cual es la mediana de la edad

mediana_edad <- quantile(df$edad, probs = 0.5, na.rm = TRUE)

#Creacion de una nueva columna con dos categorias: es > o < al percentil 50 (mediana) de la edad
df$categorias_edad <- cut (df$edad,
                           breaks = c(-Inf, mediana_edad, Inf),
                           labels = c ("Edad < p50", "Edad >= p50"),
                           right = FALSE)
df$categorias_edad <- as.factor (df$categorias_edad)


#-----Estudio de la normalidad-----
#Tabla en funcion del tratamiento
table (df$categorias_edad)

#Definimos los dataframes para tratamiento A y B
df_edad_sup <- df %>% filter (categorias_edad =="Edad < p50")  %>% select (starts_with("AQ_"))
df_edad_inf <- df %>% filter (categorias_edad == "Edad >= p50") %>% select (starts_with("AQ_"))


#Hacemos una matriz vacia para los dos tratamientos 
pvalues <- matrix(NA, nrow = 46, ncol = 2, dimnames = list(c(NULL), c("Edad < p50", "Edad >= p50")))

for (i in 1:length(variables)) {
  shapiro_result <- shapiro.test(df_edad_sup[[variables[i]]])
  pvalues[i,1] <- shapiro_result$p.value
  
  shapiro_result <- shapiro.test(df_edad_inf[[variables[i]]])
  pvalues[i,2] <- shapiro_result$p.value
}

rownames(pvalues) <- variables
format(pvalues, scientific = FALSE)

#Son todos no parametricos, hacer una tabla supongo
#------Estudio de la homogeneidad de varianzas con test de levene---

# H0: Hay homogeneidad de varianzas entre ambos grupos  -> p≥0.05
# H1: No hay homogeneidad de varianzas entre ambos grupos -> p<0.05
pvalues_edades <- matrix(NA, nrow = 46, ncol = 1, dimnames = list(c(NULL), c("valor_p")))
pvalues_edades
length (variables)
for (i in 1:length(variables)) {
  levene <- leveneTest(df[[variables[i]]], df$categorias_edad)
  pvalues_edades[i,1] <- levene$`Pr(>F)`[1]
}

rownames (pvalues_edades) <- variables
format (pvalues_edades, scientific = FALSE)


pvalues_edades <- as.data.frame (pvalues_edades)

pvalues_edades <- pvalues_edades %>%
  tibble::rownames_to_column("Variable")

pvalues_edades <- pvalues_edades %>%
  mutate(
    Prueba = "Levene Test",
    Interpretacion = ifelse(valor_p < 0.05, "NO homogeneidad de varianzas", "Homogeneidad de varianzas"),
  )
pvalues_edades

#El unico gen con NO homogeneidad de varianzas es AQ_IL6, por lo que solo en ese hay que usar Welch
#Tabla con los genes segun edad (mas grande que la mediana o menor)

tabla_categoria_edad <- df %>% 
  select (categorias_edad,starts_with("AQ_")) %>%
  tbl_summary (by = categorias_edad,
               statistic = all_continuous () ~ "{p50} ({p25} - {p75})",
               digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
  add_p ( test = list (all_continuous() ~ "t.test",
                       all_categorical () ~ "chisq.test"),
          pvalue_fun = ~ style_pvalue (.x, digits = 3)
  ) %>%
  modify_header(label = "**Gen**") %>%
  bold_p(t = 0.05, q = FALSE) %>%
  bold_labels() %>%
  modify_caption("**Tabla estadistica de los genes respecto a las categorias de edad**") 

tabla_categoria_edad <- as_gt(tabla_categoria_edad)

tabla_categoria_edad <- tabla_categoria_edad %>%
  tab_source_note(
    source_note = "Figura 4. Analisis estadisticos de los genes segun la categoria de edad. Los valores estadisticos se encuentran representados por la mediana y los rangos intercuartilicos (p25-p75). Se utilizo Welch Two Sample t-test para el analisis estadistico. P-values < 0.05 se encuentran destacados en negrita."
  )
tabla_categoria_edad
#Guardar la tabla en formato png para el word
gtsave(tabla_categoria_edad, filename = "tablaedad.png")
