
# Directorio en donde está el Script de la ACTIVIDAD:
setwd("C:/Users/Karlos/OneDrive - Universidad Rey Juan Carlos/MÁSTER BIOINFORMÁTICA/TEMARIO/MASTER/asignaturas 1er Cuatri/1_ALGORTIMOS E IA/ARCHIVOS/SCRIPTS")

#Librerías necesarias para todos los algorítmos:
library(stats)   # librería para el MDS
library(ggplot2) # librería para hacer la representación gráfica
library(readr)

clases <- read.csv("classes.csv", sep = ";" , header = FALSE)
# Cambiar los nombres de las columnas
colnames(clases) <- c("Ejemplo", "clase")
clases <- clases[-1, ]

Expresión_de_genes_df <- read.csv("gene_expression.csv", sep = ";", col.names = read_lines("column_names.txt"))
sumas <- colSums(Expresión_de_genes_df) # sumo los datos por columnas
columnascero <- names(sumas[sumas==0]) # veo cuantas sumas son == 0
print(columnascero) # veo qué datos hay que excluir
Expresión_de_genes_df_final <- Expresión_de_genes_df[, !names(Expresión_de_genes_df) %in% columnascero] # reemplazo el dataset df sin esas columnas

#######################################
# Analisis de componentes principales #
#######################################


# Guardado de las variables (genes) en un data frame
data <- data.frame(sapply(Expresión_de_genes_df_final[1:497], as.numeric))

# Calculo de componentes principales con la funcion prcomp
pca.results <- prcomp(data, center=TRUE, scale=FALSE)

# Resultado de las componentes principales
pca.df <- data.frame(pca.results$x)

# Varianza (cuadrado de la desviacion tipica)
varianzas <- pca.results$sdev^2

# Total de la varianza de los datos
total.varianza <- sum(varianzas)

# Varianza explicada por cada componente principal
varianza.explicada <- varianzas/total.varianza

# Calculamos la varianza acumulada 
varianza.acumulada <- cumsum(varianza.explicada)

# Tomamos el numero de componentes principales que explican el 90% de la varianza
n.pc <- min(which(varianza.acumulada > 0.9))

# Etiquetas de los ejes del gráfico
x_label <- paste0(paste('PC1', round(varianza.explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(varianza.explicada[2] * 100, 2)), '%')

# Representación gráfica de las primeras dos componentes principales respecto a los datos
ggplot(pca.df, aes(x=PC1, y=PC2, color=clases$clase)) +
  geom_point(size=3) +
  scale_color_manual(values=c('red', 'blue', 'green', 'orange', 'purple')) +
  labs(title='PCA TIPOS DE CANCER', x= "Dimensión 1", y= "Dimensión 2", color='TIPO') +
  theme_classic() +
  theme(panel.grid.major = element_line(color="gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title = element_text(hjust = 0.5))





# Algoritmo IMAP #


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RDRToolbox")
library(RDRToolbox)



# Guardado en un dataframe de los 500 primeros genes 
data2 <- sapply(Expresión_de_genes_df_final[1:497], as.numeric) #no ponemos dataframe!

# Algoritmo

# Funcion Isomap()
#   data -> datos (matriz) sobre los que haremos reduccion de dimensionalidad
#   dim -> dimensiones de las columnas del espacio reducido
#   k -> numero de vecinos cercanos a cada punto. A mayor k mayor computacion
#   potResiduals -> devuelve la varianza explicada por las diferentes dimensiones

#   Si se ha elegido una unica dimension devuelve una matriz
#   Si se ha elegido un vector de dimensiones devolvera una matriz por cada elemento del vector

# Calculamos isomap de 1 a 10 dimensiones y con 5 vecinos (y 30?)
isomap.results = Isomap(data=data2, dims=1:4, k=15, plotResiduals=TRUE)

# Dataframe con los puntos que queremos dibujar en el plano 2D
#     (elegiriamos otro si queremos otra dimension)
isomap.df <- data.frame(isomap.results$dim2) 

# Graficamos
ggplot(isomap.df, aes(x = X1, y = X2, color = clases$clase )) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Isomap TIPOS DE CANCER", x = 'dim1', y = 'dim2', color = "TIPO") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))



# Algoritmo tSNE #
library(Rtsne)
# Seteamos la semilla para que sea replicable el algoritmo
set.seed(1234)


# Guardado en un dataframe de los 500 primeros genes 
data3 <- sapply(Expresión_de_genes_df_final[1:497], as.numeric)

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
ggplot(tsne_result, aes(x = X1, y = X2, color = clases$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método t-SNE", x = "PC1", y = "PC2", color = "TIPO") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))


# Algoritmo LLE #

# Carga de librerias
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("lle")
install.packages("lle") #NO ME PERMITE INTALAR LLE EN MI VERSIÓN DE R (4.4.2)
library(lle)



# Guardado en un dataframe de los 500 primeros genes 
data5 <- sapply(Expresión_de_genes_df_final[1:400], as.numeric)

# Algoritmo
# Funcion LLE()
#   x: matriz sobre la cual se va a reducir la dimensionalidad
#   dim: numero de dimensiones de salida
#   k: numero de vecinos cercanos. Puede aproximarse con la funcion calc_k()

#   la salida será un dataframe con dimensión = dim

# Funcion calc_k()
#   x: matriz sobre la cual se va a reducir la dimensionalidad
#   k_min, k_max: intervalo para la busqueda de la k optima
#   plotres: representará los resultados de cada k
#   paralel: utilización del calculo paralelo del ordenador
#   cpus: numero de procesadores que se quiera utilizar (max nº procesadores maquina)

calc_k <-
  function(X,m,kmin=1,kmax=20,plotres=TRUE,parallel=FALSE,cpus=2,iLLE=FALSE){
    N <- dim(X)[1]
    if(kmax>=N ) kmax <- N - 1 # más vecino que puntos no tiene sentido
    if( .Platform$OS.type=="windows" ) dev <- "nul" else dev <- "/dev/null"
    
    if( parallel==TRUE) sfInit( parallel=TRUE, cpus=cpus ) else sfInit(parallel=FALSE )
    options("warn"=-1)
    sfLibrary( lle)
    options("warn"=0)
  }

perform_calc <- function( k, X, m, iLLE=FALSE ){
  
  N <- dim(X)[1]
  
  sink( dev )
  Y <- LLE(X,m,k,2,0,iLLE=iLLE)$Y
  sink()
  
  Dx <- as.matrix(dist(X))
  Dy <- as.matrix(dist(Y))
  
  rho <- c()
  for( i in 1:N ) rho <- c( rho, cor(Dx[i,]) )
  
  return( mean(1=rho^2) )
}

rho <- invisible( sflapply( kmin:kmax, perform_calc, X, m, iLLE ))
rho <- unclass( unlist( rho ))
sfStop()

res <- data.frame( k=c(kmin:kmax), rho= rho)

if( lotres ){
  par( mar=c(5,5,4,2)+0.1 )
  plot( res$k, res$rho, type="b", xlab="k", ylab=expression(1-rho^2), main="" )
  abline(h=min(res$rho,na.rm=TRUE), col="red")
  grid()
} else cat( "best k:",head(res$k[order(res$rho)],3), "\n\n" )
return ( res )

calc_k(X=data5, m=2, kmin=1, kmax=150)


lle.results <- LLE(data5, dim = 2, k=20)


lle.df <- data.frame(lle.results)

# Graficamos
ggplot(lle.df, aes(x = X1, y = X2, color = clases$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método LLE Types of Cancer", x = "PC1", y = "PC2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

#ME VEO OBLIGADO A DESCARTAR ESTE ALGORITMO. QUERÍA UTILIZARLO PARA VER LAS AGRUPACIONES LOCALES, YA QUE LAS CONSERVA MUY BIEN.


# Algoritmo LE  #
library(Rdimtools)
install.packages("Rdimtools")

# Guardado en una matriz (dataframe) de los 500 primeros genes 
data5 <- sapply(Expresión_de_genes_df_final[1:497], as.numeric)

# Definimos los parametros k vecinos mas cercanos como el 20% de la muestra
# Funcion do.lapeig()
#   X: matriz sobre la que se reducirá la dimensionalidad
#   Ndim: dimensiones de los datos finales
#   type: forma de generar vecinos
#   preprocess: 6 tipos de preprocesamiento
#   weighted: si queremos darle pesos

#   Y: nuevo espacio con menor dimensionalidad
#   eigenvals: vector con los valores propios de las columnas de Y



le.results <- do.lapeig(data5, type=c("proportion", 0.15), weighted=FALSE)

le.df <- data.frame(le.results$Y)

# Graficamos
ggplot(le.df, aes(x = X1, y = X2, color = clases$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método LE Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))



# Algoritmo UMAP  #


# Carga de librerias
library(uwot)  #tenéis el paquete de instalación uwot en Packages


# Guardado en un dataframe de los 500 primeros genes 
data5 <- sapply(Expresión_de_genes_df_final[1:497], as.numeric)

# Se dejan los parametros a excepcion del n_neighbours

# Funcion umap()
#     x: reduccion de la dimensionalidad
#     n_neighbours: entero qeu indica el numero de vecinos cercanos
#     n_componentes: entero que determina el tamaño del espacio de salida
#     metric: define la distancia entre puntos
#     min_dist: distancia minima permitida entre puntos
#     scale: tipos de escalado
#     verbose: tiempo hasta que se complete el calculo

#     Y resultado




umap.results <- umap(data5, n_neighbors=0.2 * nrow(data5),
                     n_components = 2, min_dist = 0.1, local_connectivity=1, ret_model = TRUE, verbose=TRUE)

umap.df <- data.frame(umap.results$embedding)

# Graficamos
ggplot(umap.df, aes(x = X1, y = X2, color = clases$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método UMAP Tipos de Cáncer", x = "X1", y = "X2", color = "Tipo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))








