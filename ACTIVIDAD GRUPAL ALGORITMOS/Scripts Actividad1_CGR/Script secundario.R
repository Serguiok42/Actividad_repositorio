
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





# Algoritmo MDS #

# Funcion cmdscale()
#   d: matriz de distancias (usaremos la funcion dist)
#   k: numero que indica el tamaño final de los datos (max num de variables)
#   eig: si calculamos autovalores de las variables. Nos sirve para el calculo 
#        de la varianza explicada, es decir, para coger las columnas de mayor
#        variabilidad
#   x.ret: para devolver la matriz de distancias que calcule el algoritmo

#   points: dataframe de tamaño k que representa las nuevas coordenadas
#   eig: vector con los autovalores para elegir el numero de dimensiones


# Utilizamos la funcion dist para calcular la matriz de distancias euclideas
# matriz NxN de distancias entre todos los puntos
distances <- dist(data, method = 'euclidean')

# Utilizamos la funcon cmdscale para realizar el MSD
mds.results <- cmdscale(distances, eig=TRUE, k=3, x.ret=TRUE)

# Calculamos la varianza explicada
varianza.explicada <- mds.results$eig/sum(mds.results$eig) * 100

# Sacamos en un dataframe los puntos del mds
mds.df <- data.frame(mds.results$points)


# Grafico
ggplot(mds.df, aes(x=X1, y=X2, color = clases$clase)) +
  geom_point(size=2) + 
  scale_color_manual(values=c("red", "blue", "green", "orange", "purple")) +
  labs(title="MDS TIPOS DE CANCER", x="Dimension 1 (X1)", y="Dimension 2 (X2)", color = "Grupo") +
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

#ME VEO OBLIGADO A DESCARTAR ESTE ALGRITMO. QUERÍA UTILIZARLO PARA VER LAS AGRUPACIONES LOCALES, YA QUE LAS CONSERVA MUY BIEN.


# Algoritmo LE  #
library(Rdimtools)


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



le.results <- do.lapeig(data5, type=c("proportion", 0.50), weighted=FALSE)

le.df <- data.frame(le.results$Y)

# Graficamos
ggplot(le.df, aes(x = X1, y = X2, color = clases$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método LE Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))




#ALGRITMO MVU

# Carga de librerias
library(Rdimtools)


# Guardado en un dataframe de los 500 primeros genes 
data <- sapply(Expresión_de_genes_df_final[1:5], as.numeric)

# Definimos k a través del argumento type de la funcion do.mvu
# Funcion do.mvu()
#   x: matriz donde reduciremos la dimensionalidad
#   ndim: fimension final de los datos
#   type: cantidad de puntos vecinos
#   preprocess: preprocesamiento de datos

# Variable Y

# A mas muestras mas tiempo de ejecucion


mvu.results <- do.mvu(data, ndim=2, type=c("proportion", 0.1))

mvu.df <- data.frame(mvu.results$Y)

# Graficamos
ggplot(mvu.df, aes(x = X1, y = X2, color = clases[1:5]$clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método MVU Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
