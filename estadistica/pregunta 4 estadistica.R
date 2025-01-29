### 4 ###

# Cargar librerías necesarias
library(dplyr)
library(gtsummary)

# Crear terciles para cada componente del PCA
Pca_data <- as.data.frame(Pca_result$x)  # Convertir PCA a DataFrame

for (pc in colnames(Pca_data)) {
  Pca_data[[paste0(pc, "_tercil")]] <- cut(
    Pca_data[[pc]], 
    breaks = quantile(Pca_data[[pc]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = c("T1", "T2", "T3"))}

# Unir los datos PCA con las variables originales
Dataset_final <- cbind(Dataset_nuevo, Pca_data)

# Función para determinar si la variable es paramétrica o no
is_parametric <- function(x) {
  shapiro.test(x)$p.value > 0.05  # Prueba de normalidad de Shapiro-Wilk
}

# Crear tabla descriptiva con gtsummary
tabla_descriptiva <- Dataset_final %>%
  select(-starts_with("PC")) %>%  # Excluir las columnas de PCA
  tbl_summary(
    by = PC1_tercil,  # Agrupar por los terciles de PC1
    statistic = all_continuous() ~ "{mean} ({sd})",  # Estadísticas si paramétrico
    type = all_continuous() ~ "continuous2",  # Permitir cálculo de mediana e IQR
    missing = "no"
  ) %>%
  add_p(
    test = all_continuous() ~ "kruskal.test", 
    pvalue_fun = ~style_pvalue(.x, digits = 3)
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("**Tabla Descriptiva**") 

# Mostrar la tabla
tabla_descriptiva

