# Actividad_repositorio
Repositorio de Github,  grupo Europa 2 lote 10.

   # Variante Asignada: NC_000023.10:g.38226614G>A

## 1. Nomenclatura de la Mutación, Gen y Locus

Mediante el programa **Mutalyzer 2** y la opción “position converter”, es posible transformar la nomenclatura de la variante a **c.**. Para obtener la nomenclatura en **p**, se debe utilizar la opción “Name checker”. Como resultado obtenemos:

- **Nomenclatura de la variante en c:** NC_000023.10(OTC_v001):c.148G>A
- **Nomenclatura de la variante en p:** NC_000023.10:p.(Gly50Arg)

Respecto al número del transcrito según los códigos del NCBI, el número seleccionado es **NM_000531.5**.

El gen en el que se localiza la mutación según el número de transcrito seleccionado es el **Homo sapiens ornithine carbamoyltransferase (OTC)**, correspondiente al locus citogenético **Xp11.4** (X = cromosoma X, p = brazo corto, 11.4 = posición exacta).

Por otro lado, esta mutación se trata de una **mutación missense**, en la cual se cambia una **Glicina** por una **Arginina**.

## 2. Frecuencia Alélica Poblacional

A continuación, se determina la frecuencia alélica poblacional general mediante la base de datos **gnomAD**.

![image](https://github.com/user-attachments/assets/7c302fdd-7ebf-4f47-9083-6a62b3d39030)

# Figura 1: Resultados de gnomAD

En la **Figura 1** se puede observar la frecuencia alélica de esta variante, tanto en distintas etnias como en la población general. La **frecuencia alélica poblacional general** de esta variante es de **0.00004915**, lo que indica que es muy rara en la población general.

## 3. Conservación y Predicción de la Patogenicidad

Para predecir la patogenicidad de las variantes se utilizaron las herramientas **SIFT** y **Polyphen2**. 

- Respecto a **SIFT**, se obtuvo un **score de 0.58**, lo que sugiere que sería una mutación tolerada.
- Sin embargo, en **PolyPhen2** se obtuvo un **score de 0.817**, lo cual indica que es una mutación probablemente dañina.

En este caso, vemos que en cada programa se obtiene un resultado distinto. Aunque ambos están diseñados para predecir el impacto de las mutaciones en la función de las proteínas, **SIFT** y **PolyPhen** usan diferentes enfoques y fuentes de información, lo que puede dar lugar a discrepancias en la clasificación de las mutaciones. Una mutación podría ser vista como "tolerada" por **SIFT** y "proba

![image](https://github.com/user-attachments/assets/85c7093c-74d5-46d8-acf7-8979ce3f03a5)

# Figura 2. Análisis en Uniprot

La variante está localizada en la **región N-terminal**, que corresponde a la **cadena funcional del OTC** en la zona mitocondrial. Concretamente, se encuentra en una **hélice alfa**, una estructura secundaria crucial para la estabilidad y el funcionamiento de la proteína.

El cambio de aminoácido de **glicina** a **arginina** podría alterar esta hélice alfa y las propiedades físico-químicas de la proteína. Como se ve en la **Figura 2**, la **estabilidad predicha** del cambio es de **0.7**, lo que significa que **Uniprot** predice que es improbable que el cambio aminoacídico desestabilice la proteína.

Al ser una mutación de cambio de nucleótido, se realizó un **alineamiento de las secuencias** mediante **Clustal Omega** para la evaluación de la conservación del aminoácido afectado. Para analizar la conservación del aminoácido afectado, se comparó la **secuencia referencia de Homo sapiens** con diferentes especies como:

- **Macaca**
- **Mus musculus**
- **Gallus gallus**
- **Danio rerio**
- **S. cerevisiae**

![image](https://github.com/user-attachments/assets/a8fb066a-790d-4669-89a9-2b7596fe8837)

# Figura 3. Cambio de Aminoácido

La **Glicina (G)** marcada en azul corresponde a la **posición 50**, donde se da el cambio de nucléotido por **Arginina** en nuestra secuencia problema.

## 4. Enfermedad Hereditaria

Una vez conocida el tipo de mutación en cuestión, el gen en el que se localiza y su locus citogenético, para conocer el posible efecto fenotípico y la enfermedad asociada a la mutación, es de gran utilidad el uso de la base de datos **OMIM**. Como resultado de la búsqueda, se obtiene una mutación en el **Xp11.4** asociada a una deficiencia en la proteína **OTC**.

Esta deficiencia específica está relacionada fenotípicamente con la enfermedad conocida como **deficiencia de ornitina transcarbamilasa (OTCD, por sus siglas en inglés: Ornithine transcarbamylase deficiency)** (1). 

**OTCD** es una enfermedad ligada al cromosoma X basada en un desorden metabólico hereditario que afecta al metabolismo normal de la urea y la producción de urea, conocida como **urea génesis**. Este desorden metabólico del ciclo de la urea puede provocar diferentes afectaciones fenotípicas, como:

- **Hiperamonemia**
- **Hipotermia**
- **Alcalosis respiratoria** (1,2)
