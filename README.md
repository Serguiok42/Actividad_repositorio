# Actividad_repositorio
Repositorio de Github,  grupo Europa 2 lote 10.
Hemos decidido usar este repositorio para mostrar la realización de la tarea grupal de genética clínica y de poblaciones (explicada posteriormente), también hemos decidido seguir usando el repositorio después de la entrega para trabajar conjuntamente en los scripts de las demas actividades grupales. 

# Análisis in-sílico de la Mutación NC_000023.10:g.38226614G>A

Este proyecto tiene como objetivo el análisis detallado de una mutación genética específica localizada en el gen **OTC** (Ornithine Carbamoyltransferase), un gen ligado al cromosoma X. La variante en estudio es la **NC_000023.10:g.38226614G>A**, la cual se ha identificado como una mutación missense que cambia la glicina por la arginina en la posición 50 de la proteína OTC.

## Propósito

El propósito de este proyecto es proporcionar una comprensión más profunda de la mutación en el gen **OTC**, su implicación en la **deficiencia de ornitina transcarbamilasa** (OTCD), y los métodos de análisis genético para evaluar la patogenicidad de la variante, su frecuencia en diferentes poblaciones y su impacto potencial en la función proteica. Este análisis puede ser útil en el diagnóstico y asesoramiento genético relacionados con esta enfermedad metabólica.

## Objetivos

1. **Estudio de la mutación**: Análisis de la nomenclatura genética de la variante NC_000023.10:g.38226614G>A en sus formas c. (c.148G>A) y p. (p.Gly50Arg).
2. **Frecuencia alélica**: Determinación de la frecuencia de esta variante en diferentes poblaciones utilizando datos de gnomAD.
3. **Predicción de la patogenicidad**: Evaluación de la mutación usando herramientas bioinformáticas como **SIFT**, **PolyPhen2** y **Uniprot**.
4. **Análisis de la conservación evolutiva**: Estudio de la conservación del aminoácido afectado a través del alineamiento de secuencias con diversas especies.
5. **Estudio fenotípico**: Relación de la mutación con la enfermedad **deficiencia de ornitina transcarbamilasa (OTCD)** y los efectos clínicos asociados.
6. **Asesoramiento genético**: Proporcionar información sobre los patrones de herencia y las probabilidades de que los descendientes hereden la mutación.

## Instrucciones de Uso

Este repositorio contiene los resultados del análisis de la mutación utilizando diversas herramientas genéticas y bioinformáticas. Para replicar o modificar el análisis, se recomienda tener acceso a las siguientes herramientas y bases de datos:

1. **Mutalyzer 2** - Para la conversión de la nomenclatura genética.
2. **gnomAD** - Para obtener datos sobre la frecuencia de la variante en diferentes poblaciones.
3. **SIFT** y **PolyPhen2** - Para predecir la patogenicidad de la mutación.
4. **Uniprot** - Para el análisis de la estructura y estabilidad de la proteína afectada.
5. **Clustal Omega** - Para realizar alineamientos de secuencias y evaluar la conservación del aminoácido.

## Herramientas y Métodos

- **Mutalyzer 2**: Para convertir la nomenclatura de la variante genética de formato genómico (g.) a formato c. y p.
- **gnomAD**: Para calcular la frecuencia alélica en diversas poblaciones.
- **SIFT** y **PolyPhen2**: Para predecir el impacto funcional de la mutación.
- **Uniprot**: Para obtener información sobre la proteína OTC y cómo la mutación podría afectarla.
- **Clustal Omega**: Para realizar un análisis de conservación del aminoácido afectado por la mutación.

## Resultados y Discusión

### Nomenclatura de la mutación:
La mutación NC_000023.10:g.38226614G>A se traduce en:
- Nomenclatura c.: NC_000023.10(OTC_v001):c.148G>A
- Nomenclatura p.: NC_000023.10:p.(Gly50Arg)

### Frecuencia alélica:
La variante tiene una frecuencia alélica poblacional general muy baja, con un valor de 0.00004915, indicando que es una mutación rara en la población.

### Predicción de la Patogenicidad:
- **SIFT** predice que la mutación es tolerada, con un score de 0.58.
- **PolyPhen2** indica que la mutación podría ser dañina, con un score de 0.817.

### Conservación de la proteína:
El cambio de glicina por arginina en la posición 50 de la proteína OTC se encuentra en una región crítica de la proteína, lo que podría alterar su estabilidad y función.

### Enfermedad asociada:
La mutación está asociada con la **deficiencia de ornitina transcarbamilasa (OTCD)**, una enfermedad metabólica ligada al cromosoma X que afecta el ciclo de la urea y puede causar hiperamonemia, alcalosis respiratoria, entre otros síntomas graves.

### Asesoramiento genético:
En un escenario donde los mellizos son varones y la madre es heterocigota para la mutación, la probabilidad de que ambos estén afectados por la enfermedad es del 25%.

## Conclusiones

Este análisis genético proporciona información valiosa sobre la mutación NC_000023.10:g.38226614G>A en el gen **OTC** y su relación con la enfermedad OTCD. El uso de diversas herramientas bioinformáticas ha permitido predecir el impacto potencial de esta mutación en la función de la proteína y su posible contribución a la aparición de la enfermedad. Además, se ha discutido el patrón de herencia de esta mutación y las probabilidades de que los descendientes hereden la enfermedad.


