# MSASA
**M**ultiple **S**equence **A**lignment using **S**imulated **A**nnealing

---

![logo](https://github.com/agdiaz/msasa_2024/blob/main/logo.png)

Projecto desarrollado por Ing. Adrián Díaz & Dra. Gabriela Minetti dentro del contexto del Trabajo Final de Maestría en Bioinformática y Biología de Sistemas de la Universidad Nacional de Quilmes en la República Argentina.

## Acerca de
Este proyecto contiene dos componentes principales: el código fuente de una implementación de SA adaptada para encontrar MSA de proteínas y un workflow escrito en Nextflow para comparar distintas herramientas de MSA existentes en el estado del arte.

### Código fuente de MSASA
El código se encuentra dentro de la carpeta `src`. Hay un archivo que contiene un script que puede ejecutarse mediante CLI.

### Workflow
El componente inicial está definido en `main.nf`.

## Instalación

Los pre-requisitos son:

- Git
- Docker
- Anaconda/Miniconda

Una vez clonado este repositorio, se debe crear un ambiente de `Conda` para poder instalar las dependencias necesarias usando los canales `bioconda` y `conda-forge`:

- python v3.12
- numpy
- pandas
- matplotlib
- tcoffee
- biopython
- nextflow

## Uso

A través de Nextflow se puede ejecutar el workflow de la siguiente manera:

```
nextflow run main.nf --input_dir ./datasets/BB11004 -with-docker -with-conda
```

La configuración de los parámetros de cada función heurística se encuentran dentro de `modules/predictors.nf`

## Sobre las heurísticas:

- Identidad: `--quality_function identity`
    - Puntaje por coincidencia: `--match_score 8.0`
    - Puntaje por no-coincidencia: `--mismatch_score 1.0`
    - Penalidad por gap: `--gap_score 2.0`
    - Temperatura inicial: `--temp 5`
    - Temperatura final: `--min_temp 0.00009`
    - Tasa de enfriamiento: `--cooling_rate 0.99`
    - Cantidad máxima de eventos sin cambiar de energía: `--max_no_changes 2000`
    - Cantidad de estados vecinos a evaluar por iteración: `--iteration_neighbors 10`
    - Cantidad máxima de cambios a efectuar en cada estado vecino: `--changes 10`
    - Experimentos: `--experiments 30`

- Coincidencias: `--quality_function coincidences`
    - Puntaje por coincidencia: `--match_score 10`
    - Puntaje por no-coincidencia: `--mismatch_score 0.75`
    - Penalidad por gap: `--gap_score 0.1`
    - Temperatura inicial: `--temp 0.90`
    - Temperatura final: `--min_temp 0.00001`
    - Tasa de enfriamiento: `--cooling_rate 0.99`
    - Cantidad máxima de eventos sin cambiar de energía: `--max_no_changes 2000`
    - Cantidad de estados vecinos a evaluar por iteración: `--iteration_neighbors 10`
    - Cantidad máxima de cambios a efectuar en cada estado vecino: `--changes 10`
    - Experimentos: `--experiments 30`

- Similitud Blosum62: ``--quality_function similarity_blosum62`
    - Puntaje por no-coincidencia: `--mismatch_score -4.0`
    - Penalidad por gap: `--gap_score -4.0`
    - Temperatura inicial: `--temp 1.0`
    - Temperatura final: `--min_temp 0.00001`
    - Tasa de enfriamiento: `--cooling_rate 0.9925`
    - Cantidad máxima de eventos sin cambiar de energía: `--max_no_changes 2000`
    - Cantidad de estados vecinos a evaluar por iteración: `--iteration_neighbors 5`
    - Cantidad máxima de cambios a efectuar en cada estado vecino: `--changes 10`
    - Experimentos: `--experiments 30`

- Similitud PAM250: ``--quality_function similarity_pam250`
    - Puntaje por no-coincidencia: `--mismatch_score -8.0`
    - Penalidad por gap: `--gap_score -8.0`
    - Temperatura inicial: `--temp 1.0`
    - Temperatura final: `--min_temp 0.00001`
    - Tasa de enfriamiento: `--cooling_rate 0.9925`
    - Cantidad máxima de eventos sin cambiar de energía: `--max_no_changes 2000`
    - Cantidad de estados vecinos a evaluar por iteración: `--iteration_neighbors 5`
    - Cantidad máxima de cambios a efectuar en cada estado vecino: `--changes 10`
    - Experimentos: `--experiments 30`

- Similitud Gonnet92: ``--quality_function similarity_gonnet`
    - Puntaje por no-coincidencia: `--mismatch_score -6.0`
    - Penalidad por gap: `--gap_score -4.0`
    - Temperatura inicial: `--temp 1.0`
    - Temperatura final: `--min_temp 0.00001`
    - Tasa de enfriamiento: `--cooling_rate 0.9925`
    - Cantidad máxima de eventos sin cambiar de energía: `--max_no_changes 2000`
    - Cantidad de estados vecinos a evaluar por iteración: `--iteration_neighbors 5`
    - Cantidad máxima de cambios a efectuar en cada estado vecino: `--changes 10`
    - Experimentos: `--experiments 30`

- Global: `--quality_function global`
    - Puntaje por coincidencia: `--match_score 10.0`
    - Puntaje por no-coincidencia: `--mismatch_score -1.0`
    - Penalidad por gap: `--gap_score -4.0`
    - Temperatura inicial: `--temp 0.80`
    - Temperatura final: `--min_temp 0.00005`
    - Tasa de enfriamiento: `--cooling_rate 0.995`
    - Cantidad máxima de eventos sin cambiar de energía: `--max_no_changes 2000`
    - Cantidad de estados vecinos a evaluar por iteración: `--iteration_neighbors 10`
    - Cantidad máxima de cambios a efectuar en cada estado vecino: `--changes 5`
    - Experimentos: `--experiments 30`

- Local: `--quality_function local`
    - Puntaje por coincidencia: `--match_score 7.0`
    - Puntaje por no-coincidencia: `--mismatch_score 1.0`
    - Penalidad por gap: `--gap_score 3.0`
    - Temperatura inicial: `--temp 0.80`
    - Temperatura final: `--min_temp 0.00009`
    - Tasa de enfriamiento: `--cooling_rate 0.995`
    - Cantidad máxima de eventos sin cambiar de energía: `--max_no_changes 2000`
    - Cantidad de estados vecinos a evaluar por iteración: `--iteration_neighbors 10`
    - Cantidad máxima de cambios a efectuar en cada estado vecino: `--changes 10`
    - Experimentos: `--experiments 30`

---

Ing. Adrián G. Díaz ([ORCID](https://orcid.org/0000-0003-0165-1318)) & Dra. Gabriela F. Minetti ([ORCID](https://orcid.org/0000-0003-1076-6766))