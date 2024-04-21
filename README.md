# MSASA
**M**ultiple **S**equence **A**lignment using **S**imulated **A**nnealing 

---

Projecto desarrollado por Ing. Adrián Díaz & Dra. Gabriela Minetti dentro del contexto del Trabajo Final de Maestría en Bioinformática y Biología de Sistemas de la Universidad Nacional de Quilmes en la República Argentina.

## Acerca de
Este proyecto contiene dos componentes principales: el código fuente de una implementación de SA adaptada para encontrar MSA de proteínas y un workflow escrito en Nextflow para comparar distintas herramientas de MSA existentes en el estado del arte.

## Instalación
Una vez clonado este repositorio, se debe crear un ambiente de `Conda` para poder instalar las dependencias necesarias usando los canales `bioconda` y `conda-forge`:

- numpy
- pandas
- matplotlib
- tcoffee
- biopython
- nextflow

## Código fuente de MSASA
El código se encuentra dentro de la carpeta `src`. Hay un archivo que contiene un script que puede ejecutarse mediante CLI. 

## Workflow
El componente inicial está definido en `main.nf`. 
