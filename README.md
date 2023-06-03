VASP Post-Processing Scripts
Este repositorio contiene una serie de scripts escritos en Python3 diseñados para analizar y post-procesar los resultados de VASP. Estos scripts son útiles para el modelado termodinámico utilizando la teoría del \textit{CHE}, por ejemplo, en el estudio de la reacción de reducción de oxígeno (ORR).

Scripts disponibles
Cada script en este repositorio se nombra de acuerdo con el archivo de salida VASP que puede leer y analizar, o el tipo de post-procesado que puede realizar. Aquí está la lista de scripts disponibles y sus funciones:

BADER.py: Este script puede ser utilizado para analizar la carga de Bader a partir de los resultados de VASP.

CHGCAR.py: Este script procesa el archivo CHGCAR de VASP, el cual contiene la densidad de carga.

DOSCAR.py: Este script puede ser utilizado para analizar el archivo DOSCAR de VASP, el cual contiene la densidad de estados.

Data.py: Este script probablemente se utiliza para gestionar y manipular los datos obtenidos de las simulaciones de VASP.

EIGENVAL.py: Este script analiza el archivo EIGENVAL de VASP, el cual contiene los autovalores de los estados electrónicos.

INCAR.py: Este script puede manipular o generar archivos INCAR, que contienen los parámetros de entrada para las simulaciones de VASP.

Logs.py: Este script puede ser utilizado para procesar y analizar archivos de registro generados por las simulaciones de VASP.

ORR.py: Este script parece estar especializado en el análisis y post-procesamiento de la reacción de reducción de oxígeno (ORR) utilizando los resultados de VASP.

OSZICAR.py: Este script analiza el archivo OSZICAR de VASP, el cual contiene la evolución de la autoconsistencia electrónica y la geometría durante la simulación.

OUTCAR.py: Este script procesa el archivo OUTCAR de VASP, que contiene una variedad de información sobre la simulación, incluyendo energías totales y fuerzas.

PARCHG.py: Este script maneja el archivo PARCHG de VASP, que contiene la densidad de carga proyectada.

PDB.py: Este script puede manipular y/o generar archivos en formato PDB, comúnmente utilizados para la visualización de estructuras moleculares.

PDB_Search.py: Este script probablemente busca en una base de datos de estructuras en formato PDB.

POSCAR.py: Este script maneja el archivo POSCAR de VASP, que contiene la información de la estructura atómica.

POTCAR.py: Este script puede manipular o generar archivos POTCAR, que contienen la información de los pseudopotenciales utilizados en las simulaciones de VASP.

PROCAR.py: Este script analiza el archivo PROCAR de VASP, el cual contiene las proyecciones de los orbitales atómicos en los estados electrónicos.

QMLEARNING.py: Este script parece estar relacionado con la implementación de técnicas de aprendizaje automático para analizar los resultados de las simulaciones de VASP.

SDF.py: Este script puede manipular y/o generar archivos en formato SDF, comúnmente utilizados para almacenar información de estructuras moleculares.

Set.py: Este script probablemente maneja conjuntos de datos de las simulaciones de VASP.

System.py: Este script puede estar relacionado con la configuración del sistema para las simulaciones de VASP.

analisys.py: Este script parece estar diseñado para realizar análisis generales de los resultados de las simulaciones de VASP.

cookbook.py: Este script puede contener recetas o flujos de trabajo comunes para trabajar con los resultados de VASP.

data2XML.py: Este script probablemente convierte los datos de las simulaciones de VASP a formato XML.

dataset_maker.py: Este script puede estar diseñado para generar conjuntos de datos a partir de los resultados de las simulaciones de VASP para su posterior análisis.

functional_base.py: Este script puede ser un módulo de base para definir las funcionalidades de otros scripts.

visual.py: Este script probablemente proporciona funcionalidades para la visualización de los datos de las simulaciones de VASP.


