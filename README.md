# RJ23001UNO

`RJ23001UNO` es un paquete de Python que contiene distintos métodos numéricos para
resolver sistemas de ecuaciones lineales y ecuaciones no lineales. Usando los métodos de: Eliminación de Gauss, Gauss-Jordan, Crammer, Jacobi, Gauss-Seidel, Bisección

## Instalación

Instalar en proyecto con:

```bash
pip install RJ23001UNO
```

## Métodos incluidos

- Eliminación de Gauss
- Gauss-Jordan
- Cramer
- Descomposición LU
- Jacobi
- Gauss-Seidel
- Bisección

## Tests

1. Instala `pytest`si no está instalado en tu equipo:
    ```bash
    pip install pytest
    ```

2. Desde la carpeta raíz (donde está `setup.py` y la carpeta `tests/`), ejecuta:
    ```bash
    pytest
    ```
3. Verás los resultados en la terminal. 

## Ejemplo de uso

```python
from RJ23001UNO import gauss_eliminacion

A = [[2,1,-1],
     [-3,-1,2],
     [-2,1,2]]
b = [8, -11, -3]

solucion = gauss_eliminacion(A, b)
print("Solución usando Gauss:", solucion)

