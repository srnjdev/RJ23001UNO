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


### Eliminación de Gauss
```python
from RJ23001UNO import gauss_eliminacion

A = [
  [2, 1, -1],
  [-3, -1, 2],
  [-2, 1, 2]
]
b = [8, -11, -3]

solucion = gauss_eliminacion(A, b)
print(solucion)  # [2.0, 3.0, -1.0]
```

### Gauss-Jordan
```python
from RJ23001UNO import gauss_jordan

A = [
  [2, 1, -1],
  [-3, -1, 2],
  [-2, 1, 2]
]
b = [8, -11, -3]

solucion = gauss_jordan(A, b)
print(solucion)  # [2.0, 3.0, -1.0]
```

### Cramer
```python
from RJ23001UNO import cramer

A = [
    [2, 1, -1],
    [-3, -1, 2],
    [-2, 1, 2]
]
b = [8, -11, -3]

solucion = cramer(A, b)
print("Solución (Cramer):", solucion)
# Ejemplo de salida esperada: [2.0, 3.0, -1.0]
```

### Descomposición LU
```python
from RJ23001UNO import descomposicion_lu

A = [
    [2, 1, -1],
    [-3, -1, 2],
    [-2, 1, 2]
]
b = [8, -11, -3]

solucion = descomposicion_lu(A, b)
print("Solución (Descomposición LU):", solucion)
# Ejemplo de salida esperada: [2.0, 3.0, -1.0]
```

### Jacobi
```python
from RJ23001UNO import jacobi

A = [
    [4, 1],
    [2, 3]
]
b = [1, 1]

x0 = [0, 0]  # Aproximación inicial
solucion, iteraciones = jacobi(A, b, x0=x0, tol=1e-9, max_iter=1000)
print("Solución (Jacobi):", solucion)
print("Iteraciones realizadas:", iteraciones)
# La solución exacta para este sistema es [0.2, 0.2]
```

### Gauss-Seidel
```python
from RJ23001UNO import gauss_seidel

A = [
    [4, 1],
    [2, 3]
]
b = [1, 1]

x0 = [0, 0]  # Aproximación inicial
solucion, iteraciones = gauss_seidel(A, b, x0=x0, tol=1e-9, max_iter=1000)
print("Solución (Gauss-Seidel):", solucion)
print("Iteraciones realizadas:", iteraciones)
# La solución exacta para este sistema es [0.2, 0.2]
```

### Bisección
```python
import math
from RJ23001UNO import biseccion

def f(x):
    return x**2 - 2  # Queremos la raíz √2 ~= 1.4142

raiz, iteraciones = biseccion(f, 0, 2, tol=1e-9, max_iter=200)
print(raiz, iteraciones)

```

