# tests/test_metodos.py

import pytest
from RJ23001UNO import (
    gauss_eliminacion,
    gauss_jordan,
    cramer,
    descomposicion_lu,
    jacobi,
    gauss_seidel,
    biseccion
)

def test_gauss_eliminacion():
    A = [[2,1,-1],
         [-3,-1,2],
         [-2,1,2]]
    b = [8, -11, -3]
    x_sol = [2.0, 3.0, -1.0]

    x = gauss_eliminacion(A, b)
    for i, val in enumerate(x_sol):
        assert pytest.approx(x[i], 1e-7) == val

def test_gauss_jordan():
    A = [[2,1,-1],
         [-3,-1,2],
         [-2,1,2]]
    b = [8, -11, -3]
    x_sol = [2.0, 3.0, -1.0]

    x = gauss_jordan(A, b)
    # Gauss-Jordan retorna la solución directamente en vector
    for i, val in enumerate(x_sol):
        assert pytest.approx(x[i], 1e-7) == val

def test_cramer():
    A = [[2,1,-1],
         [-3,-1,2],
         [-2,1,2]]
    b = [8, -11, -3]
    x_sol = [2.0, 3.0, -1.0]

    x = cramer(A, b)
    for i, val in enumerate(x_sol):
        assert pytest.approx(x[i], 1e-7) == val

def test_descomposicion_lu():
    A = [[2,1,-1],
         [-3,-1,2],
         [-2,1,2]]
    b = [8, -11, -3]
    x_sol = [2.0, 3.0, -1.0]

    x = descomposicion_lu(A, b)
    for i, val in enumerate(x_sol):
        assert pytest.approx(x[i], 1e-7) == val

def test_jacobi():
    A = [[4,1],
         [2,3]]
    b = [1,1]
    x_sol = [0.2, 0.2]  # Solución exacta para 4x + y = 1, 2x + 3y = 1

    x, it = jacobi(A, b, x0=[0,0], tol=1e-9, max_iter=1000)
    # Validamos que la solución sea cercana a la real
    for i, val in enumerate(x_sol):
        assert pytest.approx(x[i], 1e-7) == val

def test_gauss_seidel():
    A = [[4,1],
         [2,3]]
    b = [1,1]
    x_sol = [0.2, 0.2]

    x, it = gauss_seidel(A, b, x0=[0,0], tol=1e-9, max_iter=1000)
    # Validamos que la solución sea cercana a la real
    for i, val in enumerate(x_sol):
        assert pytest.approx(x[i], 1e-7) == val

def test_biseccion():
    import math
    def f(x):
        return x**2 - 2  # raíz en sqrt(2) ~ 1.41421356237
    raiz, it = biseccion(f, 0, 2, tol=1e-9, max_iter=200)
    
    # Validamos que la raíz esté cerca de sqrt(2)
    assert pytest.approx(raiz, 1e-7) == math.sqrt(2)
