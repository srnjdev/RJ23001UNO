# __init__.py

from .metodos import (
    gauss_eliminacion,
    gauss_jordan,
    cramer,
    descomposicion_lu,
    jacobi,
    gauss_seidel,
    biseccion
)

__all__ = [
    "gauss_eliminacion",
    "gauss_jordan",
    "cramer",
    "descomposicion_lu",
    "jacobi",
    "gauss_seidel",
    "biseccion"
]
