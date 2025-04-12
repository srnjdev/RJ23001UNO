"""
Archivo: metodos.py
Descripción: Implementación de métodos para resolver sistemas de ecuaciones lineales
y no lineales.
"""

import copy

def gauss_eliminacion(A, b):
    """
    Resuelve un sistema de ecuaciones lineales Ax = b usando Eliminación de Gauss.
    
    Parámetros:
    -----------
    A : list of lists (matriz de coeficientes)
    b : list (vector de términos independientes)
    
    Retorna:
    --------
    x : list (solución aproximada del sistema)
    
    Ejemplo de uso:
    ---------------
    >>> A = [[2,1,-1],
    ...      [-3,-1,2],
    ...      [-2,1,2]]
    >>> b = [8, -11, -3]
    >>> x = gauss_eliminacion(A, b)
    >>> print(x)
    [2.0, 3.0, -1.0]
    """
    A = copy.deepcopy(A)
    b = copy.deepcopy(b)
    n = len(A)
    
    # Eliminación hacia adelante
    for i in range(n):
        # Pivote
        pivote = A[i][i]
        if pivote == 0:
            raise ValueError("El pivote es 0. No se puede continuar con la eliminación de Gauss.")
        
        # Normalizar fila pivote
        for k in range(i, n):
            A[i][k] = A[i][k] / pivote
        b[i] = b[i] / pivote
        
        # Eliminar en las filas inferiores
        for j in range(i+1, n):
            factor = A[j][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            b[j] -= factor * b[i]
    
    # Sustitución hacia atrás
    x = [0]*n
    for i in range(n-1, -1, -1):
        suma = 0
        for j in range(i+1, n):
            suma += A[i][j]*x[j]
        x[i] = b[i] - suma
    return x


def gauss_jordan(A, b):
    """
    Resuelve un sistema de ecuaciones lineales Ax = b usando Gauss-Jordan.
    
    Parámetros:
    -----------
    A : list of lists (matriz de coeficientes)
    b : list (vector de términos independientes)
    
    Retorna:
    --------
    x : list (solución del sistema)
    
    Ejemplo de uso:
    ---------------
    >>> A = [[2,1,-1],
    ...      [-3,-1,2],
    ...      [-2,1,2]]
    >>> b = [8, -11, -3]
    >>> x = gauss_jordan(A, b)
    >>> print(x)
    [2.0, 3.0, -1.0]
    """
    A = copy.deepcopy(A)
    b = copy.deepcopy(b)
    n = len(A)
    
    for i in range(n):
        # Pivote
        pivote = A[i][i]
        if pivote == 0:
            raise ValueError("El pivote es 0. No se puede continuar con Gauss-Jordan.")
        
        # Normalizar la fila pivote
        for k in range(i, n):
            A[i][k] /= pivote
        b[i] /= pivote
        
        # Hacer ceros en todas las filas excepto la pivote
        for j in range(n):
            if j != i:
                factor = A[j][i]
                for k in range(i, n):
                    A[j][k] -= factor * A[i][k]
                b[j] -= factor * b[i]
                
    return b


def cramer(A, b):
    """
    Resuelve un sistema de ecuaciones lineales Ax = b usando la Regla de Cramer.
    
    Parámetros:
    -----------
    A : list of lists (matriz de coeficientes)
    b : list (vector de términos independientes)
    
    Retorna:
    --------
    x : list (solución del sistema)
    
    Ejemplo de uso:
    ---------------
    >>> A = [[2,1,-1],
    ...      [-3,-1,2],
    ...      [-2,1,2]]
    >>> b = [8, -11, -3]
    >>> x = cramer(A, b)
    >>> print(x)
    [2.0, 3.0, -1.0]
    """
    from copy import deepcopy
    
    def determinante(matriz):
        # Implementación recursiva del determinante (no es la más óptima, pero es ilustrativa)
        if len(matriz) == 1:
            return matriz[0][0]
        if len(matriz) == 2:
            return matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0]
        
        det = 0
        for c in range(len(matriz)):
            submatriz = deepcopy(matriz)
            submatriz = submatriz[1:]
            for i in range(len(submatriz)):
                submatriz[i] = submatriz[i][:c] + submatriz[i][c+1:]
            det += ((-1)**c) * matriz[0][c] * determinante(submatriz)
        return det
    
    D = determinante(A)
    if D == 0:
        raise ValueError("La matriz A es singular. No se puede usar Cramer.")
    
    n = len(A)
    x = []
    for i in range(n):
        # Reemplazar la columna i de A con b
        A_mod = deepcopy(A)
        for row in range(n):
            A_mod[row][i] = b[row]
        
        x_i = determinante(A_mod) / D
        x.append(x_i)
    
    return x


def descomposicion_lu(A, b):
    """
    Resuelve un sistema lineal Ax = b usando Descomposición LU (con método Doolittle).
    
    Parámetros:
    -----------
    A : list of lists (matriz de coeficientes)
    b : list (vector de términos independientes)
    
    Retorna:
    --------
    x : list (solución del sistema)

    Ejemplo de uso:
    ---------------
    >>> A = [[2,1,-1],
    ...      [-3,-1,2],
    ...      [-2,1,2]]
    >>> b = [8, -11, -3]
    >>> x = descomposicion_lu(A, b)
    >>> print(x)
    [2.0, 3.0, -1.0]
    """
    import copy
    n = len(A)
    A = copy.deepcopy(A)

    # Inicializar L y U como matrices identidad/cero
    L = [[0.0]*n for _ in range(n)]
    U = [[0.0]*n for _ in range(n)]
    
    # Factorización
    for i in range(n):
        # U[i][k]
        for k in range(i, n):
            suma = 0.0
            for j in range(i):
                suma += (L[i][j] * U[j][k])
            U[i][k] = A[i][k] - suma
        
        # L[k][i]
        for k in range(i, n):
            if i == k:
                L[i][i] = 1.0  # Diagonal de L en 1
            else:
                suma = 0.0
                for j in range(i):
                    suma += (L[k][j] * U[j][i])
                L[k][i] = (A[k][i] - suma) / U[i][i]
    
    # Resolver Ly = b mediante sustitución hacia adelante
    y = [0.0]*n
    for i in range(n):
        suma = 0.0
        for j in range(i):
            suma += L[i][j]*y[j]
        y[i] = (b[i] - suma)
    
    # Resolver Ux = y mediante sustitución hacia atrás
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        suma = 0.0
        for j in range(i+1, n):
            suma += U[i][j]*x[j]
        x[i] = (y[i] - suma)/U[i][i]
    
    return x


def jacobi(A, b, x0=None, tol=1e-7, max_iter=100):
    """
    Resuelve un sistema de ecuaciones lineales Ax = b usando el método iterativo de Jacobi.
    
    Parámetros:
    -----------
    A : list of lists
        Matriz de coeficientes (debe ser diagonalmente dominante o convergente).
    b : list
        Vector de términos independientes.
    x0 : list (opcional)
        Aproximación inicial para x.
    tol : float (opcional)
        Tolerancia para el criterio de parada.
    max_iter : int (opcional)
        Número máximo de iteraciones.
    
    Retorna:
    --------
    x : list (solución aproximada)
    iteraciones : int (iteraciones usadas)

    Ejemplo de uso:
    ---------------
    >>> A = [[4,1],
    ...      [2,3]]
    >>> b = [1, 1]
    >>> x_aprox, it = jacobi(A, b, x0=[0,0], tol=1e-9, max_iter=100)
    >>> print(x_aprox, it)
    """
    import math
    
    n = len(A)
    if x0 is None:
        x0 = [0]*n
    
    x = x0[:]
    for k in range(max_iter):
        x_new = [0]*n
        for i in range(n):
            suma = 0
            for j in range(n):
                if i != j:
                    suma += A[i][j]*x[j]
            x_new[i] = (b[i] - suma)/A[i][i]
        
        # Criterio de parada
        diff = max(abs(x_new[i] - x[i]) for i in range(n))
        x = x_new
        if diff < tol:
            return x, k+1
    return x, max_iter


def gauss_seidel(A, b, x0=None, tol=1e-7, max_iter=100):
    """
    Resuelve un sistema de ecuaciones lineales Ax = b usando el método iterativo de Gauss-Seidel.
    
    Parámetros:
    -----------
    A : list of lists
        Matriz de coeficientes (debe ser diagonalmente dominante o convergente).
    b : list
        Vector de términos independientes.
    x0 : list (opcional)
        Aproximación inicial para x.
    tol : float (opcional)
        Tolerancia para el criterio de parada.
    max_iter : int (opcional)
        Número máximo de iteraciones.
    
    Retorna:
    --------
    x : list (solución aproximada)
    iteraciones : int (iteraciones usadas)
    
    Ejemplo de uso:
    ---------------
    >>> A = [[4,1],
    ...      [2,3]]
    >>> b = [1, 1]
    >>> x_aprox, it = gauss_seidel(A, b, x0=[0,0], tol=1e-9, max_iter=100)
    >>> print(x_aprox, it)
    """
    import math
    
    n = len(A)
    if x0 is None:
        x0 = [0]*n
    
    x = x0[:]
    for k in range(max_iter):
        x_old = x[:]
        for i in range(n):
            suma = 0
            for j in range(n):
                if i != j:
                    suma += A[i][j]*x[j]
            x[i] = (b[i] - suma)/A[i][i]
        
        # Criterio de parada
        diff = max(abs(x[i] - x_old[i]) for i in range(n))
        if diff < tol:
            return x, k+1
    return x, max_iter


def biseccion(f, a, b, tol=1e-7, max_iter=100):
    """
    Encuentra una raíz de la ecuación f(x)=0 en el intervalo [a,b] usando el método de Bisección.
    
    Parámetros:
    -----------
    f : function
        Función f(x) de la que se busca la raíz.
    a : float
        Límite inferior del intervalo.
    b : float
        Límite superior del intervalo.
    tol : float (opcional)
        Tolerancia para el criterio de parada.
    max_iter : int (opcional)
        Número máximo de iteraciones.
    
    Retorna:
    --------
    (raiz, iteraciones)
    
    Ejemplo de uso:
    ---------------
    >>> def funcion(x):
    ...     return x**2 - 2  # Se busca la raíz de x^2 - 2 = 0
    >>> raiz, it = biseccion(funcion, 0, 2, tol=1e-9, max_iter=100)
    >>> print(raiz, it)
    """
    fa = f(a)
    fb = f(b)
    if fa*fb > 0:
        raise ValueError("La función no cambia de signo en el intervalo dado.")
    
    for i in range(max_iter):
        m = (a + b)/2
        fm = f(m)
        if abs(fm) < tol:
            return m, i+1
        
        if fa * fm < 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm
    
    return (a + b)/2, max_iter
