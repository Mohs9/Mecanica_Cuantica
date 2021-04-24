'''
Autor:  Rubí E. Ramírez Milián
Compliador:
Para compilar: python QM_P2.py
Fecha:  Sat Apr 23 17:38:09 CST 2021
Librerías: math, sympy.combinatorics.permutations, itertools
Entradas, Salidas, Resumen: Se hará un programa para calcular conmutadores. Matrices
'''

from math import sqrt
from sympy.combinatorics.permutations import Permutation
from itertools import permutations

paulix=[[0,1],[1,0]]
pauliy=[[0,-1j],[1j,0]]
pauliz=[[1,0], [0,-1]]
identidad=[[1,0],[0,1]]


#FUNCION 1#########################################################################
def kronecker(A, B):
#   Producto de Kronecker:
#   Implementación del producto de Kronecker de dos matrices A y B de
#   dimensiones mxn y pxq, respectivamete, utlizando la fórmula
#   c_{p(i-1)+k, q(j-1)+l}=a_{ij}b_{kl}. Devuelve A tensor B.
#   1. Inicializar una matriz de (mn)x(pq) con ceros en todas sus entradas
  m = len(A)
  n = len(A[0])
  p = len(B)
  q = len(B[0])

  output = []
  for row in range(m*p):
    output.append([])
    for column in range(n*q):
      output[row].append(0)

  #   2. Aplicar la fórmula para calcular cada elemento de matriz de la matriz
  #      resultante.
  for i in range(m):
    for j in range(n):
      for k in range(p):
        for l in range(q):
          alpha = p*i + k
          beta = q*j + l
          output[alpha][beta] = A[i][j]*B[k][l]
  return output

#FUNCION 2##############################################
def matrix_product(A, B):
#Implementacion del producto matricial
    filas2 = len(B)
    filas1 = len(A)
    columnas1 = len(A[0])
    columnas2 = len(B[0])
    AB=[[0 for column in range(columnas2)] for row in range(filas1)]
    for i in range(filas1):
        for j in range(columnas2):
            for n in range(columnas1):
                AB[i][j]+=A[i][n]*B[n][j]
    return AB

#FUNCION 3 ###################################################################
def matrix_sum(A, B):
#Implementacion del producto matricial
    filas2 = len(B)
    filas1 = len(A)
    columnas1 = len(A[0])
    columnas2 = len(B[0])
    E=[[0 for column in range(columnas1)] for row in range(filas1)]
    for i in range(filas1):
        for j in range(columnas1):
                E[i][j]=A[i][j]+B[i][j]
    return E

#FUNCION 4 ###################################################################
def conmutador(A, B):
#Implementacion de conmutadores
    filas2 = len(B)
    filas1 = len(A)
    columnas1 = len(A[0])
    columnas2 = len(B[0])
    C= matrix_product(A,B)
    D= matrix_product(B,A)
    E=[[0 for column in range(columnas1)] for row in range(filas1)]
    
    for i in range(filas1):
        for j in range(columnas1):
                E[i][j]=C[i][j]-D[i][j]
    return E

#FUNCION 5 ###################################################################
def fancy_print(A):
#Implementacion del producto matricial
    filas1 = len(A)
    columnas1 = len(A[0])
    for fila in A:
        for valor in fila:
            print("\t", valor, end=" ")
        print()

#Definimos las matrices sigma_n^{\alpha} con la función de Kronecker:

#Sigma_1
sigma1x=kronecker(kronecker(paulix, identidad),identidad)
sigma1y=kronecker(kronecker(pauliy, identidad),identidad)
sigma1z=kronecker(kronecker(pauliz, identidad),identidad)

#sigma_2
sigma2x=kronecker(kronecker(identidad, paulix),identidad)
sigma2y=kronecker(kronecker(identidad, pauliy),identidad)
sigma2z=kronecker(kronecker(identidad, pauliz),identidad)

#sigma_3
sigma3x=kronecker(kronecker(identidad, identidad),paulix)
sigma3y=kronecker(kronecker(identidad, identidad),pauliy)
sigma3z=kronecker(kronecker(identidad, identidad),pauliz)

#Primero multilpicamos las matrices que estan en el hamiltoniano:
sigma1x_sigma2x=matrix_product(sigma1x,sigma2x)
sigma1y_sigma2y=matrix_product(sigma1y,sigma2y)
sigma2x_sigma3x=matrix_product(sigma2x,sigma3x)
sigma2y_sigma3y=matrix_product(sigma2y,sigma3y)

#La suma de las matrices anteriores conforman la parte del hamiltoniano que nos interesa
H=matrix_sum(matrix_sum(sigma1x_sigma2x,sigma1y_sigma2y), matrix_sum(sigma2x_sigma3x,sigma2y_sigma3y) )

#Ahora conmutamos H con sigma_n^{alpha}
H_1=conmutador(H,sigma1z)
H_2=conmutador(H, sigma2z)
H_3=conmutador(H,sigma3z)

#La matriz resultante es cero como se quería probar.
fancy_print(matrix_sum(matrix_sum(H_1,H_2),H_3))