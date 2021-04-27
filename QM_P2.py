'''
Autor:  Rubí E. Ramírez Milián
Compliador:
Para compilar: python QM_P2.py
Fecha:  Sat Apr 23 17:38:09 CST 2021
Librerías: 
Entradas, Salidas, Resumen: Se hará un programa para calcular conmutadores. Matrices
'''



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
#Implementacion de  la suma de matrices
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
#Imprime matrices de forma ordenada

    filas1 = len(A)
    columnas1 = len(A[0])
    for fila in A:
        for valor in fila:
            print("\t", valor, end="  ")
        print()

#FUNCION 6 ###################################################################
def scalar_product(n, A):
#Implementacion de producto de una matriz por un escalar
    filas1 = len(A)
    columnas1 = len(A[0])
    E=[[0 for column in range(columnas1)] for row in range(filas1)]
    
    for i in range(filas1):
        for j in range(columnas1):
                E[i][j]=n*A[i][j]
    return E



#######################  PROBLEMA 2 #########################################################################
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
H=matrix_sum(matrix_sum(sigma1x_sigma2x,sigma1y_sigma2y), matrix_sum(sigma2x_sigma3x,sigma2y_sigma3y))

#Ahora conmutamos H con sigma_n^{z}
H_1=conmutador(H,sigma1z)
H_2=conmutador(H, sigma2z)
H_3=conmutador(H,sigma3z)

#La matriz resultante es cero como se quería probar.
fancy_print(matrix_sum(matrix_sum(H_1,H_2),H_3))

#FIN DEL PROBLEMA 2#######################################################################################################

#################################PROBLEMA 3 ##############################################################################

#Definimos los operadores escalera para 8 dimensiones

sigma1_mas=scalar_product(0.5,matrix_sum(sigma1x,scalar_product(1j,sigma1y)))
sigma2_mas=scalar_product(0.5,matrix_sum(sigma2x,scalar_product(1j,sigma2y)))
sigma1_menos=scalar_product(0.5,matrix_sum(sigma1x,scalar_product(-1j,sigma1y)))
sigma2_menos=scalar_product(0.5,matrix_sum(sigma2x,scalar_product(-1j,sigma2y)))
sigma3_mas=scalar_product(0.5,matrix_sum(sigma3x,scalar_product(1j,sigma3y)))
sigma3_menos=scalar_product(0.5,matrix_sum(sigma3x,scalar_product(-1j,sigma3y)))

if matrix_product(sigma1_menos, sigma2_mas)==matrix_product(sigma2_mas, sigma1_menos):
    print(True)
else: 
    print(False)

if matrix_product(sigma2_menos, sigma3_mas)==matrix_product(sigma3_mas, sigma2_menos):
    print(True)
else: 
    print(False)



########## FIN PROBLEMA 3####################################################################################################

#######################  PROBLEMA 4 #########################################################################


#Para el producto tensorial de dos estados definimos los autovectores de pauliz
up=[[1],[0]]
down=[[0],[1]]

#Definimos los valores de los cuatro estados 
up_up= kronecker(up, up)
up_down = kronecker(up,down)
down_up = kronecker(down, up)
down_down= kronecker(down, down)


#matrices de un espacio de cuatro dimensiones
SIGMA1X=kronecker(paulix,identidad)
SIGMA2X=kronecker(identidad,paulix)
SIGMA1Y=kronecker(pauliy, identidad)
SIGMA2Y=kronecker(identidad,pauliy)

#Definamos ahora los operadores escalera

SIGMA1_mas=scalar_product(0.5,matrix_sum(SIGMA1X,scalar_product(1j,SIGMA1Y)))
SIGMA2_mas=scalar_product(0.5,matrix_sum(SIGMA2X,scalar_product(1j,SIGMA2Y)))
SIGMA1_menos=scalar_product(0.5,matrix_sum(SIGMA1X,scalar_product(-1j,SIGMA1Y)))
SIGMA2_menos=scalar_product(0.5,matrix_sum(SIGMA2X,scalar_product(-1j,SIGMA2Y)))

#Definimos h_{n, n+1}
H_1_2=matrix_sum(matrix_product(SIGMA1_mas,SIGMA2_menos),matrix_product(SIGMA2_mas, SIGMA1_menos))


#Definimos h_{n, n+1}
h_1_2=matrix_sum(matrix_product(sigma1_mas,sigma2_menos),matrix_product(sigma2_mas, sigma1_menos))
h_2_3=matrix_sum(matrix_product(sigma2_mas,sigma3_menos),matrix_product(sigma2_menos, sigma3_mas))


#Comprobamos que 2(h_1_2+h_2_3) sea igual al Hamiltoniano
if scalar_product(2,matrix_sum(h_1_2,h_2_3))==H:
    print(True)
else:
    print(False)






#Imprimimos el resultado decalcular el operador h en cada uno de los cuatro estados

print('Primer Estado')
fancy_print(matrix_product(H_1_2,up_up))

print('Segundo Estado')
fancy_print(matrix_product(H_1_2,up_down))


print('Tercer Estado')
fancy_print(matrix_product(H_1_2,down_up))


print('Cuarto Estado')
fancy_print(matrix_product(H_1_2,down_down))

#Definir los estados entrelazados en un espacio de 8 dimensiones
up_up_up=kronecker(up_up,up)
up_up_down=kronecker(up_up,down)
up_down_up=kronecker(up_down,up)
down_up_up=kronecker(down_up,up)
up_down_down=kronecker(up_down,down)
down_down_up=kronecker(down_down,up)
down_up_down=kronecker(down_up,down)
down_down_down=kronecker(down_down,down)



#Visualizamos los estados puros
print('primero')
fancy_print(up_up_up)
print('segundo')
fancy_print(up_up_down)
print('tercero')
fancy_print(up_down_up)
print('cuarto')
fancy_print(up_down_down)
print('quinto')
fancy_print(down_up_up)
print('sexto')
fancy_print(down_up_down)
print('septimo')
fancy_print(down_down_up)
print('octavo')
fancy_print(down_down_down)

#Imprimir los estados

#Aplicamos el hamiltoniano a los estados.

print('primero')
fancy_print(matrix_product(H,up_up_up))
print('segundo')
fancy_print(matrix_product(H,up_up_down))
print('tercero')
fancy_print(matrix_product(H,up_down_up))
print('cuarto')
fancy_print(matrix_product(H,down_up_up))
print('quinto')
fancy_print(matrix_product(H,up_down_down))
print('sexto')
fancy_print(matrix_product(H,down_up_down))
print('septimo')
fancy_print(matrix_product(H,down_down_up))
print('octavo')
fancy_print(matrix_product(H,down_down_down))

##FIN PROBLEMA 4 ############################################################################################

#######################  PROBLEMA 6 #########################################################################


#Producto de sigma1z y sigma2z
productoria=matrix_product(sigma1z,sigma2z)

#Ahora definimos los operadores a_j y a_Daga_j

a_1= sigma1_menos
a_2= matrix_product(sigma1z, sigma2_menos)
a_3= matrix_product(productoria, sigma3_menos)

a_Daga_1= sigma1_mas
a_Daga_2= matrix_product(sigma1z, sigma2_mas)
a_Daga_3= matrix_product(productoria,sigma3_mas)

#Comprobamos que sigma1z y productoria son su propia inversa
fancy_print(matrix_product(sigma1z,sigma1z))
print('\t')
fancy_print(matrix_product(productoria,productoria))




#Comprobamos que a_1(sigma1z) sea igual a a_1
if matrix_product(a_1,sigma1z)==a_1:
    print(True)
else:
    print(False)

print('\t')

#Comprobamos que (sigma1z)a_Daga_1 sea igual a a_Daga_1
if matrix_product(sigma1z,a_Daga_1)==a_Daga_1:
    print(True)
else:
    print(False)

print('\t')

#Comprobamos que a_2(sigma2z) sea igual a a_2
if matrix_product(a_2,sigma2z)==a_2:
    print(True)
else:
    print(False)

print('\t')

#Comprobamos que (sigma2z)a_Daga_2 sea igual a a_Daga_2
if matrix_product(sigma2z,a_Daga_2)==a_Daga_2:
    print(True)
else:
    print(False)

print('\t')

#Se define por partes las sumas de la nueva función 

H1 = matrix_product(a_2,a_Daga_1)
H2 = matrix_product(a_1,a_Daga_2)
H3 = matrix_product(a_3,a_Daga_2)
H4 = matrix_product(a_2,a_Daga_3)

#El hamiltoniano en función de los nuevos operadores escalera.
H_escalera=matrix_sum(matrix_sum(H1,H2),matrix_sum(H3,H4))


#Comprobamos que 2(H_escalera) sea igual al Hamiltoniano
if scalar_product(2,H_escalera)==H:
    print(True)
else:
    print(False)


#### FIN PROBLEMA 6###########################################################################################
