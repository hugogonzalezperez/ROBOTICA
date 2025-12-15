#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - 
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure()
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.pause(1.5)
  plt.show()
  plt.close()

def matriz_T(d,th,a,al):
   
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD
#plt.ion() # modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 4:
  sys.exit("python " + sys.argv[0] + " x y <config_file.txt>")
config_file = open(sys.argv[3])

linea = config_file.readline()
while linea[0] == '#':
  linea = config_file.readline()

# valores articulares arbitrarios para la cinemática directa inicial
# Leemos el número de articulaciones que tiene el brazo
num_articulaciones=int(linea[0])

# Ángulos de rotación inicial de las articulaciones
linea=config_file.readline()
#th=[0.,0.,0.]
th = linea.split() # Separa la linea por espacios y los guarda en th
if len(th) != num_articulaciones:
  sys.exit('Error: el número de ángulos iniciales de las articulaciones no coincide con el número de articulaciones')

for i in range(len(th)):
  th[i] = float(th[i])


# Distancia inicial entre las articulaciones
linea=config_file.readline()
#a =[5.,5.,5.] # Distancia entre las articulaciones
a = linea.split()
if len(a) != num_articulaciones:
  sys.exit('Error: el número de distancias iniciales entre articulaciones no coincide con el número de articulaciones')

for i in range(len(a)):
  a[i] = float(a[i])


# Límite superior de rotación de las articulaciones
linea=config_file.readline()
limite_sup_rotacion = linea.split()
if len(limite_sup_rotacion) != num_articulaciones:
  sys.exit('Error: el número de límites superiores de rotación no coincide con el número articulaciones')

for i in range(len(limite_sup_rotacion)):
  limite_sup_rotacion[i] = float(limite_sup_rotacion[i])


# Límite inferior de rotacion de las articulaciones
linea=config_file.readline()
limite_inf_rotacion = linea.split()
if len(limite_inf_rotacion) != num_articulaciones:
  sys.exit('Error: el número de límites inferiores de rotación no coincide con el número articulaciones')

for i in range(len(limite_inf_rotacion)):
  limite_inf_rotacion[i] = float(limite_inf_rotacion[i])


# Tipo de articulacion, Rotacion = 0, Prismática = 1
linea=config_file.readline()
tipo_articulacion = linea.split()

if len(tipo_articulacion) != num_articulaciones:
  sys.exit('Error: el número de tipos de articulaciones no coincide con el número articulaciones')

for i in range(len(tipo_articulacion)):
  if tipo_articulacion[i] != '0' and tipo_articulacion[i] != '1':
    sys.exit('Tipo de articulación ' + tipo_articulacion[i] + ' incorrecto, debe ser 0 o 1, Rotación = 0, Prismática = 1')
  tipo_articulacion[i] = int(tipo_articulacion[i])


# Límite superior de articulación prismática
linea=config_file.readline()
limite_sup_prismatica = linea.split()
if len(limite_sup_prismatica) != num_articulaciones:
  sys.exit('Error: el número de límites superiores de distancia entre articulaciones no coincide con el número articulaciones')

for i in range(len(limite_sup_prismatica)):
  limite_sup_prismatica[i] = float(limite_sup_prismatica[i])


# Límite inferior de articulación prismática
linea=config_file.readline()
limite_inf_prismatica = linea.split()
if len(limite_inf_prismatica) != num_articulaciones:
  sys.exit('Error: el número de límites inferiores de distancia entre articulaciones no coincide con el número articulaciones')

for i in range(len(limite_inf_prismatica)):
  limite_inf_prismatica[i] = float(limite_inf_prismatica[i])

L = sum(a) # variable para representación gráfica
EPSILON = .01
objetivo=[float(i) for i in sys.argv[1:len(sys.argv) - 1]]

# Incremento L de 5 en 5 hasta que el punto objetivo quede dentro de la gráfica
while ((L <= objetivo[0]) | (L <= objetivo[1])):
  L += 5

O=cin_dir(th,a)
#O=zeros(len(th)+1) # Reservamos estructura en memoria
# Calculamos la posicion inicial
print ("- Posicion inicial:")
muestra_origenes(O)

dist = float("inf")
prev = 0.
iteracion = 1
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  O=[cin_dir(th,a)]
  # Para cada combinación de articulaciones:
  for i in range(len(th)):
    if tipo_articulacion[len(th) - 1 - i] == 0:
      # cálculo de la cinemática inversa:
      Px, Py = O[i][len(th) - 1 - i] # coordenadas (x, y) de la articulacion que se esta estudiando
      EFx, EFy = O[i][-1] # coordenadas (x, y) del EF
      Tx, Ty = objetivo # coordenadas (x, y) del punto objetivo

      alpha1 = atan2(EFy - Py, EFx - Px)
      alpha2 = atan2(Ty - Py, Tx - Px)
      inc_th = alpha2 - alpha1
      new_th = th[len(th) - 1 - i] + inc_th
      # Normalizar new_th (entre -pi y pi)y después comprobar que no se pase de los límites para los límites hacer un vector del mismo tamaño que th y poner ahi los límites
      if new_th > pi:
        new_th -= 2 * pi
      elif new_th < pi:
        new_th += 2 * pi

      if new_th > limite_sup_rotacion[len(th) - 1 - i]:
        new_th = limite_sup_rotacion[len(th) - 1 - i]
      elif new_th < limite_inf_rotacion[len(th) - 1 - i]:
        new_th = limite_inf_rotacion[len(th) - 1 - i]

      th[len(th) - 1 - i] = new_th
    elif tipo_articulacion[len(th) - 1 - i] == 1:
      # Calcular d = v * u, v = (Xt - Xp, Yt - Yp), u = (cos(w), sen(w)), w =  Sumarorio de los tita anteriores, incluido el de la propia articulacion
      w = sum(th[len(th) - 1 - i: + 1])
      d = np.dot(np.subtract(objetivo, O[i][-1]), [cos(w), sin(w)])
      a[len(th) - 1 - i] += d
      if a[len(th) - 1 - i] > limite_sup_prismatica[len(th) - 1 - i]:
        a[len(th) - 1 - i] = limite_sup_prismatica[len(th) - 1 - i]
      elif a[len(th) - 1 - i] < limite_inf_prismatica[len(th) - 1 - i]:
        a[len(th) - 1 - i] = limite_inf_prismatica[len(th) - 1 - i]
  
    O.append(cin_dir(th,a))

  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print ("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")
print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))
