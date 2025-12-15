#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Rob�tica Computacional -
# Grado en Ingenier�a Inform�tica (Cuarto)
# Pr�ctica: Filtros de particulas.

from math import *
from robot import *
import random
import numpy as np
import matplotlib.pyplot as plt
import sys
import select
from datetime import datetime
# ******************************************************************************
# Declaraci�n de funciones

def distancia(a,b):
  # Distancia entre dos puntos (admite poses)
  return np.linalg.norm(np.subtract(a[:2],b[:2]))

def angulo_rel(pose,p):
  # Diferencia angular entre una pose y un punto objetivo 'p'
  w = atan2(p[1]-pose[1],p[0]-pose[0])-pose[2]
  while w >  pi: w -= 2*pi
  while w < -pi: w += 2*pi
  return w

def pinta(secuencia,args):
  # Dibujar una secuencia de puntos
  t = np.array(secuencia).T.tolist()
  plt.plot(t[0],t[1],args)

def mostrar(objetivos,trayectoria,trayectreal,filtro, step_mode=False):
  # Mostrar mapa y trayectoria
  plt.ion() # modo interactivo
  plt.clf()
  # Fijar los bordes del gráfico
  objT   = np.array(objetivos).T.tolist()
  bordes = [min(objT[0]),max(objT[0]),min(objT[1]),max(objT[1])]
  centro = [(bordes[0]+bordes[1])/2.,(bordes[2]+bordes[3])/2.]
  radio  = max(bordes[1]-bordes[0],bordes[3]-bordes[2])
  plt.xlim(centro[0]-radio,centro[0]+radio)
  plt.ylim(centro[1]-radio,centro[1]+radio)
  plt.gca().set_aspect('equal', adjustable='box')
  # Representar mapa
  for p in filtro:
    dx = cos(p.orientation)*.05
    dy = sin(p.orientation)*.05
    plt.arrow(p.x,p.y,dx,dy,head_width=.05,head_length=.05,color='k')
  pinta(trayectoria,'--g')
  pinta(trayectreal,'-r')
  pinta(objetivos,'-.ob')
  p = hipotesis(filtro)
  dx = cos(p[2])*.05
  dy = sin(p[2])*.05
  plt.arrow(p[0],p[1],dx,dy,head_width=.075,head_length=.075,color='m')
  # Mostrar y comprobar pulsaciones de teclado:
  plt.show()
  # Si estamos en modo paso, no forzamos un pequeño delay continuo;
  # el avance se controlará desde el bucle principal mediante `input()`.
  plt.draw()
  plt.pause(0.01)


# GENERACION ALEATORIA DE LAS PARTICULAS DEL FILTRO
def genera_filtro(num_particulas, balizas, real, centro=[2,2], radio=3):
  # Inicializaci�n de un filtro de tama�o 'num_particulas', cuyas part�culas
  # imitan a la muestra dada y se distribuyen aleatoriamente sobre un �rea dada.

  # Creamos el array de robots
  robots = []
  for i in range(num_particulas):
    # Generamos para cada partícula un valor aleatorio entre el círculo inicial
    # y un radio. Para la orientación del robot también generamos un valor aleatorio
    x_random_value = random.random()
    y_random_value = random.random()

    if (x_random_value > 0.5):
      x_random = centro[0] + random.random() * radio
    else:
      x_random = centro[0] - random.random() * radio
    
    if (y_random_value > 0.5):
      y_random = centro[1] + random.random() * radio
    else:
      y_random = centro[1] - random.random() * radio

    orientation_random = random.random() * 2 * pi
    # Creamos el robot y calculamos su peso asociado.
    new_robot = robot()
    new_robot.set(x_random, y_random, orientation_random)
    new_robot.set_noise(.01,.01,.01)
    new_robot.measurement_prob(real.sense(balizas), balizas)
    # Añadimos el robot al conjunto.
    robots.append(new_robot)
  return robots

def dispersion(filtro):
  # Devuelve el max y min de las posiciones de las partículas del filtro
  xs = [r.x for r in filtro]
  ys = [r.y for r in filtro]
  return [max(xs)-min(xs), max(ys)-min(ys)]


# SE CALCULA LA MEDIA DE LAS POSICIONES
def peso_medio(filtro):
  total = sum(p.weight for p in filtro)
  for p in filtro:
    p.weight /= total

# ******************************************************************************

random.seed(0)

# Definici�n del robot:
P_INICIAL = [0.,4.,0.] # Pose inicial (posici�n y orientacion)
V_LINEAL  = .7         # Velocidad lineal    (m/s)
V_ANGULAR = 140.       # Velocidad angular   (�/s)
FPS       = 10.        # Resoluci�n temporal (fps)
HOLONOMICO = 0         # Robot holon�mico
GIROPARADO = 0         # Si tiene que tener vel. lineal 0 para girar
LONGITUD   = .1        # Longitud del robot

N_PARTIC  = 100         # Tama�o del filtro de part�culas
N_INICIAL = 2000       # Tama�o inicial del filtro

# Modo paso: si True, la simulación espera a que pulses Enter antes de
# ejecutar cada frame. Escribe 'c' + Enter para salir del modo paso.
STEP_MODE = True

# Definici�n de trayectorias:
trayectorias = [
    [[0,2],[4,2]],
    [[2,4],[4,0],[0,0]],
    [[2,4],[2,0],[0,2],[4,2]],
    [[2+2*sin(.4*pi*i),2+2*cos(.4*pi*i)] for i in range(5)],
    [[2+2*sin(.8*pi*i),2+2*cos(.8*pi*i)] for i in range(5)],
    [[2+2*sin(1.2*pi*i),2+2*cos(1.2*pi*i)] for i in range(5)],
    [[2*(i+1),4*(1+cos(pi*i))] for i in range(6)],
    [[2+.2*(22-i)*sin(.1*pi*i),2+.2*(22-i)*cos(.1*pi*i)] for i in range(20)],
    [[2+(22-i)/5*sin(.1*pi*i),2+(22-i)/5*cos(.1*pi*i)] for i in range(20)]
    ]

# Definici�n de los puntos objetivo:
if len(sys.argv)<2 or int(sys.argv[1])<0 or int(sys.argv[1])>=len(trayectorias):
  sys.exit(sys.argv[0]+" <indice entre 0 y "+str(len(trayectorias)-1)+">")
objetivos = trayectorias[int(sys.argv[1])]

# Definici�n de constantes:
EPSILON = .1                # Umbral de distancia
V = V_LINEAL/FPS            # Metros por fotograma
W = V_ANGULAR*pi/(180*FPS)  # Radianes por fotograma

real = robot()
real.set_noise(.01,.01,.01) # Ruido lineal / radial / de sensado
real.set(*P_INICIAL)

#inicializaci�n del filtro de part�culas y de la trayectoria
#------------------------------------------------------------------------------
initial_pos = [P_INICIAL[0], P_INICIAL[1]]
filtro_particulas = genera_filtro(N_INICIAL, objetivos, real, initial_pos, 1.5)

trayectreal = [real.pose()]
trayectoria = [hipotesis(filtro_particulas)]
#------------------------------------------------------------------------------


tiempo  = 0.
espacio = 0.
for punto in objetivos:
  while distancia(trayectoria[-1],punto) > EPSILON and len(trayectoria) <= 1000:
    # Si estamos en modo paso, esperar a que el usuario pulse Enter antes
    # de ejecutar el siguiente fotograma. Si introduce 'c' se sale al modo continuo.
    if STEP_MODE:
      try:
        cmd = input("Pulsa Enter para avanzar un frame (o escribe 'c' y Enter para continuar): ")
      except EOFError:
        cmd = ''
      if cmd.strip().lower() == 'c':
        STEP_MODE = False
    #seleccionar pose
    # Escogemos el mejor robot del filtro
    pose = hipotesis(filtro_particulas)

    # Movemos todos los robots en base al mejor robot.
    w = angulo_rel(pose,punto)
    if w > W:  w =  W
    if w < -W: w = -W
    v = distancia(pose,punto)
    if (v > V): v = V
    if (v < 0): v = 0
    if HOLONOMICO:
      if GIROPARADO and abs(w) > .01:v = 0
      real.move(w,v)
    else:
      real.move_triciclo(w,v,LONGITUD)
      # Movemos las part�culas del filtro
      for i in range(len(filtro_particulas)):
        filtro_particulas[i].move_triciclo(w, v, LONGITUD)

    # Seleccionar hip�tesis de localizaci�n y actualizar la trayectoria
    # Añadimos a la trayectoria la mejor calculada y la real.
    trayectoria.append(hipotesis(filtro_particulas))
    trayectreal.append(real.pose())
    mostrar(objetivos,trayectoria,trayectreal,filtro_particulas, step_mode=STEP_MODE)

    # remuestreo
    filtro_particulas = resample(filtro_particulas, N_PARTIC)
    #recalculamos el peso de cada robot (la función los pone a 1)
    for particle in filtro_particulas:
      particle.measurement_prob(real.sense(objetivos), objetivos)

    espacio += v
    tiempo  += 1

if len(trayectoria) > 1000:
  print ("<< ! >> Puede que no se haya alcanzado la posicion final.")
print ("Recorrido: "+str(round(espacio,3))+"m / "+str(tiempo/FPS)+"s" )
print ("Error medio de la trayectoria: "+str(round(sum(\
    [distancia(trayectoria[i],trayectreal[i])\
    for i in range(len(trayectoria))])/tiempo,3))+"m" )
input()

