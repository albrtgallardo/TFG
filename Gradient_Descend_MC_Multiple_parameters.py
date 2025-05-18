import numpy as np
import autograd as ag
import matplotlib.pyplot as plt
from collections import Counter
import itertools
import pandas as pd
import time
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox


def ini_vec(inf_range, sup_range):
    return round(np.random.uniform(inf_range,sup_range), 1)

def gradient(v, function):
    grad_f = ag.grad(function)
    return (grad_f(np.array(v)))

def gradient_descend(inf_range, sup_range, n_iter, learning_rate, tolerance, function):
    v = ini_vec(inf_range,sup_range)
    for i in range(n_iter):
        g = gradient(v, function)
        v_new = v - learning_rate*g

        if abs(v_new - v) < tolerance:
            break

        v = v_new
        
    return v

def generar_combinaciones(diccionario):
    claves = diccionario.keys()
    valores = diccionario.values()
    combinaciones = list(itertools.product(*valores))
    return [dict(zip(claves, c)) for c in combinaciones]
    


if __name__ == "__main__":
    T4, T3, T2, T1, T0 = 0.0, 0.0, 0.0, 0.0, 0.0
    input("Comprueba el nombre del archivo csv que se va a guardar")

    Gradient_parameters = {
        "N_iter": [100, 500, 1000, 2000],
        "Min_dist": [0.1, 0.2, 0.25, 0.5, 1],
        "N_mc": [50, 100, 150, 200],
        "Toler": [0.1, 0.05, 0.01, 0.0001],
        "Learning rate": [0.1, 0.05, 0.01, 0.0001],
    }

    opcion = input('1.Lennard-Jones, 2.Higgs, 3.Hidrogeno: ')
    

    if opcion == "1":
        Gradient_parameters = {
            "N_iter": [100, 500, 1000, 2000],
            "N_mc": [50, 100, 150, 200],
            "Toler": [0.1, 0.05, 0.01, 0.0001],
            "Learning rate": [0.1, 0.05, 0.01, 0.0001],
        }
        true_min = 1.1225
        def f(x):
            if x == 0:
                x = 0.01
            y = 4 * 10 * ((0.25/x)**12 - (0.25/x)**6)
            return y

        inf_range = 0.8
        sup_range = 3

    elif opcion == "2":
        true_min = -1
        true_min2 = 1
        def f(x):
            y = -1 * x**2 + 0.5 * x**4
            return y

        inf_range = -2
        sup_range = 2

    elif opcion == "3":
        true_min = 1.1
        true_min2 = -1.06
        def f(x):
            y = 3*(x**4) - 7*(x**2) - 0.5*x + 6
            return y

        inf_range = -2
        sup_range = 2
    else:
        print(f'Selecciona una opcion valida')

    combinaciones = generar_combinaciones(Gradient_parameters)

    contador = 0 
    resultados = []    
    

    for params in combinaciones:
        contador += 1
        print(f'===============Iteracion {contador} de {len(combinaciones)}===================')

        n_iter = params['N_iter']
        num_mc = params['N_mc']
        tolerance = params['Toler']
        learning_rate = params['Learning rate']


    
        soluciones = []
        if opcion == '1':
            print(f"Parametros a introducir: n_MC {num_mc}, n_iter {n_iter}, toler {tolerance}, learning rate {learning_rate}")
            
            start_time = time.time()
            for n in range(num_mc):
                solucion = gradient_descend(inf_range, sup_range, n_iter,  learning_rate, tolerance, f)
                soluciones.append(round(solucion,4))
                
            minimo = Counter(soluciones).most_common(1)[0][0]
            print(f'El valor {minimo} se ha repetido {Counter(soluciones).most_common(1)[0][1]}')
            
            elapsed_time = time.time() - start_time          
            resultado = params.copy()
            resultado["Tiempo ejecucion"] = round(elapsed_time, 4)
            resultado["Primer minimo"] = round(minimo, 4)
            resultado["Δ True Min"] = abs(true_min - minimo)
            
        elif opcion == '2':
            min_dist = params['Min_dist']
            print(f"Parametros a introducir: n_MC {num_mc}, n_iter {n_iter}, min_dist {min_dist}, toler {tolerance}")
            
            start_time = time.time()
            for n in range(num_mc):
                solucion = gradient_descend(inf_range, sup_range, n_iter,  learning_rate, tolerance, f)
                soluciones.append(round(solucion,4))

            cont = Counter(soluciones)
            mas_comunes = cont.most_common()
            primer_min = mas_comunes[0][0]

            for i in range(len(mas_comunes)-1):
                if (abs(primer_min - mas_comunes[i+1][0])) > min_dist:
                    segundo_min = mas_comunes[i+1][0]
                    segundo_min_val = i+1
                    break
            else:
                segundo_min = None
                
            print(f'El primer minimo {primer_min} se ha repetido {Counter(soluciones).most_common(1)[0][1]}')
            print(f'El segundo minimo {segundo_min} se ha repetido {Counter(soluciones).most_common()[segundo_min_val][1]}')
            
            elapsed_time = time.time() - start_time          
            resultado = params.copy()
            resultado["Tiempo ejecucion"] = round(elapsed_time, 4)
            resultado["Primer minimo"] = round(primer_min, 4)
            resultado["Δ True Min"] = min(abs(true_min - primer_min), abs(true_min2 - primer_min))
            
            resultado["Segundo minimo"] = segundo_min
            resultado["Δ True Min2"] = min(abs(true_min - segundo_min), abs(true_min2 - segundo_min))
            
                        
        elif opcion == '3':
            min_dist = params['Min_dist']
            print(f"Parametros a introducir: n_MC {num_mc}, n_iter {n_iter}, min_dist {min_dist}, toler {tolerance}")
            
            start_time = time.time()
            for n in range(num_mc):
                solucion = gradient_descend(inf_range, sup_range, n_iter,  learning_rate, tolerance, f)
                soluciones.append(round(solucion,4))

            cont = Counter(soluciones)
            mas_comunes = cont.most_common()
            primer_min = mas_comunes[0][0]

            for i in range(len(mas_comunes)-1):
                if (abs(primer_min - mas_comunes[i+1][0])) > min_dist:
                    segundo_min = mas_comunes[i+1][0]
                    segundo_min_val = i+1
                    break
            else:
                segundo_min = None
                
            print(f'El primer minimo {primer_min} se ha repetido {Counter(soluciones).most_common(1)[0][1]}')
            print(f'El segundo minimo {segundo_min} se ha repetido {Counter(soluciones).most_common()[segundo_min_val][1]}')
            
            elapsed_time = time.time() - start_time          
            resultado = params.copy()
            resultado["Tiempo ejecucion"] = round(elapsed_time, 4)
            resultado["Primer minimo"] = round(primer_min, 4)
            resultado["Δ True Min"] = min(abs(true_min - primer_min), abs(true_min2 - primer_min))
            
            resultado["Segundo minimo"] = segundo_min
            resultado["Δ True Min2"] = min(abs(true_min - segundo_min), abs(true_min2 - segundo_min))


        print(f'El resultado es {resultado}')

        resultados.append(resultado)

    if resultados:
        df = pd.DataFrame(resultados)
        nombre_archivo = "Gradient_Descend_Lennard_Jones.csv" 
    
        df.to_csv(nombre_archivo, index=False)
        print(f" Resultados guardados en: {nombre_archivo}")
    else:
        print("No se generaron resultados para guardar.")
    
    
