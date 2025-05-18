import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import QasmSimulator
from qiskit.circuit import Parameter
from qiskit.visualization import plot_histogram
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import pandas as pd
import time

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

from qiskit.quantum_info import Statevector
import itertools
from sympy import symbols, Eq, solve


def qubit_states_generator(n_qubits):
    qubit_states = []
    
    if n_qubits < 6: #Cortafuegos
        #Para cada iteracion posible de los valores 01 con una longitud n_qubits se guarda en una lista
        for comb in itertools.product('01', repeat=n_qubits):
            qubit_states.append(''.join(comb))

    return qubit_states


def Z_values_generator(qubit_states, n_qubits):
    
    Z_values = {}
    
    for i in range(len(qubit_states)):
        string = qubit_states[i]
        Z = []
        #parte no combinada de los terminos de Z
        for j in range(len(string)):
            if string[j] == '0':
                Z.append(+1)
            else:
                Z.append(-1)
                
        Z = list(reversed(Z)) #invierto la lista, recordar que los qubits se leen de derecha a izquierda

        
        ############# parte combinada de los terminos de mas a menos ################
        
        #creacion del indice para los terminos cruzados 
        #Esta parte crea la variable combination, que ira sirviendo de indice para crear los terminados cruzados
        combination = '01'
        for j in range(n_qubits-2):
            combination += str(j+2)
            #Para 3 qubits devolera 012, para 4 0123, para 5 01234 etc.


        for i in range(n_qubits-1):
        #Este bucle es para tener todas las combinaciones  de longitud 2 hasta longitud n_qubits
            
            for j in itertools.combinations(combination,i+2):
                j = list(j) #pasamos a una lista los valores de la combinacion generada por itertools
                New_Z = 1 #Se crea el valor que vamos a añadir
                #Ej: si tenemos ('012',2) nos devolvera ['0','1'], ['0','2'] y ['1','2']


                for k in range(len(j)):
                    #Para determinar el valor de ese nuevo resultado, se multiplican los valores que lo forman
                    New_Z *= Z[int(j[k])]
                    
                Z.append(New_Z)
                #finalmente se añade el nuevo resultado

        Z_values[string] = Z

    return Z_values



def solver(Z_values, qubit_states, x_vals, f):
    symbols_list = []
    equations = []
    
    
    for i in range(len(Z_values[qubit_states[0]])+1):
        symbols_list.append(symbols(f"a{i}"))
        # creo una lista de todos los strings que van a ser las variables
        
        
    for i in range(len(symbols_list)): # i es igual a la cantidad de incognitas que voy a crear quizas se pueda cambiar por n_qubits**2
        Ecuacion = symbols_list[0] #Esto lo hago asi, por que el primer simbolo a0 va con la identidad asi que no necestio hacerle ninguna operacion
        for j in range(len(symbols_list)-1):
            Ecuacion += symbols_list[j+1]*Z_values[qubit_states[i]][j]
            #ahora añado el nuevo termino a la ecuacion
        
        eq = Eq(Ecuacion, f(x_vals[i]))#genero la ecuacion de sympy, con su valor respectivo f(x)
        equations.append(eq)           #añado la ecuacion nueva
    
    #esta parte se encagra de pasar a una sola lista de floats la solucion para acceder de manera mas facil
    solution_dic = solve(equations, symbols_list)
    solution = []
    for k in range(len(solution_dic)):
        solution.append(float(solution_dic[symbols_list[k]]))

    return solution


def Z_lengths_iterations(n):
    Z_lenghts = []
    
    Z_lenghts.append(n)
    
    combination = '01'
    for j in range(n-2):
        combination += str(j+2)
        
    
    for i in range(n-1):      
        Z_lenghts.append(len(list(itertools.combinations(combination,i+2))))


    return Z_lenghts


def QAOA_circuit(n_qubits, solution, gamma, beta, p=1):    
    qc = QuantumCircuit(n_qubits)
    qc.h(range(n_qubits))
    #hadammard en todos

    Z_len = Z_lengths_iterations(n_qubits)
    
    for _ in range(p):
        for i in range(Z_len[0]):
            qc.rz(2 * gamma * solution[i+1],i)
        qc.barrier()
        #aplico las puertas Z a los qubits correspondientes
    
        #esta parte es simplemente para tener combination como lo tenia en la parte del Z_values_generator
        combination = '01'
        for j in range(n_qubits-2):
            combination += str(j+2)
    
        posicion_aux = n_qubits+1 #esta va a ser una variable auxiliar, simplemente para poder referirnos a los valores de solucion correctamente
        for i in range(n_qubits-1): # determino las combinaciones que tengo para esos qubits.
            lista = list(itertools.combinations(combination,i+2)) #listo las combinaciones de i puertas Z
            
            #hacemos la parte de los terminos de varios Z de manera general:
            for j in range(len(lista)):
                for k in range(len(lista[0])-1):
                    qc.cx(int(lista[j][k]), int(lista[j][k+1]))
                
                qc.rz(2 * gamma * solution[posicion_aux],int(lista[j][k+1]))
                posicion_aux += 1
                for k in range(len(lista[0])-1):
                    k = len(lista[0])-2 -k # es necesario este cambio para que las rotaciones se deshagan en el orden correcto
                    qc.cx(int(lista[j][k]), int(lista[j][k+1]))
                    
                qc.barrier()
    
        qc.barrier()
    
    
        for q in range(n_qubits):
            qc.rx(2 * beta, q)
        
    qc.measure_all()
    return qc

def QAOA(qubit_states, Z_values, n_qubits, x_vals, beta, gamma, p, f):
    #qubit_states = qubit_states_generator(n_qubits)
    #primero me creo la lista con las combinaciones de los qubits

    #Z_values = Z_values_generator(qubit_states, n_qubits)
    #ahora creo un diccionario con los valores de las puertas Z para cada estado

    solution = solver(Z_values, qubit_states, x_vals, f)
    #calculo la solucion del hamiltoniano y la guardo en una lista

    qc_QAOA = QAOA_circuit(n_qubits, solution,gamma, beta, p) 


    backend = QasmSimulator()
    job = backend.run(qc_QAOA, shots=2048)
    counts = job.result().get_counts()

    
    prob_state = max(counts, key=counts.get)
    prob_val_position = qubit_states.index(prob_state)
    prob_val = x_vals[prob_val_position]
    
    return(prob_val)

def QAOA_MC(n_MC, inf_range, sup_range, qubit_states, Z_values, n_qubits, beta, gamma, p, tolerance, learning_rate, f):
    soluciones = []
    
    for m in range(n_MC):
        #genero una soluciona a partir de un valor aleatorio
        rand_val = round(np.random.uniform(inf_range,sup_range), 1)
        x_vals = [rand_val]
        

        #Creo la lista con los valores iniciales
        for a in range(int((2**n_qubits)/2)):
            x_vals.insert(0, round(rand_val - (0.05*(a+1)),2))
            x_vals.append(round(rand_val +(0.05*(a+1)),2))
        x_vals.pop(2**n_qubits)
        
        QAOA_solution = QAOA(qubit_states, Z_values, n_qubits, x_vals, beta, gamma, p, f)
        
        #genero otra solucion de entre unos valores, cuyo rango esta centrado en la solucion anterior
        x_vals = [QAOA_solution]
        for a in range(int((2**n_qubits)/2)):
            x_vals.insert(0, round(QAOA_solution - (0.05*(a+1)),2))
            x_vals.append(round(QAOA_solution +(0.05*(a+1)),2))
        x_vals.pop(2**n_qubits)
        QAOA_solution2 = QAOA(qubit_states, Z_values, n_qubits, x_vals, beta, gamma, p, f)

        
        #learning_rate = 0.5
        
        
        
        n_points = n_qubits**2
        center_index = n_points // 2  # índice donde estará QAOA_solution // division entera

        for i in range(30):
            if abs(QAOA_solution - QAOA_solution2) < tolerance:
                soluciones.append(QAOA_solution2)
                break

            QAOA_solution = QAOA_solution2
            mu = learning_rate / (1+ (i))
            x_vals = [
                round(QAOA_solution + mu * (j - center_index), 2)
                for j in range(n_points)
            ]
            QAOA_solution2 = QAOA(qubit_states, Z_values, n_qubits, x_vals, beta, gamma, p, f)

    return soluciones

def generar_combinaciones(diccionario):
    claves = diccionario.keys()
    valores = diccionario.values()
    combinaciones = list(itertools.product(*valores))
    return [dict(zip(claves, c)) for c in combinaciones]


if __name__ == "__main__":
    T4, T3, T2, T1, T0 = 0.0, 0.0, 0.0, 0.0, 0.0
    input("Comprueba el nombre del archivo csv que se va a guardar")
    
    QAOA_parameters = {
        "N_qubits": [2, 3, 4],
        "beta": [0.5, 0.3, 0.1],
        "gamma": [0.5, 0.25, 0.1, 0.05, 0.01],
        "Prof_p": [2, 3, 4],
        "N_mc": [50, 100],
        "Toler": [0.05],
        "Learning rate": [0.5],
    }
    
    opcion = input('1.Lennard-Jones, 2.Higgs, 3.Hidrogeno: ')
    
    if opcion == "1":
        true_min = 1.1225
        def f(x):
            if x == 0:
                x = 0.01
            y = 4 * 10 * ((0.25/x)**12 - (0.25/x)**6)
            return -y
        
        inf_range = 0.8
        sup_range = 3

    elif opcion == "2":
        true_min = -1
        true_min2 = 1
        def f(x):
            y = -1 * x**2 + 0.5 * x**4
            return -y
            
        inf_range = -2
        sup_range = 2

    elif opcion == "3":
        true_min = 1.1
        true_min2 = -1.06
        def f(x):
            y = 3*(x**4) - 7*(x**2) - 0.5*x + 6
            return -y
            
        inf_range = -2
        sup_range = 2
    else:
        print(f'Selecciona una opcion valida')
        
    combinaciones = generar_combinaciones(QAOA_parameters)
    
    contador = 0 
    resultados = []
    for params in combinaciones:
        
        contador += 1
        print(f'===============Iteracion {contador} de {len(combinaciones)}===================')
 
        n_qubits = params["N_qubits"]
        beta = params["beta"]
        gamma = params["gamma"]
        p = params["Prof_p"]
        n_MC = params["N_mc"]
        toler = params["Toler"]
        learn = params["Learning rate"]     


        print(f"Parametros a introducir: n_MC {n_MC}, n_qubits {n_qubits}, beta {beta}, gamma {gamma}, p {p}, toler {toler}, learn {learn}")
        
        qubit_states = qubit_states_generator(n_qubits)
        Z_values = Z_values_generator(qubit_states, n_qubits)


        soluciones_lista = []
        start_time = time.time()
       
        soluciones_lista = QAOA_MC(n_MC, inf_range, sup_range, qubit_states, Z_values, n_qubits, beta, gamma, p, toler, learn, f)
        
        elapsed_time = time.time() - start_time
        
        x = np.arange(inf_range, sup_range, 0.1)
        
        conteos, bordes = np.histogram(soluciones_lista, bins=x)
        conteos = list(conteos)
        
        max_conteos = conteos.index(max(conteos))
        primer_min = bordes[max_conteos]

        if opcion == '1':
            resultado = params.copy()
            resultado["Tiempo ejecucion"] = round(elapsed_time, 4)
            resultado["Primer minimo"] = round(primer_min, 4)
            resultado["Δ True Min"] = abs(true_min - primer_min)
    
            print(f'Minimo es {primer_min} y el minimo real {true_min}')
        
        elif opcion == '2':
            resultado = params.copy()
            resultado["Tiempo ejecucion"] = round(elapsed_time, 4)
            resultado["Primer minimo"] = round(primer_min, 4)
            resultado["Δ True Min"] = min(abs(true_min - primer_min), abs(true_min2 - primer_min))

            print(f'Minimo es {primer_min} y el minimo real ±1')

            conteos.pop(max_conteos)
            max_conteos2 = conteos.index(max(conteos))
            if abs(primer_min - bordes[max_conteos2]) < 0.25:
                conteos.pop(max_conteos2)
                max_conteos2 = conteos.index(max(conteos))
                
            segundo_min = bordes[max_conteos2]

            resultado["Segundo minimo"] = segundo_min
            resultado["Δ True Min2"] = min(abs(true_min - segundo_min), abs(true_min2 - segundo_min))
            
            print(f'Segundo minimo es {segundo_min} y el minimo real ±1')

        elif opcion == '3':
            resultado = params.copy()
            resultado["Tiempo ejecucion"] = round(elapsed_time, 4)
            resultado["Primer minimo"] = round(primer_min, 4)
            resultado["Δ True Min"] = min(abs(true_min - primer_min), abs(true_min2 - primer_min))
    
            print(f'Minimo es {primer_min} y el minimo real {true_min}')

            conteos.pop(max_conteos)
            max_conteos2 = conteos.index(max(conteos))
            if abs(primer_min - bordes[max_conteos2]) < 0.25:
                conteos.pop(max_conteos2)
                max_conteos2 = conteos.index(max(conteos))
                
            segundo_min = bordes[max_conteos2]

            resultado["Segundo minimo"] = segundo_min
            resultado["Δ True Min2"] = min(abs(true_min - segundo_min), abs(true_min2 - segundo_min))
            
            print(f'Segundo minimo es {segundo_min} y el minimo real {true_min2}')


        print(f'El resultado es {resultado}')

        resultados.append(resultado)
    
    if resultados:
        df = pd.DataFrame(resultados)
        nombre_archivo = "QAOA_Hidrogeno.csv" 
    
        df.to_csv(nombre_archivo, index=False)
        print(f" Resultados guardados en: {nombre_archivo}")
    else:
        print("No se generaron resultados para guardar.")


    
    
