import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import QasmSimulator
from qiskit.circuit import Parameter
from qiskit.visualization import plot_histogram
from scipy.optimize import minimize
import matplotlib.pyplot as plt

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
        
        #print(Ecuacion) #descomentar esta linea para ver las ecuacuiones
        eq = Eq(Ecuacion, f(x_vals[i]))#genero la ecuacion de sympy, con su valor respectivo f(x)
        equations.append(eq)           #añado la ecuacion nueva
    
    #esta parte se encagra de pasar a una sola lista de floats la solucion para acceder de manera mas facil
    solution_dic = solve(equations, symbols_list)
    solution = []
    for k in range(len(solution_dic)):
        solution.append(float(solution_dic[symbols_list[k]]))
    #print(solution) #descomentar esta linea para ver las soluciones finalmente

    return solution


def Z_lengths_iterations(n):
    Z_lenghts = []
    
    Z_lenghts.append(n)
    
    combination = '01'
    for j in range(n-2):
        combination += str(j+2)
        
    
    for i in range(n-1):      
        Z_lenghts.append(len(list(itertools.combinations(combination,i+2))))
        #print(i)


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
            #print(lista) #eliminar este print para ver las combinaciones que hace
            
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
    progress["value"] = 0
    progress["maximum"] = n_MC

    
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
                #print(f':La solucion es {QAOA_solution2}')
                soluciones.append(QAOA_solution2)
                break

            QAOA_solution = QAOA_solution2
            mu = learning_rate / (1+ (i))
            print(mu,QAOA_solution)
            x_vals = [
                round(QAOA_solution + mu * (j - center_index), 2)
                for j in range(n_points)
            ]
            print(x_vals)
            QAOA_solution2 = QAOA(qubit_states, Z_values, n_qubits, x_vals, beta, gamma, p, f)

        progress["value"] += 1
        ventana.update_idletasks()
        
    return soluciones
            
#====================FUNCIONES PARA LA INTERFAZ ===========================#


def fijar_valores():
    if usar_valores_predefinidos.get():
        beta_val.delete(0, tk.END)
        gamma_val.delete(0, tk.END)
        p_val.delete(0, tk.END)
        toler_val.delete(0, tk.END)
        learning_val.delete(0, tk.END)
        beta_val.insert(0, "0.3")
        gamma_val.insert(0, "0.05")
        p_val.insert(0, "3")
        toler_val.insert(0, "0.05")
        learning_val.insert(0, "0.5")
        beta_val.config(state='disabled')
        gamma_val.config(state='disabled')
        p_val.config(state='disabled')
        toler_val.config(state='disabled')
        learning_val.config(state='disabled')
    else:
        beta_val.config(state='normal')
        gamma_val.config(state='normal')
        p_val.config(state='normal')
        toler_val.config(state='normal')
        learning_val.config(state='normal')
        beta_val.delete(0, tk.END)
        gamma_val.delete(0, tk.END)
        p_val.delete(0, tk.END)
        toler_val.delete(0, tk.END)
        learning_val.delete(0, tk.END)


############################ BOTON DE PLOTEO #####################################

def plot_button():
    opcion = opcion_var.get()
    if opcion == "Cualquiera":
        try:
            inf_range = float(inf.get())
            sup_range = float(sup.get())
            
            T4 = float(V_T4.get()) 
            T3 = float(V_T3.get())
            T2 = float(V_T2.get())
            T1 = float(V_T1.get())
            T0 = float(V_T0.get())
        
            
            x = np.arange(inf_range, sup_range, 0.01)
            y = T4*(x**4) + T3*(x**3) + T2*(x**2) + T1*x + T0
            plt.figure()
            plt.plot(x, y)
            plt.grid()
            plt.show()
    
        except:
            print("Error: falta algun campo")
            
    elif opcion == "Lennard-Jones":
        try:
            inf_range = 0.8
            sup_range = 3

            x = np.arange(inf_range, sup_range, 0.01)
            y =  4 * ((1/x)**12 - (1/x)**6)
            plt.figure()
            plt.plot(x, y)
            plt.grid()
            plt.show()
    
        except:
            print("Error: falta algun campo")   
    
    elif opcion == "Higgs":
        try:
            inf_range = -2
            sup_range = 2

            x = np.arange(inf_range, sup_range, 0.01)
            y = -1 * x**2 + 0.5 * x**4
            plt.figure()
            plt.plot(x, y)
            plt.grid()
            plt.show()
    
        except:
            print("Error: falta algun campo")

    elif opcion == "Hidrogeno":
        try:
            inf_range = -2
            sup_range = 2

            x = np.arange(inf_range, sup_range, 0.01)
            y = 3*(x**4) - 7*(x**2) - 0.5*x + 6
            plt.figure()
            plt.plot(x, y)
            plt.grid()
            plt.show()
    
        except:
            print("Error: falta algun campo")
    

################################# BOTON DE EJECUCION ######################################
def ejecutar_MC():
    global T4, T3, T2, T1, T0
    n_qubits = int(entry_qubits.get())
    n_MC = int(entry_mc_iteraciones.get())

    toler = float(toler_val.get())
    learn = float(learning_val.get())
    
    beta = float(beta_val.get())
    gamma = float(gamma_val.get())
    p = int(p_val.get())

    opcion = opcion_var.get()
    
    if opcion == "Cualquiera":
        def f(x):
            y = T4*(x**4) + T3*(x**3) +T2*(x**2) + T1*x + T0
            return -y
            
        inf_range = 0.8
        sup_range = 3
    
        T4 = float(V_T4.get()) 
        T3 = float(V_T3.get())
        T2 = float(V_T2.get())
        T1 = float(V_T1.get())
        T0 = float(V_T0.get())

    elif opcion == "Lennard-Jones":
        def f(x):
            if x == 0:
                x = 0.01
            y = 4 * 10 * ((0.25/x)**12 - (0.25/x)**6)
            return -y
        
        inf_range = 0.8
        sup_range = 3

    elif opcion == "Higgs":
        def f(x):
            y = -1 * x**2 + 0.5 * x**4
            return -y
            
        inf_range = -2
        sup_range = 2

    elif opcion == "Hidrogeno":
        def f(x):
            y = 3*(x**4) - 7*(x**2) - 0.5*x + 6
            return -y
            
        inf_range = -2
        sup_range = 2

    
    qubit_states = qubit_states_generator(n_qubits)
    Z_values = Z_values_generator(qubit_states, n_qubits)

    start_time = time.time()    
    soluciones_lista = []
    
    soluciones_lista = QAOA_MC(n_MC, inf_range, sup_range, qubit_states, Z_values, n_qubits, beta, gamma, p, toler, learn, f)
    print(soluciones_lista)
    elapsed_time = time.time() - start_time
    tiempo_label.config(text=f"Tiempo de ejecucion: {elapsed_time:.2f} segundos")

    x = np.arange(inf_range, sup_range, 0.1)
    plt.figure()
    plt.hist(soluciones_lista, bins=x, color='skyblue', edgecolor='black')
    plt.title("Distribución de soluciones (QAOA + MC)")
    plt.xlabel("Valores de x")
    plt.ylabel("Frecuencia")
    plt.show(block=False)

    conteos, bordes = np.histogram(soluciones_lista, bins=x)
    conteos = list(conteos)

    max_conteos = conteos.index(max(conteos))
    primer_minimo = bordes[max_conteos]
    resultado_label.config(text=f"Primer minimo: {primer_minimo}")

    conteos.pop(max_conteos)
    
    max_conteos2 = conteos.index(max(conteos))
    if abs(primer_minimo - bordes[max_conteos2]) < 0.5:
        conteos.pop(max_conteos2)
        max_conteos2 = conteos.index(max(conteos))
    
    resultado_label2.config(text=f"Segundo minimo: {bordes[max_conteos2]}")


######################## FUNCIONES DEL MENU DESPLEGABLE ###########################
""" 
def no_function():
    function_def.config(text="Función Polinomio: ax^n + ... + c")
    f = g

def lennar_jones_function():
    function_def.config(text="Función Lennard-Jones: 4ε[(σ/r)^12 - (σ/r)^6]")
    f = lennard_jones

def higgs_function():
    function_def.config(text="Función Higgs: V(ϕ)=μ^2 ϕ^2 + λϕ^4")
    f = higgs

def hidrogeno_function():
    function_def.config(text="Hidrogeno")
    f = hidrogeno
"""
'''
def opcion_seleccionada(event):
    opcion = opcion_var.get()
    if opcion == "Cualquiera":
        def f(x):
            y = T4*(x**4) + T3*(x**3) +T2*(x**2) + T1*x + T0
            return -y

    elif opcion == "Lennard-Jones":
        def f(x):
            y = 4 * 10 * ((0.25/x)**12 - (0.25/x)**6)
            return -y

    elif opcion == "Higgs":
        def f(x):
            y = -1 * x**2 + 0.5 * x**4
            return -y

    elif opcion == "Hidrogeno":
        def f(x):
            y = 3*(x**4) - 7*(x**2) - 0.5*x + 6
            return -y
''' 

'''
def g(x):
    y = T4*(x**4) + T3*(x**3) +T2*(x**2) + T1*x + T0
    return -y

def lennard_jones(x):
    y = 4 * 10 * ((0.25/x)**12 - (0.25/x)**6)
    return -y

def higgs(x):
    y = -1 * x**2 + 0.5 * x**4
    return -y

def hidrogeno(x):
    y = 3*(x**4) - 7*(x**2) - 0.5*x + 6
    return -y
'''

def opcion_seleccionada(event):
    opcion = opcion_var.get()
    if opcion == "Cualquiera":
        function_def.config(text="Función Polinomio: ax^n + ... + c")
        inf.config(state='normal')
        sup.config(state='normal')
        inf.delete(0, tk.END)
        sup.delete(0, tk.END)

        V_T4.config(state='normal')
        V_T3.config(state='normal')
        V_T2.config(state='normal')
        V_T1.config(state='normal')
        V_T0.config(state='normal')
        V_T4.delete(0, tk.END)
        V_T3.delete(0, tk.END)
        V_T2.delete(0, tk.END)
        V_T1.delete(0, tk.END)
        V_T0.delete(0, tk.END)
        
        #no_function()
    elif opcion == "Lennard-Jones":
        function_def.config(text="Función Lennard-Jones: 4ε[(σ/r)^12 - (σ/r)^6]")
        inf.config(state='disabled')
        sup.config(state='disabled')
        inf.delete(0, tk.END)
        sup.delete(0, tk.END)

        V_T4.config(state='disabled')
        V_T3.config(state='disabled')
        V_T2.config(state='disabled')
        V_T1.config(state='disabled')
        V_T0.config(state='disabled')
        V_T4.delete(0, tk.END)
        V_T3.delete(0, tk.END)
        V_T2.delete(0, tk.END)
        V_T1.delete(0, tk.END)
        V_T0.delete(0, tk.END)

        #lennar_jones_function()
    elif opcion == "Higgs":
        function_def.config(text="Función Higgs: V(ϕ)=μ^2 ϕ^2 + λϕ^4")
        inf.config(state='disabled')
        sup.config(state='disabled')
        inf.delete(0, tk.END)
        sup.delete(0, tk.END)

        V_T4.config(state='disabled')
        V_T3.config(state='disabled')
        V_T2.config(state='disabled')
        V_T1.config(state='disabled')
        V_T0.config(state='disabled')
        V_T4.delete(0, tk.END)
        V_T3.delete(0, tk.END)
        V_T2.delete(0, tk.END)
        V_T1.delete(0, tk.END)
        V_T0.delete(0, tk.END)
        
        #higgs_function()

    elif opcion == "Hidrogeno":
        function_def.config(text="Hidrogeno")
        inf.config(state='disabled')
        sup.config(state='disabled')
        inf.delete(0, tk.END)
        sup.delete(0, tk.END)

        V_T4.config(state='disabled')
        V_T3.config(state='disabled')
        V_T2.config(state='disabled')
        V_T1.config(state='disabled')
        V_T0.config(state='disabled')
        V_T4.delete(0, tk.END)
        V_T3.delete(0, tk.END)
        V_T2.delete(0, tk.END)
        V_T1.delete(0, tk.END)
        V_T0.delete(0, tk.END)
        
        #hidrogeno_function()

####################################INICIO DEL CODIGO################################
if __name__ == "__main__":
    T4, T3, T2, T1, T0 = 0.0, 0.0, 0.0, 0.0, 0.0

    ventana = tk.Tk()
    ventana.title("QAOA con Monte Carlo")
    ventana.geometry("800x400")
    
    # ==================== SECCIÓN: Parámetros del Algoritmo ====================
    frame_parametros = tk.LabelFrame(ventana, text="Parámetros del algoritmo", padx=10, pady=10)
    frame_parametros.grid(row=0, column=0, padx=10, pady=10, sticky="nsew", columnspan=3)
    
    tk.Label(frame_parametros, text="Número de iteraciones (MC):").grid(row=0, column=0, sticky='w')
    entry_mc_iteraciones = tk.Entry(frame_parametros)
    entry_mc_iteraciones.grid(row=0, column=1, padx=5)
    
    tk.Label(frame_parametros, text="Número de qubits:").grid(row=0, column=2, sticky='w')
    entry_qubits = tk.Entry(frame_parametros)
    entry_qubits.grid(row=0, column=3, padx=5)
    
    tk.Label(frame_parametros, text="Beta:").grid(row=1, column=0, sticky='w')
    beta_val = tk.Entry(frame_parametros)
    beta_val.grid(row=1, column=1, padx=5)
    
    tk.Label(frame_parametros, text="Gamma:").grid(row=1, column=2, sticky='w')
    gamma_val = tk.Entry(frame_parametros)
    gamma_val.grid(row=1, column=3, padx=5)
    
    tk.Label(frame_parametros, text="Nivel de profundidad p:").grid(row=1, column=4, sticky='w')
    p_val = tk.Entry(frame_parametros)
    p_val.grid(row=1, column=5, padx=5)

    tk.Label(frame_parametros, text="Tolerancia").grid(row=2, column=2, sticky='w')
    toler_val = tk.Entry(frame_parametros)
    toler_val.grid(row=2, column=3, padx=5)    

    tk.Label(frame_parametros, text="Learning rate").grid(row=2, column=4, sticky='w')
    learning_val = tk.Entry(frame_parametros)
    learning_val.grid(row=2, column=5, padx=5)    
    
    usar_valores_predefinidos = tk.BooleanVar()
    check = tk.Checkbutton(frame_parametros, text="Usar valores predefinidos", variable=usar_valores_predefinidos, command=fijar_valores)
    check.grid(row=2, column=0, columnspan=3, pady=10, sticky='w')
    
    # ==================== SECCIÓN: Rango de búsqueda ====================
    frame_rangos = tk.LabelFrame(ventana, text="Rango de búsqueda", padx=10, pady=10)
    frame_rangos.grid(row=1, column=0, padx=10, pady=5, sticky="nsew")
    
    tk.Label(frame_rangos, text="Límite inferior:").grid(row=0, column=0, sticky='w')
    inf = tk.Entry(frame_rangos)
    inf.grid(row=0, column=1, padx=5)
    
    tk.Label(frame_rangos, text="Límite superior:").grid(row=1, column=0, sticky='w')
    sup = tk.Entry(frame_rangos)
    sup.grid(row=1, column=1, padx=5)
    
    # ==================== SECCIÓN: Coeficientes de la función ====================
    frame_coef = tk.LabelFrame(ventana, text="Coeficientes del polinomio", padx=10, pady=10)
    frame_coef.grid(row=1, column=1, padx=10, pady=5, sticky="nsew", columnspan=2)
    
    V_T4 = tk.Entry(frame_coef, width=7)
    V_T4.grid(row=0, column=0)
    tk.Label(frame_coef, text="x⁴").grid(row=0, column=1)
    
    V_T3 = tk.Entry(frame_coef, width=7)
    V_T3.grid(row=0, column=2)
    tk.Label(frame_coef, text="x³").grid(row=0, column=3)
    
    V_T2 = tk.Entry(frame_coef, width=7)
    V_T2.grid(row=0, column=4)
    tk.Label(frame_coef, text="x²").grid(row=0, column=5)
    
    V_T1 = tk.Entry(frame_coef, width=7)
    V_T1.grid(row=0, column=6)
    tk.Label(frame_coef, text="x").grid(row=0, column=7)
    
    V_T0 = tk.Entry(frame_coef, width=7)
    V_T0.grid(row=0, column=8)
    tk.Label(frame_coef, text="Constante").grid(row=0, column=9)


    # ==================== SECCIÓN: Elección de tipo de función ====================
    frame_opcion = tk.LabelFrame(ventana, text="Tipo de función", padx=10, pady=10)
    frame_opcion.grid(row=2, column=0, padx=10, pady=5, sticky="nsew", columnspan=3)

    opciones = ["Cualquiera", "Lennard-Jones", "Higgs", "Hidrogeno"]
    opcion_var = tk.StringVar(value="Cualquiera")
    menu_combo = ttk.Combobox(frame_opcion, textvariable=opcion_var, values=opciones, state="readonly")
    menu_combo.grid(row=0, column=0, columnspan=2, pady=10)
    menu_combo.current(0)

    menu_combo.bind("<<ComboboxSelected>>", opcion_seleccionada)

    function_def = ttk.Label(frame_opcion, text="Función Polinomio: ax^4 + ... + c")
    function_def.grid(row=0, column=2, padx=10, sticky="w")
        
    # ==================== SECCIÓN: Acciones ====================
    frame_acciones = tk.Frame(ventana)
    frame_acciones.grid(row=3, column=0, columnspan=3, pady=15)
    
    boton_grafica = tk.Button(frame_acciones, text="Graficar función", command=plot_button, width=20)
    boton_grafica.grid(row=0, column=0, padx=10)
    
    boton_ejecutar = tk.Button(frame_acciones, text="Ejecutar QAOA-Monte Carlo", command=ejecutar_MC, width=20)
    boton_ejecutar.grid(row=0, column=1, padx=10)
    
    # ==================== Barra de progreso ====================
    progress = ttk.Progressbar(ventana, orient="horizontal", length=600, mode="determinate")
    progress.grid(row=4, column=0, columnspan=3, pady=10)

    tiempo_label = tk.Label(ventana, text="", font=("Arial", 12), fg="black")
    tiempo_label.grid(row=5, column=0, columnspan=3, pady=10)
    
    resultado_label = tk.Label(ventana, text="", font=("Arial", 12), fg="black")
    resultado_label.grid(row=6, column=0, columnspan=3, pady=10)
    
    resultado_label2 = tk.Label(ventana, text="", font=("Arial", 12), fg="black")
    resultado_label2.grid(row=7, column=0, columnspan=3, pady=10)
    
    
    
    ventana.mainloop()

