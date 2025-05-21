import numpy as np
import autograd as ag
import matplotlib.pyplot as plt
import time

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

import sympy as sp



def ini_vec(inf_range, sup_range):
    return np.random.uniform(inf_range, sup_range)

def newton_method(inf_range, sup_range, n_iter, tolerance, function, x):
    k = ini_vec(inf_range, sup_range)

    
    f_prime = sp.diff(function, x)
    f_second = sp.diff(function, x, 2)

    for i in range(n_iter):
        if f_second.subs(x, k) == 0:
            print("Error en la segunda derivada, division entre 0")
            break

        f_prime_val = f_prime.subs(x, k)
        f_second_val = f_second.subs(x, k)
        
        frac = float(f_prime_val / f_second_val)
        k_new = (k - frac)

        if abs(k_new - k) < tolerance:
            break
            
        k = k_new

    return k

#====================FUNCIONES PARA LA INTERFAZ ===========================#

############################ BOTON DE PLOTEO ####################################

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

############################### BOTON DE EJECUCION ################################

def ejecutar_newton_method_MC():
    global T4, T3, T2, T1, T0
      
    n_iter = int(entry_iteraciones.get())
    tolerance = float(entry_tolerance.get())
    num_mc = int(entry_mc.get())
    min_dist = float(entry_min_dist.get())
    
    soluciones = []

    x = sp.symbols('x')
    opcion = opcion_var.get()
    
    if opcion == "Cualquiera":
        def f(x):
            y = T4*(x**4) + T3*(x**3) +T2*(x**2) + T1*x + T0
            return -y

        function = T4*(x**4) + T3*(x**3) +T2*(x**2) + T1*x + T0
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

        function = y = 4 * 10 * ((0.25/x)**12 - (0.25/x)**6)
        inf_range = 0.8
        sup_range = 3

    elif opcion == "Higgs":
        def f(x):
            y = -1 * x**2 + 0.5 * x**4
            return -y

        function = -1 * x**2 + 0.5 * x**4
        inf_range = -2
        sup_range = 2

    elif opcion == "Hidrogeno":
        def f(x):
            y = 3*(x**4) - 7*(x**2) - 0.5*x + 6
            return -y

        function = 3*(x**4) - 7*(x**2) - 0.5*x + 6
        inf_range = -2
        sup_range = 2

    
    start_time = time.time()    

    progress["value"] = 0
    progress["maximum"] = num_mc
    for n in range(num_mc):
        minimo = newton_method(inf_range, sup_range, n_iter, tolerance, function, x)
        progress["value"] += 1
        ventana.update_idletasks()
        print(n)
        # Evitar duplicados (comparando con tolerancia)
        if all(abs(minimo - s) > min_dist for s in soluciones):
            soluciones.append(minimo)
            if len(soluciones) > 2:
                progress["value"] = num_mc
                ventana.update_idletasks()
                break


    elapsed_time = time.time() - start_time
    tiempo_label.config(text=f"Tiempo de ejecucion: {elapsed_time:.2f} segundos")

    soluciones.sort()
    resultado_label.config(text=f"Se encontraron {len(soluciones)} mínimos:\n" + "\n".join([f"{s:.5f}" for s in soluciones]))

######################## FUNCIONES DEL MENU DESPLEGABLE ###########################

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
    ventana.title("Newton Method")
    ventana.geometry("1600x800")

    # ==================== SECCIÓN: Parámetros del Algoritmo ====================
    frame_parametros = tk.LabelFrame(ventana, text="Parámetros del algoritmo", padx=10, pady=10)
    frame_parametros.grid(row=0, column=0, padx=10, pady=10, sticky="nsew", columnspan=3)
    
    tk.Label(frame_parametros, text="Número de iteraciones:").grid(row=0, column=0, sticky='w')
    entry_iteraciones = tk.Entry(frame_parametros)
    entry_iteraciones.grid(row=0, column=1, padx=5)
    
    tk.Label(frame_parametros, text="Tolerancia:").grid(row=0, column=2, sticky='w')
    entry_tolerance = tk.Entry(frame_parametros)
    entry_tolerance.grid(row=0, column=3, padx=5)

    tk.Label(frame_parametros, text="Numero Montecarlo:").grid(row=0, column=4, sticky='w')
    entry_mc = tk.Entry(frame_parametros)
    entry_mc.grid(row=0, column=5, padx=5)

    tk.Label(frame_parametros, text="Distancia entre minimos:").grid(row=1, column=0, sticky='w')
    entry_min_dist = tk.Entry(frame_parametros)
    entry_min_dist.grid(row=1, column=1, padx=5)

    
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
    opcion_menu = ttk.Combobox(frame_opcion, textvariable=opcion_var, values=opciones, state="readonly")
    opcion_menu.grid(row=0, column=0, padx=10, pady=5)
    opcion_menu.bind("<<ComboboxSelected>>", opcion_seleccionada)

    function_def = ttk.Label(frame_opcion, text="Función Polinomio: ax^4 + ... + c")
    function_def.grid(row=0, column=2, padx=10, sticky="w")

    # ==================== SECCIÓN: Acciones ====================
    frame_acciones = tk.Frame(ventana)
    frame_acciones.grid(row=3, column=0, columnspan=3, pady=15)
    
    boton_grafica = tk.Button(frame_acciones, text="Graficar función", command=plot_button, width=20)
    boton_grafica.grid(row=0, column=0, padx=10)
    
    boton_ejecutar = tk.Button(frame_acciones, text="Ejecutar Newton Method", command=ejecutar_newton_method_MC, width=20)
    boton_ejecutar.grid(row=0, column=1, padx=10)

    # ==================== Barra de progreso ====================
    progress = ttk.Progressbar(ventana, orient="horizontal", length=600, mode="determinate")
    progress.grid(row=4, column=0, columnspan=3, pady=10)

    tiempo_label = tk.Label(ventana, text="", font=("Arial", 12), fg="black")
    tiempo_label.grid(row=5, column=0, columnspan=3, pady=10)
    
    resultado_label = tk.Label(ventana, text="", font=("Arial", 12), fg="black")
    resultado_label.grid(row=6, column=0, columnspan=3, pady=10)

    ventana.mainloop()
    