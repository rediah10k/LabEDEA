import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox

def euler(a, b, c, d, X0, t_fin, dt, num_trayectorias):
    """
    Método de Euler para resolver una ecuación estocástica.
    Coeficiente 'a' de la ecuación.
    Coeficiente 'b' de la ecuación.
    Coeficiente 'c' de la ecuación.
    Coeficiente 'd' de la ecuación.
    X0: Condición inicial.
    t_fin: Tiempo final.
    dt: Paso de tiempo.
    num_trayectorias: Número de trayectorias a graficar.
    """
    npuntos = int(t_fin / dt)
    sqrtdt = np.sqrt(dt)

    X = np.zeros((num_trayectorias, npuntos + 1))
    X[:, 0] = X0

    # Planteamiento de la EDE con los valores parametrizados
    for t in range(npuntos):
        dBt = sqrtdt * np.random.normal(0, 1, num_trayectorias)
        X[:, t] = X[:, t-1] + (a * X[:, t] + b) * dt + (c * X[:, t] + d) * dBt

    t = np.linspace(0, t_fin, npuntos + 1)
    
    term1 = (2 * b * c * d) / (a * (a - c**2))
    term2 = (2 * c * d * X0) / (a - c**2)
    term3 = (2 * b * d) / (a * c)
    term4 = (d**2) / (c**2)
    
    # Cálculo del valor esperado teórico
    E_X_teorico = (b / a) * (np.exp(a * t) - 1) + X0 * np.exp(a * t)
    E_X2_teorico = (np.exp(a * t) * (term1 + term2)) + term3 - term4 + np.exp(c**2 * t)
    Var_teorica = E_X2_teorico - E_X_teorico**2
    Var_teorica=np.sqrt(Var_teorica)
    # Graficar trayectorias simuladas
    for i in range(num_trayectorias):
        plt.plot(t, X[i, :], color='orange', linewidth=0.5)

    #Var_teorico = E_X2_teorico - np.sqrt(E_X_teorico)
    # Graficar valor esperado teórico y bandas de varianza teórica
    plt.plot(t, E_X_teorico, color='red', linewidth=2, label='Valor esperado')
    plt.plot(t, Var_teorica, color='blue', linewidth=2, label='Varianza del proceso')
    
    #ESPACIO PARA GRAFICAR LA VARIANZA

    plt.xlabel('Tiempo')
    plt.ylabel('$X_t$')
    plt.title('Trayectorias Simuladas y Valor Esperado Teórico')
    plt.legend()
    plt.grid(True)
    plt.show()


def solicitar_parametros_iniciales():
    ventana_parametros = tk.Tk()
    ventana_parametros.title("Parámetros Iniciales")
    
    labels = ['a', 'b', 'c', 'd', 'X0', 't_fin', 'dt', 'num_trayectorias']
    entradas = {}

    for i, label in enumerate(labels):
        tk.Label(ventana_parametros, text=label).grid(row=i, column=0, pady=5)
        entrada = tk.Entry(ventana_parametros)
        entrada.grid(row=i, column=1, pady=5, padx=20)
        entradas[label] = entrada

    def guardar_parametros():
        global a, b, c, d, X0, t_fin, dt, num_trayectorias
        try:
            a = float(entradas['a'].get())
            b = float(entradas['b'].get())
            c = float(entradas['c'].get())
            d = float(entradas['d'].get())
            X0 = float(entradas['X0'].get())
            t_fin = float(entradas['t_fin'].get())
            dt = float(entradas['dt'].get())
            num_trayectorias = int(entradas['num_trayectorias'].get())

            # Validación de que ningún campo esté vacío
            if any(not entrada.get() for entrada in entradas.values()):
                raise ValueError("Todos los campos deben estar llenos")
            
            if t_fin <= 0:
                raise ValueError("El tiempo final no puede ser cero")
            
            if dt < 0 or dt >= 1:
                raise ValueError("El paso del tiempo debe ser un valor pequeño entre 0 y 1")

            if num_trayectorias <= 0:
                raise ValueError("El número de trayectorias debe ser mayor que cero")

            messagebox.showinfo("Éxito", "Parámetros guardados con éxito")
            ventana_parametros.destroy()
            menu()

        except ValueError:
            messagebox.showerror("Error", "Por favor, ingrese todos los parámetros correctamente.")

    tk.Button(ventana_parametros, text="Guardar", command=guardar_parametros).grid(row=len(labels), columnspan=2, pady=10)

    ventana_parametros.mainloop()

def menu():
    ventana = tk.Tk()
    ventana.title("EDE autónoma")
    ventana.geometry("500x300")

    def cerrar_programa():
        messagebox.showinfo("Adios!", "Has decidido salir del programa.")
        ventana.destroy() 
    
    def alterar_parametros():
        parametros_ventana = tk.Toplevel(ventana)
        parametros_ventana.title("Alterar Parámetros")
        
        labels = ['a', 'b', 'c', 'd', 'X0', 't_fin', 'dt', 'num_trayectorias']
        entradas = {}

        for i, label in enumerate(labels):
            tk.Label(parametros_ventana, text=label).grid(row=i, column=0, pady=5)
            entrada = tk.Entry(parametros_ventana)
            entrada.grid(row=i, column=1, pady=5, padx=20)
            entradas[label] = entrada

        def guardar_parametros():
            global a, b, c, d, X0, t_fin, dt, num_trayectorias
            try:
                a = float(entradas['a'].get())
                b = float(entradas['b'].get())
                c = float(entradas['c'].get())
                d = float(entradas['d'].get())
                X0 = float(entradas['X0'].get())
                t_fin = float(entradas['t_fin'].get())
                dt = float(entradas['dt'].get())
                num_trayectorias = int(entradas['num_trayectorias'].get())

                # Validación de que ningún campo esté vacío
                if any(not entrada.get() for entrada in entradas.values()):
                    raise ValueError("Todos los campos deben estar llenos")
                
                if t_fin <= 0:
                    raise ValueError("El tiempo final no puede ser cero")
            
                if dt < 0 or dt >= 1:
                    raise ValueError("El paso del tiempo debe ser un valor pequeño entre 0 y 1")

                if num_trayectorias <= 0:
                    raise ValueError("El número de trayectorias debe ser mayor que cero")

                messagebox.showinfo("Éxito", "Parámetros actualizados con éxito")
                parametros_ventana.destroy()

            except ValueError:
                messagebox.showerror("Error", "Por favor, ingrese todos los parámetros correctamente.")

        tk.Button(parametros_ventana, text="Guardar", command=guardar_parametros).grid(row=len(labels), columnspan=2, pady=10)

    def graficar():
        euler(a, b, c, d, X0, t_fin, dt, num_trayectorias)

    label_menu = tk.Label(ventana, text="Menú", font=("Arial", 16, "bold"))
    label_menu.pack(pady=(10, 0))

    label_description = tk.Label(ventana, text="Seleccione una opción.", font=("Arial", 11))
    label_description.pack(pady=(10,20))

    boton_cambiar_parametros = tk.Button(ventana, font=('Arial', 10), text="1. Modificar parametros de la aplicacion", command=alterar_parametros)
    boton_cambiar_parametros.pack(pady=10)

    boton_graficar = tk.Button(ventana, font=('Arial', 10), text="2. Graficar Trayectorias", command=graficar)
    boton_graficar.pack(pady=10)

    boton_salir = tk.Button(ventana, font=('Arial', 10), text="3. Salir", command=cerrar_programa)
    boton_salir.pack(pady=10)

    ventana.mainloop()

def iniciar_app():
    solicitar_parametros_iniciales()

iniciar_app()