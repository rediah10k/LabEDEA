import pygame
import sys
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
import cmath

# Variables globales
m_value = c_value = k_value = sigma_value = x0_value = v0_value = dt_value = tf_value = tr_value = None

def validate_entries():
    entries = [entry_m, entry_c, entry_k, entry_sigma, entry_x0, entry_v0, entry_dt, entry_tf, entry_tr]
    for entry in entries:
        try:
            val = float(entry.get())
            if val < 0:
                return False
        except ValueError:
            return False
    return True

def retrieve_input():
    global m_value, c_value, k_value, sigma_value, x0_value, v0_value, dt_value, tf_value, tr_value
    if not validate_entries():
        messagebox.showerror("Error", "Todos los campos deben tener valores numéricos no negativos.")
        return

    m_value = float(entry_m.get())
    c_value = float(entry_c.get())
    k_value = float(entry_k.get())
    sigma_value = float(entry_sigma.get())
    x0_value = float(entry_x0.get())
    v0_value = float(entry_v0.get())
    dt_value = float(entry_dt.get())
    tf_value = float(entry_tf.get())
    tr_value = int(entry_tr.get())
    
    x_tiempoPosicion, x_esperado, x_tiempo, v_esperado = hacer_graficas()
    iniciar_animacion(x_tiempoPosicion, x_esperado, x_tiempo, v_esperado)

def hacer_graficas():
    def funcionPaso(x, v):
        return -k_value * x / m_value - c_value * v / m_value

    def metodo_euler(x0, v0, f, dt, tf):
        xTiempo = np.arange(0, tf + dt, dt)
        xValue = np.zeros(len(xTiempo))
        xValue[0] = x0
        vValue = np.zeros(len(xTiempo))
        vValue[0] = v0
        for i in range(1, len(xTiempo)):
            xValue[i] = xValue[i - 1] + vValue[i - 1] * dt
            vValue[i] = vValue[i - 1] + funcionPaso(xValue[i - 1], vValue[i - 1]) * dt + (sigma_value / m_value) * np.sqrt(dt) * np.random.normal(0, 1)
        return xTiempo, xValue, vValue

    # Generar las trayectorias
    trayectorias = [metodo_euler(x0_value, v0_value, funcionPaso, dt_value, tf_value) for _ in range(tr_value)]
    
    # Generación de gráficas con plano cartesiano y grilla de una unidad
    plt.figure(figsize=(10, 5))
    
    def valorMedioPosicion(x0, v0, c, k, m, tf, dt):
        x_tiempo = np.arange(0, tf + dt, dt)
        x_esperado = np.zeros(len(x_tiempo))
        gamma = c / m
        omega_squared = - k / m + (gamma ** 2) / 4
        
        omega = cmath.sqrt(omega_squared)
        for i in range(len(x_tiempo)):
            t = x_tiempo[i]
            exp_term = cmath.exp(-gamma * t / 2)
            cos_term = cmath.cosh(omega * t)
            sin_term = cmath.sinh(omega * t)
            term1 = (((m*v0+c*x0)/(m*x0)) - (gamma / 2)) / omega
            x_esperado[i] = float((exp_term * x0 * ( cos_term + term1 * sin_term)).real)
    
        return x_tiempo, x_esperado
    
    x_tiempoPosicion, x_esperado = valorMedioPosicion(x0_value, v0_value, c_value, k_value, m_value, tf_value, dt_value)
    
    def valorMedioVelocidad(x0, v0, c, k, m, tf, dt):
        x_tiempo = np.arange(0, tf + dt, dt)
        v_esperado = np.zeros(len(x_tiempo))
        gamma = c / m
        omega_squared = -k / m + (gamma ** 2) / 4
        
        omega = cmath.sqrt(omega_squared)
        for i in range(len(x_tiempo)):
            t = x_tiempo[i]
            exp_term = cmath.exp(-gamma * t / 2)
            cos_term = cmath.cosh(omega * t)
            sin_term = cmath.sinh(omega * t)
            term1 = (((k*x0)/(m*v0)) + (gamma / 2)) / omega
            v_esperadoFloat = float((exp_term * v0 * ( cos_term - term1 * sin_term)).real)
            v_esperado[i] = v_esperadoFloat
        
        return x_tiempo, v_esperado

    x_tiempo, v_esperado = valorMedioVelocidad(x0_value, v0_value, c_value, k_value, m_value, tf_value, dt_value)
    
    # Trayectorias de posición con plano cartesiano
    plt.subplot(2, 1, 1)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    for xTiempo, xValue, vValue in trayectorias:
        plt.plot(xTiempo, xValue, color="red")
        
    plt.plot(x_tiempoPosicion, x_esperado, label='Posición Esperada', color='black', lw=2)
    plt.xlabel('Tiempo')
    plt.ylabel('Posición')
    plt.title('Trayectorias de Posición')
    plt.legend()

    # Trayectorias de velocidad con plano cartesiano
    plt.subplot(2, 1, 2)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    for xTiempo, xValue, vValue in trayectorias:
        plt.plot(xTiempo, vValue, color="blue")
    plt.plot(x_tiempo, v_esperado, label='Velocidad Esperada', color='black', lw=2)
    plt.xlabel('Tiempo')
    plt.ylabel('Velocidad')
    plt.title('Trayectorias de Velocidad')
    plt.legend()
    plt.tight_layout()
    
    plt.show()
    
    return x_tiempoPosicion, x_esperado, x_tiempo, v_esperado

def iniciar_animacion(x_tiempoPosicion, x_esperado, x_tiempo, v_esperado):
    pygame.init()
    pygame.font.init()

    # Obtener el tamaño de la pantalla
    info_obj = pygame.display.Info()
    anchoP, altoP = info_obj.current_w, info_obj.current_h

    # Configuración del eje de anclaje
    margen_superior = 50  # píxeles desde el borde superior
    margen_inferior = 50  # píxeles desde el borde inferior
    espacio_disponible = altoP - margen_superior - margen_inferior

    # Ajustar el escalado para que siempre se vea el sistema completo
    escala = espacio_disponible / 80  # Factor de escala para ajustar los valores de posición y velocidad a la pantalla de Pygame
    desplazamiento = margen_superior + espacio_disponible / 2  # Desplazamiento para mantener la masa en el centro

    freq_vibracion = 30  # Frecuencia de vibración del punto de anclaje en Hz
    amplitud_vibracion = 2  # Amplitud de vibración del punto de anclaje en píxeles

    # Pantalla
    screen = pygame.display.set_mode((anchoP, altoP), pygame.FULLSCREEN)
    clock = pygame.time.Clock()

    # Secciones
    animacion = pygame.Surface((anchoP, altoP))
    fuente = pygame.font.Font(None, 20)

    # Colores pastel
    fondo_color = (248, 248, 255)
    barra_color = (219, 112, 147)
    resorte_color = (32, 178, 170)
    masa_color = (240, 230, 140)
    texto_color = (105, 105, 105)

    tiempo_inicial = pygame.time.get_ticks()

    y_anim = desplazamiento + x0_value * escala
    v_anim = v0_value

    x_esperado_list = []
    v_esperado_list = []

    # Calcular valores esperados de posición y velocidad para la animación
    for i in range(int(tf_value / dt_value)):
        tiempo_actual = i * dt_value
        x_esperado_list.append(np.interp(tiempo_actual, x_tiempoPosicion, x_esperado))
        v_esperado_list.append(np.interp(tiempo_actual, x_tiempo, v_esperado))

    while True:
        clock.tick(60)  # Framerate
        
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()

        animacion.fill(fondo_color)  # colorea la pantalla de blanco
        
        # Calcular el tiempo transcurrido
        tiempo_actual = pygame.time.get_ticks()
        Tiempo = (tiempo_actual - tiempo_inicial) / 1000  # Tiempo en segundos

        if Tiempo > tf_value:
            break

        # Posición vibratoria del punto de anclaje
        y_anclaje_vibracion = margen_superior + amplitud_vibracion * np.sin(2 * np.pi * freq_vibracion * Tiempo)

        # Obtener la posición y velocidad esperada en el tiempo actual
        indice_actual = int(Tiempo / dt_value)
        if indice_actual < len(x_esperado_list):
            y_esperado = desplazamiento + x_esperado_list[indice_actual] * escala
            v_esperado = v_esperado_list[indice_actual]

            # Actualizar la posición y velocidad de la masa según los valores esperados
            y_anim = y_esperado
            v_anim = v_esperado

        # Asegurarse de que la masa no supere la posición del anclaje
        if y_anim < y_anclaje_vibracion:
            y_anim = y_anclaje_vibracion
            v_anim *= -0.5  # Rebote con pérdida de energía

        # Dibuja la barra del punto de anclaje
        barra_rect = pygame.Rect(anchoP // 2 - 40, y_anclaje_vibracion, 80, 10)
        pygame.draw.rect(animacion, barra_color, barra_rect)

        # Dibuja el resorte
        pygame.draw.line(animacion, resorte_color, (anchoP // 2, y_anclaje_vibracion + 10), (anchoP // 2, y_anim), 5)

        # Dibuja la masa
        masa_rect = pygame.Rect(anchoP // 2 - 20, y_anim, 40, 40)
        pygame.draw.rect(animacion, masa_color, masa_rect)

        # Dibujar el peso en la masa
        peso_texto = fuente.render(f'{m_value}', True, texto_color)
        texto_rect = peso_texto.get_rect(center=masa_rect.center)
        animacion.blit(peso_texto, texto_rect)

        # Dibujar la posición esperada al lado de la masa
        posicion_esperada_texto = fuente.render(f'Posicion trayectoria: {(x_esperado_list[indice_actual]):.2f}', True, texto_color)
        animacion.blit(posicion_esperada_texto, (anchoP // 2 + 50, y_anim))

        # Dibuja la sección de datos en tiempo real
        tiempo_texto = fuente.render(f'Tiempo: {Tiempo:.2f} s', True, texto_color)
        velocidad_texto = fuente.render(f'Velocidad: {v_anim:.2f} m/s', True, texto_color)
        posicion_texto = fuente.render(f'Posición: {(y_anim - desplazamiento) / escala:.2f} m', True, texto_color)

        animacion.blit(tiempo_texto, (10, 10))
        animacion.blit(velocidad_texto, (10, 30))
        animacion.blit(posicion_texto, (10, 50))

        # Dibuja el plano cartesiano
        pygame.draw.line(animacion, (0, 0, 0), (0, desplazamiento), (anchoP, desplazamiento), 1)
        pygame.draw.line(animacion, (0, 0, 0), (anchoP // 2, 0), (anchoP // 2, altoP), 1)
        for i in range(0, anchoP, 50):
            pygame.draw.line(animacion, (200, 200, 200), (i, 0), (i, altoP), 1)
        for j in range(0, altoP, 50):
            pygame.draw.line(animacion, (200, 200, 200), (0, j), (anchoP, j), 1)

        screen.blit(animacion, (0, 0))

        pygame.display.update()

    pygame.quit()

# Interfaz gráfica con Tkinter
window = tk.Tk()
window.title("Simulación de Sistema de Masa-Resorte-Amortiguador")

# Aplicar estilo ttk
style = ttk.Style()
style.theme_use('clam')

# Crear frame para contener los widgets
frame = ttk.Frame(window, padding="10 10 10 10")
frame.pack(fill=tk.BOTH, expand=True)

# Crear etiquetas y campos de entrada
ttk.Label(frame, text="Masa (m):").grid(column=0, row=0, sticky=tk.W)
entry_m = ttk.Entry(frame)
entry_m.grid(column=1, row=0, sticky=(tk.W, tk.E))
    
ttk.Label(frame, text="Constante de amortiguamiento (c):").grid(column=0, row=1, sticky=tk.W)
entry_c = ttk.Entry(frame)
entry_c.grid(column=1, row=1, sticky=(tk.W, tk.E))

ttk.Label(frame, text="Constante del resorte (k):").grid(column=0, row=2, sticky=tk.W)
entry_k = ttk.Entry(frame)
entry_k.grid(column=1, row=2, sticky=(tk.W, tk.E))

ttk.Label(frame, text="Sigma:").grid(column=0, row=3, sticky=tk.W)
entry_sigma = ttk.Entry(frame)
entry_sigma.grid(column=1, row=3, sticky=(tk.W, tk.E))

ttk.Label(frame, text="Posición Inicial (x0):").grid(column=0, row=4, sticky=tk.W)
entry_x0 = ttk.Entry(frame)
entry_x0.grid(column=1, row=4, sticky=(tk.W, tk.E))

ttk.Label(frame, text="Velocidad Inicial (v0):").grid(column=0, row=5, sticky=tk.W)
entry_v0 = ttk.Entry(frame)
entry_v0.grid(column=1, row=5, sticky=(tk.W, tk.E))

ttk.Label(frame, text="Paso (dt):").grid(column=0, row=6, sticky=tk.W)
entry_dt = ttk.Entry(frame)
entry_dt.grid(column=1, row=6, sticky=(tk.W, tk.E))

ttk.Label(frame, text="Número de Trayectorias:").grid(column=0, row=7, sticky=tk.W)
entry_tr = ttk.Entry(frame)
entry_tr.grid(column=1, row=7, sticky=(tk.W, tk.E))

ttk.Label(frame, text="Tiempo Final (tf):").grid(column=0, row=8, sticky=tk.W)
entry_tf = ttk.Entry(frame)
entry_tf.grid(column=1, row=8, sticky=(tk.W, tk.E))

button = tk.Button(window, text="Simular", command=retrieve_input)
button.pack()

exit_button = tk.Button(window, text="Salir", command=window.quit)
exit_button.pack(side=tk.BOTTOM)

# Ajustar el espaciado de las columnas
for child in frame.winfo_children():
    child.grid_configure(padx=5, pady=5)

# Iniciar el bucle principal de Tkinter
window.mainloop()
