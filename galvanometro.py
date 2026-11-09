
# ==========================================
# SIMULADOR DE GALVANOMETRO
# PROYECTO DE ECUACIONES DIFERENCIALES
# ==========================================

# Importar librerias
import math
import matplotlib.pyplot as plt


# ==========================================
# 1. PARAMETROS DEL GALVANOMETRO
# ==========================================

# Momento de inercia
J = 0.000001

# Coeficiente de amortiguamiento
B = 0.000014

# Constante del resorte
KAPPA = 0.0001

# Constante del torque electromagnetico
G = 0.01


# ==========================================
# 2. CORRIENTE ELECTRICA
# ==========================================

# Corriente aplicada al galvanometro
I0 = 0.005


# ==========================================
# 3. CONDICIONES INICIALES
# ==========================================

# Angulo inicial de la aguja
theta = 0.0

# Velocidad angular inicial
velocidad = 0.0


# ==========================================
# 4. CONFIGURACION DE LA SIMULACION
# ==========================================

# Tiempo total de simulacion
tiempo_final = 2.0

# Paso de tiempo
dt = 0.001


# ==========================================
# 5. FUNCION DE CORRIENTE ELECTRICA
# ==========================================

def corriente(t):

    return I0


# ==========================================
# 6. FUNCION DE LA ECUACION DIFERENCIAL
# ==========================================

def aceleracion(t, theta, velocidad):

    # Torque producido por la corriente
    fuerza = G * corriente(t)

    # Calcular aceleracion angular
    resultado = (
        fuerza - B * velocidad - KAPPA * theta
    ) / J

    return resultado


# ==========================================
# 7. METODO RUNGE-KUTTA DE CUARTO ORDEN
# ==========================================

def runge_kutta(theta, velocidad, t, dt):

    # PRIMER CALCULO
    k1_theta = velocidad

    k1_velocidad = aceleracion(
        t, theta, velocidad
    )


    # SEGUNDO CALCULO
    k2_theta = velocidad + (
        dt / 2
    ) * k1_velocidad

    k2_velocidad = aceleracion(
        t + dt / 2,
        theta + (dt / 2) * k1_theta,
        velocidad + (dt / 2) * k1_velocidad
    )


    # TERCER CALCULO
    k3_theta = velocidad + (
        dt / 2
    ) * k2_velocidad

    k3_velocidad = aceleracion(
        t + dt / 2,
        theta + (dt / 2) * k2_theta,
        velocidad + (dt / 2) * k2_velocidad
    )


    # CUARTO CALCULO
    k4_theta = velocidad + dt * k3_velocidad

    k4_velocidad = aceleracion(
        t + dt,
        theta + dt * k3_theta,
        velocidad + dt * k3_velocidad
    )


    # ACTUALIZAR EL ANGULO
    nuevo_theta = theta + (dt / 6) * (
        k1_theta +
        2 * k2_theta +
        2 * k3_theta +
        k4_theta
    )


    # ACTUALIZAR LA VELOCIDAD
    nueva_velocidad = velocidad + (dt / 6) * (
        k1_velocidad +
        2 * k2_velocidad +
        2 * k3_velocidad +
        k4_velocidad
    )


    # RETORNAR LOS NUEVOS VALORES
    return nuevo_theta, nueva_velocidad


# ==========================================
# 8. SIMULACION DEL GALVANOMETRO
# ==========================================

# Listas para guardar los resultados
tiempos = [0]

angulos = [theta]

velocidades = [velocidad]


# Calcular numero de pasos
n_pasos = int(tiempo_final / dt)


# Repetir los calculos
for i in range(n_pasos):

    # Tiempo actual
    t = i * dt

    # Aplicar el metodo Runge-Kutta
    theta, velocidad = runge_kutta(
        theta,
        velocidad,
        t,
        dt
    )

    # Guardar los resultados
    tiempos.append((i + 1) * dt)

    angulos.append(theta)

    velocidades.append(velocidad)

# ==========================================
# 9. CALCULAR CARACTERISTICAS DEL SISTEMA
# ==========================================

# Frecuencia natural
omega_n = math.sqrt(KAPPA / J)


# Factor de amortiguamiento
zeta = B / (2 * math.sqrt(KAPPA * J))


# Angulo final teorico
theta_final = (G * I0) / KAPPA

# ==========================================
# 10. SOLUCION ANALITICA
# ==========================================

angulos_analiticos = []


# Verificar si el sistema es subamortiguado
if zeta < 1:

    # Frecuencia amortiguada
    omega_d = omega_n * math.sqrt(1 - zeta**2)


    # Calcular la solucion analitica
    for t in tiempos:

        valor = theta_final * (
            1 - math.exp(-zeta * omega_n * t) * (
                math.cos(omega_d * t) +
                (zeta / math.sqrt(1 - zeta**2)) *
                math.sin(omega_d * t)
            )
        )

        angulos_analiticos.append(valor)


    # Calcular el error maximo
    error_maximo = max(
        abs(a - b)
        for a, b in zip(angulos, angulos_analiticos)
    )

else:

    error_maximo = None

# ==========================================
# 11. MOSTRAR RESULTADOS
# ==========================================

print("\n===================================")
print("      SIMULADOR DE GALVANOMETRO")
print("===================================")

print("\nPARAMETROS DEL SISTEMA")

print("Momento de inercia:", J)
print("Amortiguamiento:", B)
print("Constante del resorte:", KAPPA)
print("Constante del torque:", G)

print("\nCORRIENTE ELECTRICA")

print("Corriente aplicada:", I0, "A")

print("\nRESULTADOS MATEMATICOS")

print("Frecuencia natural:", omega_n, "rad/s")

print("Factor de amortiguamiento:", zeta)

print("Angulo final teorico:", theta_final, "rad")

print("Angulo final numerico:", angulos[-1], "rad")


print("\nVALIDACION NUMERICA")

if error_maximo is not None:

    print("Error maximo:", error_maximo, "rad")

else:

    print("No se calculo la comparacion analitica.")

print("\n===================================")

# ==========================================
# 12. GRAFICA DEL MOVIMIENTO DE LA AGUJA
# ==========================================

plt.figure(figsize=(10, 5))


# Graficar solucion numerica
plt.plot(
    tiempos,
    angulos,
    label="Solucion numerica (RK4)"
)


# Graficar solucion analitica
if angulos_analiticos:

    plt.plot(
        tiempos,
        angulos_analiticos,
        "--",
        label="Solucion analitica"
    )


# Mostrar angulo final
plt.axhline(
    y=theta_final,
    color="gray",
    linestyle=":",
    label="Angulo final"
)


# Titulos y etiquetas
plt.title("Simulacion del movimiento del galvanometro")

plt.xlabel("Tiempo (segundos)")

plt.ylabel("Angulo de la aguja (radianes)")


# Mostrar leyenda y cuadricula
plt.legend()

plt.grid()


# Guardar la grafica
plt.savefig("grafica_galvanometro.png")


# Mostrar la grafica
plt.show()
