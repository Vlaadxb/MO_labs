import requests
import numpy as np 
import matplotlib.pyplot as plt

url = "https://api.open-elevation.com/api/v1/lookup?locations=48.164214,24.536044|48.164983,24.534836|48.165605,24.534068|48.166228,24.532915|48.166777,24.531927|48.167326,24.530884|48.167011,24.530061|48.166053,24.528039|48.166655,24.526064|48.166497,24.523574|48.166128,24.520214|48.165416,24.517170|48.164546,24.514640|48.163412,24.512980|48.162331,24.511715|48.162015,24.509462|48.162147,24.506932|48.161751,24.504244|48.161197,24.501793|48.160580,24.500537|48.160250,24.500106"

response = requests.get(url)
data = response.json()

results = data["results"]

n = len(results)
print("Кількість вузлів:", n)

print("\nТабуляція вузлів:")
print("№ | Latitude | Longitude | Elevation (m)")

with open("results.txt", "w") as file:

    file.write("№ | Latitude | Longitude | Elevation (m)\n")

    for i, point in enumerate(results):
        line = f"{i:2d} | {point['latitude']:.6f} | {point['longitude']:.6f} | {point['elevation']:.2f}"
        
        print(line)
        file.write(line + "\n")
        
import numpy as np

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000  # радіус Землі в метрах

    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)

    a = np.sin(dphi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2) ** 2
    return 2 * R * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

coords = [(p["latitude"], p["longitude"]) for p in results]
elevations = [p["elevation"] for p in results]

distances = [0]

for i in range(1, n):
    d = haversine(*coords[i - 1], *coords[i])
    distances.append(distances[-1] + d)

print("\nТабуляція (відстань, висота):")
print("№ | Distance (m) | Elevation (m)")

for i in range(n):
    print(f"{i:2d} | {distances[i]:10.2f} | {elevations[i]:8.2f}")
    
plt.figure()

plt.plot(distances, elevations)

plt.title("Профіль маршруту: Заросляк — Говерла")
plt.xlabel("Кумулятивна відстань (м)")
plt.ylabel("Висота (м)")

plt.grid()

plt.show()

from scipy.interpolate import CubicSpline

cs = CubicSpline(distances, elevations)

coeffs = cs.c

print("\nКоефіцієнти сплайнів (a, b, c, d):\n")

for i in range(len(distances) - 1):
    a = coeffs[3][i]
    b = coeffs[2][i]
    c = coeffs[1][i]
    d = coeffs[0][i]

    print(f"Інтервал {i}:")
    print(f"a = {a:.4f}, b = {b:.4f}, c = {c:.4f}, d = {d:.4f}\n")
    
# ---------------------------------

def thomas_algorithm(a, b, c, d):
    n = len(d)
    
    a = [float(x) for x in a]
    b = [float(x) for x in b]
    c = [float(x) for x in c]
    d = [float(x) for x in d]

    for i in range(1, n):
        m = a[i] / b[i - 1]
        b[i] = b[i] - m * c[i - 1]
        d[i] = d[i] - m * d[i - 1]

    x = [0.0] * n
    x[-1] = d[-1] / b[-1]

    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]

    return x

print("\nПеревірка методу прогонки:")

a_test = [0, 1, 1]
b_test = [4, 4, 4]
c_test = [1, 1, 0]
d_test = [7, 8, 7]

solution = thomas_algorithm(a_test, b_test, c_test, d_test)

print("Розв'язок тестової системи:")
for i, value in enumerate(solution):
    print(f"x{i + 1} = {value:.4f}")
    
# ---------------------------------

x = distances
y = elevations

n_points = len(x)
h = [x[i + 1] - x[i] for i in range(n_points - 1)]

if n_points < 3:
    print("\nНедостатньо точок для обчислення коефіцієнтів c_i")
else:
    
    m = n_points - 2

    a = [0.0] * m   
    b = [0.0] * m   
    c = [0.0] * m   
    d = [0.0] * m   

    for i in range(m):
        a[i] = h[i] if i > 0 else 0.0
        b[i] = 2 * (h[i] + h[i + 1])
        c[i] = h[i + 1] if i < m - 1 else 0.0
        d[i] = 3 * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i])

    c_inner = thomas_algorithm(a, b, c, d)

    c_coeffs = [0.0] + c_inner + [0.0]

    print("\nКоефіцієнти c_i кубічних сплайнів:")
    for i, value in enumerate(c_coeffs):
        print(f"c[{i}] = {value:.6f}")
        
# ---------------------------------

a_coeffs = []
b_coeffs = []
d_coeffs = []

for i in range(n_points - 1):
    a_i = y[i]
    b_i = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * c_coeffs[i] + c_coeffs[i + 1]) / 3
    d_i = (c_coeffs[i + 1] - c_coeffs[i]) / (3 * h[i])

    a_coeffs.append(a_i)
    b_coeffs.append(b_i)
    d_coeffs.append(d_i)

print("\nКоефіцієнти кубічних сплайнів:")

for i in range(n_points - 1):
    print(f"\nІнтервал {i}:")
    print(f"a[{i}] = {a_coeffs[i]:.6f}")
    print(f"b[{i}] = {b_coeffs[i]:.6f}")
    print(f"c[{i}] = {c_coeffs[i]:.6f}")
    print(f"d[{i}] = {d_coeffs[i]:.6f}")
    
# ---------------------------------

# ---------------------------------
# 10. Окремі графіки для 10, 15 і 20 вузлів
# ---------------------------------

plt.figure(figsize=(8, 5))
plt.plot(distances[:10], elevations[:10], marker='o')
plt.title("Графік для 10 вузлів")
plt.xlabel("Кумулятивна відстань (м)")
plt.ylabel("Висота (м)")
plt.grid(True)

plt.figure(figsize=(8, 5))
plt.plot(distances[:15], elevations[:15], marker='o')
plt.title("Графік для 15 вузлів")
plt.xlabel("Кумулятивна відстань (м)")
plt.ylabel("Висота (м)")
plt.grid(True)

plt.figure(figsize=(8, 5))
plt.plot(distances[:20], elevations[:20], marker='o')
plt.title("Графік для 20 вузлів")
plt.xlabel("Кумулятивна відстань (м)")
plt.ylabel("Висота (м)")
plt.grid(True)

plt.show()

# ---------------------------------
# 12. Побудова сплайна і похибки
# ---------------------------------

def spline(x_val):
    for i in range(n_points - 1):
        if x[i] <= x_val <= x[i + 1]:
            dx = x_val - x[i]
            return (a_coeffs[i] +
                    b_coeffs[i] * dx +
                    c_coeffs[i] * dx**2 +
                    d_coeffs[i] * dx**3)
    return None

x_dense = np.linspace(x[0], x[-1], 500)

y_real = np.interp(x_dense, x, y)  
y_spline = [spline(val) for val in x_dense]

error = np.abs(y_real - y_spline)

plt.figure(figsize=(8,5))
plt.plot(x_dense, y_real, label="Реальні дані")
plt.plot(x_dense, y_spline, label="Сплайн", linestyle="--")

plt.scatter(x, y, color='red', label="Вузли")

plt.title("Сплайн-інтерполяція")
plt.xlabel("Відстань (м)")
plt.ylabel("Висота (м)")
plt.legend()
plt.grid(True)

plt.figure(figsize=(8,5))
plt.plot(x_dense, error)

plt.title("Похибка апроксимації")
plt.xlabel("Відстань (м)")
plt.ylabel("Похибка")
plt.grid(True)

plt.show()