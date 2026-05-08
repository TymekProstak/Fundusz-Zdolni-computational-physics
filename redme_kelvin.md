# MPI Vortex Sheet Simulation

Symulacja numeryczna ewolucji powierzchni wirowej z wykorzystaniem metody wirów dyskretnych oraz równoległości MPI.

Projekt dotyczy numerycznego badania rozwoju niestabilności powierzchni wirowej, w szczególności zjawiska podobnego do niestabilności Kelvina-Helmholtza. Powierzchnia wirowa jest reprezentowana jako zbiór punktów wirowych, które poruszają się zgodnie z prędkością indukowaną przez pozostałe punkty.

Do całkowania w czasie używany jest klasyczny schemat Rungego-Kutty czwartego rzędu. Obliczenia prędkości zostały zrównoleglone z użyciem MPI.

## Najważniejsze elementy projektu

- dyskretyzacja powierzchni wirowej,
- regularizowane jądro prędkości,
- całkowanie w czasie metodą RK4,
- równoległe obliczanie prędkości z użyciem MPI,
- automatyczna redyskretyzacja powierzchni wirowej,
- zapis wyników do plików tekstowych,
- zapis podstawowych wielkości diagnostycznych.

## Struktura repozytorium

```text
.
├── main.cpp
├── README.md
└── results/
    ├── wyniki_*.txt
    ├── predkosc_*.txt
    └── diagnostyka_*.txt
```

Folder `results/` jest tworzony automatycznie podczas działania programu.

## Wymagania

Do kompilacji wymagany jest kompilator C++ z obsługą standardu C++17 oraz biblioteka MPI.

Przykładowa instalacja na Ubuntu/Debian:

```bash
sudo apt update
sudo apt install build-essential openmpi-bin libopenmpi-dev
```

## Kompilacja

```bash
mpic++ -O3 -std=c++17 main.cpp -o vortex_sheet
```

## Uruchomienie

Program można uruchomić z domyślnymi parametrami:

```bash
mpirun -np 4 ./vortex_sheet
```

Domyślne parametry symulacji:

```text
L = 1.0
u1 = 2.0
u2 = 1.0
czas = 300.0
delta = 0.1
```

Można też podać parametry ręcznie:

```bash
mpirun -np 4 ./vortex_sheet L u1 u2 czas delta
```

Przykład:

```bash
mpirun -np 4 ./vortex_sheet 1.0 2.0 1.0 300.0 0.1
```

## Znaczenie parametrów

| Parametr | Znaczenie |
|---|---|
| `L` | długość początkowej powierzchni wirowej |
| `u1` | prędkość po jednej stronie powierzchni |
| `u2` | prędkość po drugiej stronie powierzchni |
| `czas` | całkowity czas symulacji |
| `delta` | parametr regularyzacji jądra prędkości |

Całkowita cyrkulacja jest liczona jako:

```text
G = u1 * L - u2 * L
```

Gęstość cyrkulacji jest następnie przypisywana punktom wirowym wzdłuż początkowej powierzchni.

## Model numeryczny

Powierzchnia wirowa jest reprezentowana przez punkty:

```cpp
struct wir {
    double x;
    double y;
    double epsilon;
    double gamma;
};
```

gdzie:

- `x`, `y` oznaczają położenie punktu wirowego,
- `epsilon` jest współrzędną parametryzującą powierzchnię,
- `gamma` jest gęstością cyrkulacji.

Dla każdego punktu obliczana jest prędkość indukowana przez całą powierzchnię wirową. Następnie położenia punktów są aktualizowane metodą Rungego-Kutty czwartego rzędu.

## Redyskretyzacja

W trakcie symulacji odległości między sąsiednimi punktami mogą rosnąć. Jeżeli odległość między dwoma sąsiednimi punktami przekroczy wartość krytyczną, dodawany jest nowy punkt pośredni.

Nowy punkt otrzymuje uśrednione wartości:

```text
x       = (x_i + x_{i-1}) / 2
y       = (y_i + y_{i-1}) / 2
epsilon = (epsilon_i + epsilon_{i-1}) / 2
gamma   = (gamma_i + gamma_{i-1}) / 2
```

Dzięki temu powierzchnia wirowa zachowuje dokładniejszą reprezentację w miejscach, gdzie zaczyna się mocniej deformować.

## Równoległość MPI

Każdy proces MPI oblicza prędkości tylko dla części punktów wirowych. Po zakończeniu lokalnych obliczeń dane są synchronizowane między wszystkimi procesami.

Do synchronizacji używane jest:

```cpp
MPI_Allgatherv(...)
```

Każdy proces po synchronizacji posiada pełny stan powierzchni wirowej i może przejść do następnego etapu metody RK4.

## Pliki wynikowe

Program zapisuje pliki do folderu `results/`.

### Pliki `wyniki_*.txt`

Przechowują położenie powierzchni wirowej w kolejnych krokach czasowych.

Format pliku:

```text
x_0 x_1 x_2 ... x_N y_0 y_1 y_2 ... y_N
```

Najpierw zapisywane są wszystkie współrzędne `x`, a następnie wszystkie współrzędne `y`.

Przykładowe nazwy:

```text
results/wyniki_0.1_0.txt
results/wyniki_0.1_10.txt
results/wyniki_0.1_20.txt
```

### Pliki `predkosc_*.txt`

Przechowują próbki pola prędkości na osi początkowej powierzchni wirowej.

Format pliku:

```text
vx_0 vx_1 vx_2 ... vx_N vy_0 vy_1 vy_2 ... vy_N
```

### Plik `diagnostyka_*.txt`

Przechowuje podstawowe wielkości diagnostyczne:

```text
krok liczba_wirow moment_x moment_y bezwladnosc
```

## Przykładowy przebieg pracy

Kompilacja:

```bash
mpic++ -O3 -std=c++17 main.cpp -o vortex_sheet
```

Uruchomienie na 4 procesach:

```bash
mpirun -np 4 ./vortex_sheet
```

Uruchomienie z własnymi parametrami:

```bash
mpirun -np 4 ./vortex_sheet 1.0 2.0 1.0 300.0 0.1
```

Po zakończeniu symulacji wyniki znajdują się w folderze:

```text
results/
```

## Przykładowa wizualizacja wyników

Pliki `wyniki_*.txt` można łatwo zwizualizować w Pythonie. Minimalny przykład odczytu jednego pliku:

```python
import matplotlib.pyplot as plt

filename = "results/wyniki_0.1_100.txt"

with open(filename, "r") as f:
    values = [float(v) for v in f.read().split()]

n = len(values) // 2

x = values[:n]
y = values[n:]

plt.plot(x, y)
plt.scatter(x, y, s=8)
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()
```

Do animacji można wykorzystać kolejne pliki `wyniki_*.txt` jako klatki czasowe symulacji.

## Uwagi

Program zapisuje wyniki tylko z procesu o randze 0.

Obliczenia prędkości są wykonywane równolegle przez wszystkie procesy MPI.

Parametr `delta` wpływa zarówno na regularizację jądra prędkości, jak i na początkową liczbę punktów dyskretyzacji.
