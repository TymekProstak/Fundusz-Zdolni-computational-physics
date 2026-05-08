# MPI Vortex Sheet Simulation

Symulacja numeryczna niestabilności Kelvina-Helmholtza z wykorzystaniem metody powierzchni wirowej, całkowania Rungego-Kutty czwartego rzędu oraz obliczeń równoległych MPI.

Projekt dotyczy numerycznego badania niestabilności rozwijającej się na granicy dwóch warstw płynu poruszających się z różnymi prędkościami. W takim układzie ścinanie na granicy warstw prowadzi do rozwoju zaburzeń, zwijania powierzchni wirowej oraz powstawania struktur wirowych.

W projekcie powierzchnia wirowa jest reprezentowana przez dyskretne punkty wirowe. Każdy punkt porusza się z prędkością indukowaną przez wszystkie pozostałe punkty. Po dyskretyzacji problem ma charakter problemu typu N-ciał, ponieważ każdy element oddziałuje z każdym innym elementem.

## Cel projektu

Celem projektu było odtworzenie i rozwinięcie wcześniejszych obliczeń numerycznych dostępnych w literaturze przedmiotu dotyczących niestabilności Kelvina-Helmholtza z wykorzystaniem współczesnych komputerów, obliczeń równoległych oraz większej precyzji numerycznej.

W szczególności badane były:

- przestrzenna ewolucja powierzchni wirowej,
- wpływ parametru regularyzacji `delta` na rozwój struktur wirowych,
- rozkład Fouriera indukowanej prędkości,
- przepływ przez powierzchnię `y=0`,
- porównanie wyników dla różnych zaburzeń początkowych,
- wpływ precyzji obliczeń na stabilność i jakość wyników.

Główną różnicą względem wcześniejszych rezultatów było wykorzystanie obliczeń w typie `double` zamiast `float` oraz zastosowanie obliczeń równoległych. Pozwoliło to na wydłużenie czasu symulacji, zmniejszenie skali dyskretyzacji oraz dokładniejsze badanie rozwoju struktur wirowych w porównaniu z oryginalnymi pracami numerycznymi z lat 80. i 90.

Dodatkowo, względem bazowego podejścia, badana była nie tylko czasowa ewolucja powierzchni wirowej, ale także jej widmo Fouriera oraz przepływ przez początkową powierzchnię rozdziału.

## Tło fizyczne

Niestabilność Kelvina-Helmholtza pojawia się na granicy dwóch warstw płynu poruszających się z różnymi prędkościami. Małe zaburzenie granicy między warstwami może zostać wzmocnione przez pole prędkości indukowane przez samą powierzchnię wirową. W wyniku tego powierzchnia zaczyna się zwijać i tworzyć charakterystyczne struktury wirowe.

W modelu przyjęto następujące założenia:

- brak lepkości,
- brak napięcia powierzchniowego,
- zachowanie cyrkulacji,
- periodyczne warunki brzegowe,
- początkowe zaburzenie sinusoidalne albo losowe,
- regularyzacja małych skal przez parametr `delta`.

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

Całkowita cyrkulacja jest liczona jako:

```text
G = u1 * L - u2 * L
```

gdzie:

- `L` jest długością początkowej powierzchni,
- `u1` jest prędkością po jednej stronie powierzchni,
- `u2` jest prędkością po drugiej stronie powierzchni.

Po zdyskretyzowaniu powierzchni, dla każdego punktu obliczana jest prędkość indukowana przez wszystkie pozostałe punkty. Następnie położenia punktów są aktualizowane metodą Rungego-Kutty czwartego rzędu.

## Regularyzacja

W idealnym, nielepkim modelu powierzchni wirowej niestabilność może rozwijać się w dowolnie małych skalach. Numerycznie prowadzi to do problemów, ponieważ coraz mniejsze struktury wymagałyby coraz gęstszej siatki i coraz większej dokładności.

Dlatego w modelu używany jest parametr regularyzacji `delta`.

Parametr `delta` ogranicza wpływ najmniejszych skal i można go interpretować jako nadanie powierzchni wirowej małej, ale skończonej grubości. Dla większego `delta` rozwój małych struktur jest silniej tłumiony. Dla mniejszego `delta` powierzchnia może tworzyć ostrzejsze i bardziej złożone struktury wirowe.

## Redyskretyzacja

W trakcie symulacji powierzchnia wirowa ulega rozciąganiu. Odległości między sąsiednimi punktami mogą więc rosnąć. Jeżeli odległość między dwoma sąsiednimi punktami przekroczy wartość krytyczną, dodawany jest nowy punkt pośredni.

Nowy punkt otrzymuje uśrednione wartości:

```text
x       = (x_i + x_{i-1}) / 2
y       = (y_i + y_{i-1}) / 2
epsilon = (epsilon_i + epsilon_{i-1}) / 2
gamma   = (gamma_i + gamma_{i-1}) / 2
```

Dzięki temu powierzchnia wirowa zachowuje dokładniejszą reprezentację w miejscach, gdzie zaczyna się mocniej deformować.

## Obliczenia równoległe

Obliczanie prędkości dla różnych punktów jest od siebie niezależne. Z tego powodu główna część obliczeń może być wykonywana równolegle.

W projekcie wykorzystano bibliotekę MPI. Każdy proces oblicza prędkości dla własnego fragmentu tablicy punktów wirowych. Po zakończeniu lokalnych obliczeń procesy wymieniają dane, tak aby każdy proces znał pełny stan powierzchni przed kolejnym etapem metody RK4.

Do synchronizacji używane jest:

```cpp
MPI_Allgatherv(...)
```

Taki schemat pozwala przyspieszyć najbardziej kosztowny etap programu, czyli obliczanie oddziaływań między punktami.

## Struktura repozytorium

```text
.
├── main.cpp
├── README.md
├── scripts/
│   └── analyze_vortex_results.py
├── results/
│   ├── wyniki_*.txt
│   ├── predkosc_*.txt
│   └── diagnostyka_*.txt
└── figures/
    ├── plots/
    ├── animations/
    └── comparisons/
```

Folder `results/` jest tworzony przez program C++.

Folder `figures/` jest tworzony przez skrypt Pythona.

## Wymagania

### C++ i MPI

Do kompilacji wymagany jest kompilator C++ z obsługą standardu C++17 oraz biblioteka MPI.

Przykładowa instalacja na Ubuntu/Debian:

```bash
sudo apt update
sudo apt install build-essential openmpi-bin libopenmpi-dev
```

### Python

Do analizy i generowania wykresów wymagane są:

```bash
pip install numpy matplotlib pillow
```

## Kompilacja

```bash
mpic++ -O3 -std=c++17 main.cpp -o vortex_sheet
```

## Uruchomienie symulacji

Program można uruchomić z domyślnymi parametrami:

```bash
mpirun -np 4 ./vortex_sheet
```

Domyślne parametry:

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

## Pliki wynikowe

Program zapisuje wyniki do folderu `results/`.

### Pliki `wyniki_*.txt`

Pliki przechowują położenie powierzchni wirowej w kolejnych krokach czasowych.

Format:

```text
x_0 x_1 x_2 ... x_N y_0 y_1 y_2 ... y_N
```

Najpierw zapisywane są wszystkie współrzędne `x`, a następnie wszystkie współrzędne `y`.

Przykład:

```text
results/wyniki_0.1_0.txt
results/wyniki_0.1_10.txt
results/wyniki_0.1_20.txt
```

### Pliki `predkosc_*.txt`

Pliki przechowują próbki pola prędkości na osi `y=0`.

Format:

```text
vx_0 vx_1 vx_2 ... vx_N vy_0 vy_1 vy_2 ... vy_N
```

Na ich podstawie wykonywana jest analiza Fouriera oraz obliczany jest przepływ przez `y=0`.

### Pliki `diagnostyka_*.txt`

Pliki diagnostyczne zawierają podstawowe wielkości globalne:

```text
krok liczba_wirow moment_x moment_y bezwladnosc
```

## Analiza wyników

Do analizy wyników służy skrypt:

```text
scripts/analyze_vortex_results.py
```

Uruchomienie podstawowe:

```bash
python3 scripts/analyze_vortex_results.py \
    --input results \
    --output figures \
    --L 1.0 \
    --u1 2.0 \
    --u2 1.0 \
    --times 0 1 2 4.5 \
    --compare-time 4.5 \
    --fps 12
```

Skrypt generuje:

- wykresy powierzchni wirowej w kilku chwilach czasu,
- widma Fouriera w kilku chwilach czasu,
- wykres przepływu przez `y=0`,
- animację GIF powierzchni wirowej i widma,
- porównania między różnymi wartościami `delta`.

## Przykładowe wyniki

Po uruchomieniu skryptu w folderze `figures/` pojawią się pliki podobne do:

```text
figures/plots/sheet_snapshots_delta_0.1.png
figures/plots/spectrum_snapshots_delta_0.1.png
figures/plots/flow_delta_0.1.png
figures/animations/vortex_delta_0.1.gif
figures/comparisons/comparison_sheet_t_4.5.png
figures/comparisons/comparison_spectrum_t_4.5.png
figures/comparisons/comparison_flow_delta.png
```

Przykład umieszczenia wyników w README:

```markdown
![Ewolucja powierzchni wirowej](figures/plots/sheet_snapshots_delta_0.1.png)

![Widmo Fouriera](figures/plots/spectrum_snapshots_delta_0.1.png)

![Przepływ przez y=0](figures/plots/flow_delta_0.1.png)

![Animacja powierzchni wirowej](figures/animations/vortex_delta_0.1.gif)
```

## Porównanie wpływu regularyzacji

Dla wielu wartości `delta` można uruchomić analizę:

```bash
python3 scripts/analyze_vortex_results.py \
    --input results \
    --output figures \
    --deltas 0.5 0.25 0.05 0.01 \
    --times 0 1 2 4.5 \
    --compare-time 4.5
```

Skrypt wygeneruje porównania:

```text
figures/comparisons/comparison_sheet_t_4.5.png
figures/comparisons/comparison_spectrum_t_4.5.png
figures/comparisons/comparison_flow_delta.png
```

Przykładowe wstawienie do README:

```markdown
![Porównanie powierzchni wirowej dla różnych delta](figures/comparisons/comparison_sheet_t_4.5.png)

![Porównanie widma Fouriera dla różnych delta](figures/comparisons/comparison_spectrum_t_4.5.png)

![Porównanie przepływu dla różnych delta](figures/comparisons/comparison_flow_delta.png)
```

## Interpretacja wykresów

### Powierzchnia wirowa

Wykres powierzchni wirowej pokazuje przestrzenną ewolucję punktów wirowych. Dla mniejszych wartości `delta` mogą pojawiać się mniejsze i ostrzejsze struktury. Dla większych wartości `delta` rozwój małych skal jest silniej tłumiony.

W praktyce parametr `delta` kontroluje to, jak drobne struktury mogą pojawić się w rozwiązaniu. Większa wartość regularyzacji wygładza powierzchnię wirową i opóźnia powstawanie mniejszych struktur. Mniejsza wartość regularyzacji pozwala na szybszy rozwój małych skal, przez co rozwiązanie ma bardziej złożony, lokalnie turbulentny charakter.

### Widmo Fouriera

Widmo Fouriera pokazuje, które skale przestrzenne są obecne w rozwiązaniu. Niskie mody odpowiadają dużym strukturom przestrzennym, a wysokie mody odpowiadają małym skalom.

Dla większych wartości `delta` maksimum widma zwykle przypada na charakterystyczną skalę głównych wirów widocznych w danej chwili czasu. Oznacza to, że w rozwiązaniu dominuje jedna większa struktura lub kilka struktur o podobnej skali.

Dla małych wartości `delta` rozwiązanie ma bardziej złożony charakter. W początkowej fazie energia może pojawiać się w wyższych modach, co odpowiada rozwojowi mniejszych struktur. Następnie, wraz z postępem zwijania powierzchni, mniejsze wiry zaczynają łączyć się i zawijać wokół większych ośrodków wirowych. W takim przypadku maksimum widma może przesuwać się od modów odpowiadających małym zaburzeniom w stronę modów związanych z większymi strukturami, które powstają przez grupowanie i zwijanie mniejszych wirów.

W przypadku bardziej turbulentnego charakteru rozwiązania widmo nie opisuje już tylko jednej dominującej skali. Zamiast tego pokazuje obecność wielu skal jednocześnie: dużych wirów organizujących ruch oraz mniejszych struktur powstających wewnątrz nich.

### Przepływ przez `y=0`

Przepływ przez `y=0` jest liczony na podstawie składowej pionowej prędkości `Uy` próbkowanej na osi `y=0`.

W skrypcie liczona jest głównie wielkość:

```text
integral |Uy| dx
```

czyli miara intensywności przepływu przez początkową powierzchnię rozdziału. Wielkość ze znakiem:

```text
integral Uy dx
```

jest również zapisywana na wykresie, ale w układzie periodycznym może się częściowo znosić.

## Porównanie z wynikami referencyjnymi

Na końcu analizy można dodać porównanie wyników uzyskanych w tej implementacji z wcześniejszymi wynikami referencyjnymi lub z wcześniejszą wersją wizualizacji. Najwygodniej zrobić to przez wstawienie jednego obrazu porównawczego, na przykład zestawienia starej grafiki i nowego wykresu wygenerowanego przez skrypt.

Przykładowa struktura plików:

```text
figures/comparisons/reference_vs_current.png
```

Przykładowe wstawienie do README:

```markdown
![Porównanie wyniku referencyjnego z aktualną symulacją](figures/comparisons/reference_vs_current.png)
```

Taki obraz można przygotować ręcznie, łącząc zrzut wcześniejszej grafiki z aktualnym wykresem wygenerowanym przez skrypt. Dzięki temu porównanie jest bardziej czytelne niż automatyczne zestawianie surowych plików wynikowych, szczególnie jeżeli stare dane były zapisane w innym formacie albo pochodziły bezpośrednio z artykułu lub wcześniejszej prezentacji.

## Najważniejsze obserwacje

Typowe wyniki są zgodne z oczekiwaniami dla niestabilności Kelvina-Helmholtza:

- powierzchnia wirowa zaczyna się zwijać i tworzyć struktury wirowe,
- parametr `delta` kontroluje tempo i skalę rozwoju wirów,
- mniejsze `delta` pozwala na rozwój drobniejszych struktur,
- większe `delta` działa wygładzająco,
- widmo Fouriera odzwierciedla obecność małych i dużych skal,
- dla dużych wartości `delta` widmo jest zwykle związane z charakterystyczną skalą głównych wirów,
- dla małych wartości `delta` energia widma może obejmować wiele modów, ponieważ rozwiązanie zawiera jednocześnie większe struktury i drobniejsze zaburzenia,
- z czasem mniejsze struktury mogą zawijać się wokół większych ośrodków wirowych, co wpływa na przesuwanie maksimum widma,
- przepływ przez `y=0` opisuje intensywność mieszania między warstwami,
- zaburzenie początkowe wpływa na charakter późniejszych struktur.

## Bibliografia

[1] Robert Krasny,  
**Desingularization of Periodic Vortex Sheet Roll-up**,  
Journal of Computational Physics, Volume 65, Issue 2, 1986, pp. 292-313.  
DOI: `10.1016/0021-9991(86)90210-X`

[2] Vladimir Parezanović, Jean-Charles Laurentie, Carine Fourment, Joel Delville, Jean-Paul Bonnet, Andreas Spohn, Thomas Duriez, Laurent Cordier, Bernd R. Noack, Markus Abel, Marc Segond, Tamir Shaqarin, Steven L. Brunton,  
**Mixing Layer Manipulation Experiment: From Open-Loop Forcing to Closed-Loop Machine Learning Control**,  
Flow, Turbulence and Combustion, Volume 94, 2015, pp. 155-173.  
DOI: `10.1007/s10494-014-9581-1`

[3] Kelvin-Helmholtz instability,  
Wikipedia, dostęp: 31.05.2024.  
`https://en.wikipedia.org/wiki/Kelvin%E2%80%93Helmholtz_instability`
