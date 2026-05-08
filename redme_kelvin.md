# MPI Vortex Sheet Simulation

Numerical simulation of the Kelvin-Helmholtz instability using the vortex sheet method, fourth-order Runge-Kutta time integration, and parallel MPI computations.

The project focuses on the numerical investigation of an instability developing at the interface between two fluid layers moving with different velocities. In such a system, shear at the interface leads to the growth of perturbations, vortex sheet roll-up, and the formation of vortex structures.

In this project, the vortex sheet is represented by discrete vortex points. Each point moves with the velocity induced by all the remaining points. After discretization, the problem has the character of an N-body problem, because every element interacts with every other element.

## Project Goal

The goal of the project was to reproduce and extend earlier numerical results available in the literature on the Kelvin-Helmholtz instability, using modern computers, parallel computation, and increased numerical precision.

In particular, the project investigates:

- spatial evolution of the vortex sheet,
- influence of the regularization parameter `delta` on the development of vortex structures,
- Fourier spectrum of the induced velocity,
- flow through the surface `y=0`,
- comparison of results for different initial perturbations,
- influence of numerical precision on stability and quality of the results.

The main difference compared with earlier results was the use of `double` precision instead of `float`, together with parallel computation. This made it possible to extend the simulation time, reduce the discretization scale, and study the development of vortex structures more accurately in comparison with original numerical works from the 1980s and 1990s.

Additionally, compared with the basic approach, the project studies not only the time evolution of the vortex sheet, but also its Fourier spectrum and the flow through the initial separation surface.

## Physical Background

The Kelvin-Helmholtz instability appears at the interface between two fluid layers moving with different velocities. A small perturbation of the interface can be amplified by the velocity field induced by the vortex sheet itself. As a result, the sheet begins to roll up and form characteristic vortex structures.

The following assumptions are used in the model:

- no viscosity,
- no surface tension,
- conservation of circulation,
- periodic boundary conditions,
- sinusoidal or random initial perturbation,
- regularization of small scales by the parameter `delta`.

## Numerical Model

The vortex sheet is represented by points:

```cpp
struct wir {
    double x;
    double y;
    double epsilon;
    double gamma;
};
```

where:

- `x`, `y` denote the position of a vortex point,
- `epsilon` is the coordinate parametrizing the sheet,
- `gamma` is the circulation density.

The total circulation is computed as:

```text
G = u1 * L - u2 * L
```

where:

- `L` is the initial length of the sheet,
- `u1` is the velocity on one side of the sheet,
- `u2` is the velocity on the other side of the sheet.

After discretizing the sheet, the induced velocity is computed for every point from all remaining points. Then the point positions are updated using the fourth-order Runge-Kutta method.

## Regularization

In an ideal inviscid vortex sheet model, the instability may develop at arbitrarily small scales. Numerically, this creates difficulties because smaller and smaller structures would require increasingly dense discretization and higher accuracy.

Therefore, the model uses the regularization parameter `delta`.

The parameter `delta` limits the influence of the smallest scales and can be interpreted as assigning a small but finite thickness to the vortex sheet. For larger `delta`, the development of small structures is more strongly damped. For smaller `delta`, the sheet can form sharper and more complex vortex structures.

## Rediscretization

During the simulation, the vortex sheet stretches. As a result, the distances between neighboring points may increase. If the distance between two neighboring points exceeds a critical value, a new intermediate point is added.

The new point receives averaged values:

```text
x       = (x_i + x_{i-1}) / 2
y       = (y_i + y_{i-1}) / 2
epsilon = (epsilon_i + epsilon_{i-1}) / 2
gamma   = (gamma_i + gamma_{i-1}) / 2
```

This allows the vortex sheet to maintain a more accurate representation in regions where it becomes strongly deformed.

## Parallel Computation

Computing the velocity for different points is independent. For this reason, the main computational part can be parallelized.

The project uses the MPI library. Each process computes velocities for its own part of the vortex point array. After local computations are finished, the processes exchange data so that each process knows the full state of the sheet before the next RK4 stage.

Synchronization is performed using:

```cpp
MPI_Allgatherv(...)
```

This approach accelerates the most expensive part of the program, namely computing interactions between points.

## Repository Structure

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

The `results/` folder is created by the C++ program.

The `figures/` folder is created by the Python script.

## Requirements

### C++ and MPI

Compilation requires a C++ compiler with C++17 support and an MPI library.

Example installation on Ubuntu/Debian:

```bash
sudo apt update
sudo apt install build-essential openmpi-bin libopenmpi-dev
```

### Python

For analysis and plot generation, the following packages are required:

```bash
pip install numpy matplotlib pillow
```

## Compilation

```bash
mpic++ -O3 -std=c++17 main.cpp -o vortex_sheet
```

## Running the Simulation

The program can be run with default parameters:

```bash
mpirun -np 4 ./vortex_sheet
```

Default parameters:

```text
L = 1.0
u1 = 2.0
u2 = 1.0
time = 300.0
delta = 0.1
```

Parameters can also be passed manually:

```bash
mpirun -np 4 ./vortex_sheet L u1 u2 time delta
```

Example:

```bash
mpirun -np 4 ./vortex_sheet 1.0 2.0 1.0 300.0 0.1
```

## Output Files

The program saves results to the `results/` folder.

### `wyniki_*.txt` Files

These files store the position of the vortex sheet at consecutive time steps.

Format:

```text
x_0 x_1 x_2 ... x_N y_0 y_1 y_2 ... y_N
```

First all `x` coordinates are written, followed by all `y` coordinates.

Example:

```text
results/wyniki_0.1_0.txt
results/wyniki_0.1_10.txt
results/wyniki_0.1_20.txt
```

### `predkosc_*.txt` Files

These files store samples of the velocity field along the axis `y=0`.

Format:

```text
vx_0 vx_1 vx_2 ... vx_N vy_0 vy_1 vy_2 ... vy_N
```

They are used to perform Fourier analysis and to compute the flow through `y=0`.

### `diagnostyka_*.txt` Files

Diagnostic files contain basic global quantities:

```text
step number_of_vortices moment_x moment_y inertia
```

## Result Analysis

The results are analyzed using the script:

```text
scripts/analyze_vortex_results.py
```

Basic usage:

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

The script generates:

- plots of the vortex sheet at several time instants,
- Fourier spectra at several time instants,
- a plot of the flow through `y=0`,
- a GIF animation of the vortex sheet and its spectrum,
- comparisons between different values of `delta`.

## Example Results

After running the script, files similar to the following will appear in the `figures/` folder:

```text
figures/plots/sheet_snapshots_delta_0.1.png
figures/plots/spectrum_snapshots_delta_0.1.png
figures/plots/flow_delta_0.1.png
figures/animations/vortex_delta_0.1.gif
figures/comparisons/comparison_sheet_t_4.5.png
figures/comparisons/comparison_spectrum_t_4.5.png
figures/comparisons/comparison_flow_delta.png
```

Example of embedding results in the README:

```markdown
![Vortex sheet evolution](figures/plots/sheet_snapshots_delta_0.1.png)

![Fourier spectrum](figures/plots/spectrum_snapshots_delta_0.1.png)

![Flow through y=0](figures/plots/flow_delta_0.1.png)

![Vortex sheet animation](figures/animations/vortex_delta_0.1.gif)
```

## Comparison of the Regularization Effect

For multiple values of `delta`, the analysis can be run as:

```bash
python3 scripts/analyze_vortex_results.py \
    --input results \
    --output figures \
    --deltas 0.5 0.25 0.05 0.01 \
    --times 0 1 2 4.5 \
    --compare-time 4.5
```

The script will generate comparisons:

```text
figures/comparisons/comparison_sheet_t_4.5.png
figures/comparisons/comparison_spectrum_t_4.5.png
figures/comparisons/comparison_flow_delta.png
```

Example of embedding comparison plots in the README:

```markdown
![Comparison of the vortex sheet for different delta values](figures/comparisons/comparison_sheet_t_4.5.png)

![Comparison of Fourier spectra for different delta values](figures/comparisons/comparison_spectrum_t_4.5.png)

![Comparison of flow for different delta values](figures/comparisons/comparison_flow_delta.png)
```

## Interpretation of Plots

### Vortex Sheet

The vortex sheet plot shows the spatial evolution of vortex points. For smaller values of `delta`, smaller and sharper structures may appear. For larger values of `delta`, the development of small scales is more strongly damped.

In practice, the parameter `delta` controls how fine the structures in the solution can become. A larger regularization value smooths the vortex sheet and delays the formation of smaller structures. A smaller regularization value allows small scales to develop faster, which makes the solution more complex and locally turbulent.

### Fourier Spectrum

The Fourier spectrum shows which spatial scales are present in the solution. Low modes correspond to large spatial structures, while high modes correspond to small scales.

For larger values of `delta`, the maximum of the spectrum usually corresponds to the characteristic scale of the main vortices visible at a given time. This means that the solution is dominated by one larger structure or by several structures of similar scale.

For small values of `delta`, the solution has a more complex character. In the initial phase, energy may appear in higher modes, which corresponds to the development of smaller structures. Later, as the vortex sheet rolls up, smaller vortices begin to merge and wrap around larger vortex centers. In such a case, the maximum of the spectrum may shift from modes corresponding to small perturbations toward modes associated with larger structures formed by the grouping and rolling-up of smaller vortices.

For solutions with a more turbulent character, the spectrum no longer describes only one dominant scale. Instead, it shows the simultaneous presence of many scales: large vortices organizing the motion and smaller structures forming inside them.

### Flow Through `y=0`

The flow through `y=0` is computed from the vertical velocity component `Uy` sampled along the axis `y=0`.

The script mainly computes:

```text
integral |Uy| dx
```

which is a measure of the intensity of flow through the initial separation surface. The signed quantity:

```text
integral Uy dx
```

is also plotted, but in a periodic system it may partially cancel out.

## Comparison with Reference Results

At the end of the analysis, a comparison between the results obtained in this implementation and earlier reference results or an older visualization can be added. The most convenient way is to insert a single comparison image, for example a side-by-side figure showing the old graphic and the new plot generated by the script.

Example file structure:

```text
figures/comparisons/reference_vs_current.png
```

Example of embedding it in the README:

```markdown
![Comparison of the reference result with the current simulation](figures/comparisons/reference_vs_current.png)
```

Such an image can be prepared manually by combining a screenshot of the earlier graphic with the current plot generated by the script. This makes the comparison clearer than automatic comparison of raw output files, especially if the old data were stored in a different format or came directly from an article or an earlier presentation.

## Main Observations

Typical results are consistent with the expected behavior of the Kelvin-Helmholtz instability:

- the vortex sheet begins to roll up and form vortex structures,
- the parameter `delta` controls the rate and scale of vortex development,
- smaller `delta` allows finer structures to develop,
- larger `delta` has a smoothing effect,
- the Fourier spectrum reflects the presence of both small and large scales,
- for large values of `delta`, the spectrum is usually related to the characteristic scale of the main vortices,
- for small values of `delta`, the spectral energy may cover many modes because the solution contains both larger structures and smaller perturbations,
- over time, smaller structures may wrap around larger vortex centers, which influences the shift of the spectral maximum,
- the flow through `y=0` describes the intensity of mixing between the layers,
- the initial perturbation affects the character of later structures.

## Bibliography

[1] Robert Krasny,  
**Desingularization of Periodic Vortex Sheet Roll-up**,  
Journal of Computational Physics, Volume 65, Issue 2, 1986, pp. 292-313.  
DOI: `10.1016/0021-9991(86)90210-X`

[2] Vladimir Parezanović, Jean-Charles Laurentie, Carine Fourment, Joel Delville, Jean-Paul Bonnet, Andreas Spohn, Thomas Duriez, Laurent Cordier, Bernd R. Noack, Markus Abel, Marc Segond, Tamir Shaqarin, Steven L. Brunton,  
**Mixing Layer Manipulation Experiment: From Open-Loop Forcing to Closed-Loop Machine Learning Control**,  
Flow, Turbulence and Combustion, Volume 94, 2015, pp. 155-173.  
DOI: `10.1007/s10494-014-9581-1`

[3] Kelvin-Helmholtz instability,  
Wikipedia, accessed: 31.05.2024.  
`https://en.wikipedia.org/wiki/Kelvin%E2%80%93Helmholtz_instability`
