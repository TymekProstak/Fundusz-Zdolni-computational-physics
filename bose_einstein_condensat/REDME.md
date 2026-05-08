# Bose-Einstein Condensate Gross-Pitaevskii Simulation

Numerical simulation of a two-dimensional Bose-Einstein condensate using the Gross-Pitaevskii equation.

The project focuses on the numerical investigation of vortex structures in a rotating Bose-Einstein condensate placed in external trapping potentials. The simulation includes imaginary-time evolution for obtaining approximate stationary states and real-time evolution for studying the later dynamics of the condensate after modifying the potential.

The condensate is described by a complex wave function. Its squared modulus represents the particle density, while the phase of the wave function is related to the local flow and vortex structure.

---

## Project Goal and Main Contribution

The goal of the project was to investigate vortex formation and phase-imprinting effects in a two-dimensional rotating Bose-Einstein condensate described by the Gross-Pitaevskii equation.

**The main contribution of this work is a numerical study of potential-induced phase imprinting in a rotating Bose-Einstein condensate.** In particular, the project investigates whether a vortex and density structure prepared by an external potential can remain visible in the real-time dynamics of the condensate after this potential is switched off.

**The key physical idea is that the external potential can imprint structural information into the condensate wave function.** This means that the potential does not only shape the instantaneous density distribution, but can also influence the later evolution of the condensate after the potential has been removed.

Compared with a basic Gross-Pitaevskii simulation, this project focuses not only on obtaining stationary vortex states, but also on analyzing how the parameters of the external potential and the condensate determine the later vortex dynamics. The simulations study the relation between the optical lattice depth, lattice wave numbers, nonlinear interaction strength, angular velocity, and the resulting local and global motion of the condensate.

In particular, the project studies:

- density structure of the condensate,
- vortex formation in a rotating frame,
- influence of angular velocity `omega`,
- influence of nonlinear interaction strength `beta`,
- influence of the optical-lattice depth `V0`,
- influence of the optical-lattice wave numbers `kappa1` and `kappa2`,
- imaginary-time evolution used to obtain stationary or prepared states,
- real-time evolution after changing or switching off part of the potential,
- angular momentum of the condensate,
- relation between local vortex structures and global rotation of the condensate,
- whether the later vortex structure preserves information about the initially prepared state.

**The results suggest that potential-induced phase imprinting can persist in a rotating Bose-Einstein condensate even after the imprinting potential has been switched off.** This imprint remains visible in the later real-time evolution through the local vortex arrangement and the global rotation of the condensate density pattern.

Therefore, the physical contribution of the project can be summarized as follows:

**An external potential can prepare a characteristic vortex and density structure in a rotating condensate, and this structure can continue to influence the condensate dynamics after the potential is removed.**
---

## Physical Background

A Bose-Einstein condensate is an ultracold quantum gas in which a macroscopic number of bosons occupy the same quantum state. In the mean-field approximation, the condensate can be described by a single complex wave function:

```math
\psi(x,y,t)
```

The particle density is given by:

```math
\rho(x,y,t)=|\psi(x,y,t)|^2
```

The wave function is normalized so that:

```math
\int |\psi(x,y,t)|^2\,dx\,dy = 1
```

The evolution of the condensate is modeled using the Gross-Pitaevskii equation.

---

## Gross-Pitaevskii Equation

In dimensionless form, the rotating two-dimensional Gross-Pitaevskii equation can be written as:

```math
i\frac{\partial \psi}{\partial t}
=
\left[
-\frac{1}{2}\nabla^2
+
V(x,y)
+
\beta |\psi|^2
-
\Omega L_z
\right]\psi
```

where:

- `psi` is the condensate wave function,
- `V(x,y)` is the external trapping potential,
- `beta` controls the strength of nonlinear interactions,
- `Omega` is the angular velocity of the rotating frame,
- `L_z` is the angular momentum operator,
- `nabla^2` is the two-dimensional Laplace operator.

The angular momentum operator is:

```math
L_z
=
-i
\left(
x\frac{\partial}{\partial y}
-
y\frac{\partial}{\partial x}
\right)
```

The nonlinear term:

```math
\beta |\psi|^2
```

models mean-field interactions between particles in the condensate.

---

## Potentials

Two main types of potentials are used in the simulation.

### Harmonic Potential

The harmonic trapping potential is:

```math
V_{\mathrm{harm}}(x,y)
=
x^2+y^2
```

It confines the condensate near the center of the computational domain.

### Optical Lattice Potential

The optical lattice is modeled as:

```math
V_{\mathrm{lat}}(x,y)
=
V_0
\left[
\sin^2(\kappa_1 x)
+
\sin^2(\kappa_2 y)
\right]
```

where:

- `V0` is the lattice depth,
- `kappa1` is the wave number in the `x` direction,
- `kappa2` is the wave number in the `y` direction.

The total potential can include both harmonic confinement and the optical lattice:

```math
V(x,y)
=
V_{\mathrm{harm}}(x,y)
+
V_{\mathrm{lat}}(x,y)
```

depending on the selected simulation setup.

---

## Initial State

The initial wave function can be chosen in the form:

```math
\psi_0(x,y)
=
\exp
\left[
-p(x^2+y^2)
+
ik\sqrt{x^2+y^2}
\right]
```

where:

- `p` controls the initial spatial width of the condensate,
- `k` controls the initial radial phase profile.

After initialization, the wave function is normalized numerically:

```math
\psi
\leftarrow
\frac{\psi}{\sqrt{\int |\psi|^2\,dx\,dy}}
```

The phase term defines the initial phase profile of the wave function. In the main phase-imprinting analysis, however, the important effect comes from preparing the condensate in a structured potential and then observing its later real-time evolution after the potential is modified.

---

## Imaginary-Time Evolution

Imaginary-time evolution is used to obtain an approximate stationary state of the condensate.

The idea is to substitute:

```math
t = -i\tau
```

which transforms the Schrödinger-type evolution into a relaxation-like equation:

```math
\frac{\partial \psi}{\partial \tau}
=
-
H[\psi]\psi
```

High-energy components decay faster, so after many iterations the wave function approaches a low-energy stationary state.

Because the norm would otherwise decay during imaginary-time propagation, the wave function is normalized after every iteration:

```math
\int |\psi|^2\,dx\,dy = 1
```

In this project, imaginary-time evolution is used to prepare the condensate in a chosen external potential before analyzing vortex structures or starting real-time evolution.

---

## Real-Time Evolution

Real-time evolution is used to study the dynamics of the condensate after a prepared state has been obtained.

In this mode, the equation keeps its Schrödinger-like form:

```math
i\frac{\partial \psi}{\partial t}
=
H[\psi]\psi
```

This allows observing how the density and vortex structure evolve after the potential is changed, for example after switching off the optical lattice.

A typical workflow is:

1. Compute a prepared state using imaginary-time evolution.
2. Save the final wave function.
3. Use this wave function as the initial condition for real-time evolution.
4. Modify the potential, for example by removing the optical lattice.
5. Observe the later density and vortex dynamics.

---

## Phase Imprinting and Real-Time Evolution

In this project, phase imprinting is understood as preparing a characteristic condensate structure using an external potential and then observing how this structure affects later real-time evolution.

The condensate is first evolved in a selected potential, for example in a harmonic trap or in an optical lattice. This potential shapes the density and vortex structure of the wave function. After this prepared state is obtained, the system can be evolved in real time after modifying or switching off part of the potential.

A typical phase-imprinting workflow is:

1. Prepare the condensate state in a selected external potential.
2. Obtain a characteristic density and vortex structure.
3. Start real-time evolution from this prepared state.
4. Modify the potential, for example by switching off the optical lattice.
5. Observe whether the later vortex pattern preserves information about the initially prepared structure.

In this sense, the initial potential imprints information into the condensate wave function. The later dynamics is not independent of the preparation stage: similar initial vortex structures may lead to similar later geometries.

The rotation and deformation of the vortex pattern can also be compared with the angular momentum of the condensate:

```math
\langle L_z \rangle
=
\int
\psi^*
L_z
\psi
\,dx\,dy
```

The simulations therefore connect three elements:

- the potential used to prepare the condensate,
- the vortex and density structure obtained after preparation,
- the later real-time evolution after the potential is changed.

---

## Numerical Discretization

The simulation is performed on a two-dimensional Cartesian grid:

```math
x \in [x_{\min},x_{\max}]
```

```math
y \in [y_{\min},y_{\max}]
```

with:

```math
N_x \times N_y
```

grid points.

The spatial steps are:

```math
\Delta x =
\frac{x_{\max}-x_{\min}}{N_x-1}
```

```math
\Delta y =
\frac{y_{\max}-y_{\min}}{N_y-1}
```

The Laplace operator is approximated using finite differences. The kinetic energy operator is:

```math
T
=
-\frac{1}{2}\nabla^2
```

The Hamiltonian used in the simulation has the form:

```math
H[\psi]
=
T
+
V(x,y)
+
\beta|\psi|^2
-
\Omega L_z
```

The nonlinear part depends on the current wave function, so the Hamiltonian must be updated during the iteration.

---

## Angular Momentum

The angular momentum is computed as the expectation value of the angular momentum operator:

```math
\langle L_z \rangle
=
\int
\psi^*
L_z
\psi
\,dx\,dy
```

This quantity is used to characterize the rotating state and to compare configurations with different vortex structures.

In the numerical implementation, the integral is approximated on the discrete grid:

```math
\langle L_z \rangle
\approx
\Delta x \Delta y
\sum_{i,j}
\psi_{i,j}^{*}
\left(L_z\psi\right)_{i,j}
```

---

## Repository Structure

```text
.
├── bec_simulation.py
├── redme_bec.md
├── results/
│   └── run_name/
│       ├── config.json
│       ├── psi_ground.npy
│       ├── psi_final.npy
│       ├── psi_history.npy
│       ├── angular_momentum.npy
│       ├── potential.npy
│       ├── density_ground.png
│       ├── density_final.png
│       ├── potential.png
│       └── density_animation.gif
└── figures/
    ├── density_ground.png
    ├── density_final.png
    ├── potential.png
    └── density_animation.gif
```

The `results/` directory contains numerical output from individual runs.

Each run should have its own folder. This prevents overwriting older results and makes it easier to compare different parameter sets.

---

## Requirements

The simulation requires Python and the following packages:

```bash
pip install numpy scipy matplotlib pillow
```

The main numerical dependencies are:

- `numpy` for array operations,
- `scipy` for sparse matrices and linear solvers,
- `matplotlib` for visualization,
- `pillow` for GIF animation export.

---

## Running the Simulation

The simulation should be run from the command line.

Example imaginary-time run:

```bash
python3 bec_simulation.py \
    --mode imaginary \
    --nx 401 \
    --ny 401 \
    --domain 12 \
    --beta 500 \
    --omega 3.0 \
    --potential lattice \
    --V0 15 \
    --kappa1 1.0471975512 \
    --kappa2 1.0471975512 \
    --steps 1000 \
    --output results/beta500_omega3_lattice
```

Example real-time run:

```bash
python3 bec_simulation.py \
    --mode real \
    --load results/beta500_omega3_lattice/psi_ground.npy \
    --steps 5000 \
    --dt 0.001 \
    --turn-off-lattice \
    --output results/real_time_lattice_off
```

---

## Parameters

Important parameters:

```text
--mode              imaginary or real
--nx                number of grid points in x direction
--ny                number of grid points in y direction
--domain            half-width of the square computational domain
--beta              nonlinear interaction strength
--omega             angular velocity
--potential         harmonic, lattice, combined, or free
--V0                optical lattice depth
--kappa1            optical lattice wave number in x direction
--kappa2            optical lattice wave number in y direction
--steps             number of simulation iterations
--dt                time step
--output            output directory
--load              input wave function for real-time evolution
--save-every        interval for saving frames
```

The parameters of every run should be saved automatically to:

```text
config.json
```

This makes each result reproducible.

---

## Output Files

Each simulation run saves output to its own folder.

Typical files are:

```text
config.json
psi_ground.npy
psi_final.npy
psi_history.npy
angular_momentum.npy
potential.npy
density_ground.png
density_final.png
potential.png
density_animation.gif
```

### `config.json`

Stores all parameters used in the run.

### `psi_ground.npy`

Stores the final wave function after imaginary-time evolution.

### `psi_final.npy`

Stores the final wave function after real-time evolution.

### `psi_history.npy`

Stores selected wave functions saved during evolution.

### `angular_momentum.npy`

Stores the angular momentum values computed during the simulation.

### `potential.npy`

Stores the external potential used in the simulation.

### `density_*.png`

Images of the condensate density:

```math
|\psi(x,y)|^2
```

### `density_animation.gif`

Animation of the density evolution.

---

## Generated Results

Example visualizations produced by the project:

### External Potential

![External potential](figures/potential.png)

### Ground-State Density

![Ground-state density](figures/density_ground.png)

### Final Density

![Final density](figures/density_final.png)

### Density Evolution

![Density evolution](figures/density_animation.gif)

---

## Typical Workflow

A typical workflow consists of two stages.

### Stage 1: Preparation in Imaginary Time

First, imaginary-time evolution is used to obtain a prepared condensate state:

```bash
python3 bec_simulation.py \
    --mode imaginary \
    --beta 500 \
    --omega 3.0 \
    --potential lattice \
    --V0 15 \
    --kappa1 1.0471975512 \
    --kappa2 1.0471975512 \
    --steps 1000 \
    --output results/ground_state
```

This produces a wave function saved as:

```text
results/ground_state/psi_ground.npy
```

### Stage 2: Real-Time Evolution

Then, this wave function can be used as the initial condition for real-time evolution:

```bash
python3 bec_simulation.py \
    --mode real \
    --load results/ground_state/psi_ground.npy \
    --steps 5000 \
    --dt 0.001 \
    --turn-off-lattice \
    --output results/real_time
```

This stage shows how the vortex and density structures evolve after the potential is modified.

---

## Main Observations

The simulations show how vortex structures depend on physical and numerical parameters.

Typical observations include:

- larger angular velocity can lead to stronger vortex formation,
- the nonlinear interaction strength `beta` changes the size and density profile of the condensate,
- the optical lattice depth `V0` modifies the spatial density structure,
- different values of `kappa1` and `kappa2` change the geometry of the lattice,
- imaginary-time evolution is useful for preparing approximate stationary states,
- real-time evolution is useful for observing later vortex dynamics,
- angular momentum is connected with the rotation of the vortex structure,
- the initially prepared density and vortex geometry can influence the later real-time evolution,
- phase imprinting can be interpreted as preserving information about the initial potential-induced structure in the later condensate dynamics.

---

## Notes on Code Organization

For a clean repository version, the simulation should not store parameters directly in the filename. Instead, the recommended approach is:

1. Provide parameters from the command line.
2. Create a separate output folder for every run.
3. Save all parameters to `config.json`.
4. Save numerical arrays as `.npy` files.
5. Save plots and animations as `.png` and `.gif` files.

This makes the results easier to reproduce and compare.

---

## Bibliography

[1] Y. Zeng, Y. Zhang,  
**Numerical simulation of vortex dynamics in Bose-Einstein condensates**,  
available online: `https://web.mst.edu/~zhangyanz/Papers/Zeng-Zhang-2009.pdf`

[2] A. J. Leggett,  
**Bose-Einstein condensation in the alkali gases: Some fundamental concepts**,  
Reviews of Modern Physics.

[3] L. Pitaevskii, S. Stringari,  
**Bose-Einstein Condensation and Superfluidity**,  
Oxford University Press.

[4] N. P. Proukakis, D. W. Snoke, P. B. Littlewood,  
**Universal themes of Bose-Einstein condensation**,  
Nature Reviews Physics.
