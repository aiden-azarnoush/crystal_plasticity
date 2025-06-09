# Crystal Plasticity Framework

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Theoretical Background](#theoretical-background)
   - [Kinematics](#1-kinematics)
   - [Velocity Gradient](#2-velocity-gradient)
   - [Flow Rule](#3-flow-rule)
   - [Stress-Strain Relation](#4-stress-strain-relation)
   - [Hardening Models](#5-hardening-models)
4. [Numerical Algorithm](#numerical-algorithm)

## Introduction

This repository contains an implementation of an explicit point integration Crystal Plasticity model using the Busso (1990) flow rule modified by Luscher et al. The code simulates the mechanical response of single crystal materials (particularly copper with FCC crystal structure) under large deformation conditions.

Over the past several decades, the study of anisotropic finite plasticity of single crystals has been widely investigated due to its impact in advanced metal forming processes and various industries, particularly microelectronics, which integrate single crystal materials into devices due to their advantageous properties. There is a demand to understand the mechanical and structural response of single crystal materials at the macroscopic/microscopic scale under large deformation.

The repository was originally developed by Aiden Azarnoush (mazarnou@purdue.edu) and includes contributions from collaborators at the School of Engineering of Matter Transport and Energy.

## Features

- **Three Hardening Models:**
  1. Asaro and Needleman
  2. Bassani and Wu
  3. Luscher et al.
- **Customizable Crystal Orientation:** Users can define crystal orientation or select from common directions (<100>, <110>, <310>, <580>)
- **Comprehensive Visualizations:** Generates plots for:
  - Stress vs. strain
  - Resolved shear stress vs. shear strain
  - Shear strain rate vs. time
  - Back stress vs. time

## Theoretical Background

### 1. Kinematics

The formulation is based on the multiplicative decomposition of the deformation gradient into elastic and plastic parts:

$$\mathbf{F} = \mathbf{F}_e \mathbf{F}_p$$

where $\mathbf{F}_p$ represents crystallographic slipping along specific slip systems, and $\mathbf{F}_e$ includes elastic distortion and rigid body rotations.

The Green-Lagrange strain tensor is defined as:

$$\mathbf{E} = \frac{1}{2}(\mathbf{F}^T\mathbf{F} - \mathbf{I})$$

When pushed into the intermediate configuration:

$$\tilde{\mathbf{E}} = \mathbf{F}_p^{-T}\mathbf{E}\mathbf{F}_p^{-1}$$

This can be decomposed additively into elastic and plastic components:

$$\tilde{\mathbf{E}} = \tilde{\mathbf{E}}_e + \tilde{\mathbf{E}}_p$$

where:

$$\tilde{\mathbf{E}}_e = \frac{1}{2}(\mathbf{F}_e^T\mathbf{F}_e - \mathbf{I})$$

$$\tilde{\mathbf{E}}_p = \frac{1}{2}(\mathbf{I} - \mathbf{F}_p^{-T}\mathbf{F}_p^{-1})$$

### 2. Velocity Gradient

The spatial velocity gradient is defined as:

$$\mathbf{L} = \dot{\mathbf{F}}\mathbf{F}^{-1}$$

In the intermediate configuration, it decomposes into:

$$\tilde{\mathbf{L}} = \tilde{\mathbf{L}}_e + \tilde{\mathbf{L}}_p$$

where:

$$\tilde{\mathbf{L}}_e = \mathbf{F}_e^{-1}\dot{\mathbf{F}}_e$$

$$\tilde{\mathbf{L}}_p = \dot{\mathbf{F}}_p\mathbf{F}_p^{-1}$$

The plastic velocity gradient is formulated in terms of slip systems:

$$ \tilde{L_p} = \sum_{\alpha=1}^{N_s} \dot{\gamma}^{\alpha} m_{0}^{\alpha} \otimes n_0^{\alpha}$$


where $\dot{\gamma}^\alpha$ is the slip rate, $N_s$ is the total number of active slip systems (12 for FCC copper), and $\mathbf{m}_0^\alpha$ and $\mathbf{n}_0^\alpha$ are unit normal to the slip planes and slip directions, respectively.

### 3. Flow Rule

The viscoplastic flow rule for the single crystal slip rate is given by:

$$\dot{\gamma}^\alpha = \dot{\gamma_0} \text{sgn}(\tau^\alpha) \exp\left[-\frac{E_0}{kT}\left\langle 1 - \left\langle \frac{|\tau^\alpha| - S_\rho^\alpha}{S_l^\alpha} \right\rangle^{m_1} \right\rangle^{m_2} \right]$$


where:
- $E_0$ is the thermal barrier energy
- $\dot{\gamma}_0$ is the reference strain rate
- $m_1$ and $m_2$ are parameters for the shape of the energy barrier
- $S_l^\alpha$ is the intrinsic lattice resistance to slip
- $S_\rho^\alpha$ is the evolving slip resistance
- $\tau^\alpha$ is the Taylor-Schmid resolved shear stress

The resolved shear stress is defined as:

$$\tau^\alpha = \boldsymbol{\sigma} : (\mathbf{m}^\alpha \otimes \mathbf{n}^\alpha)$$

where $\mathbf{m}^\alpha$ and $\mathbf{n}^\alpha$ are the slip plane normal and slip direction pushed to the current configuration:

$$\mathbf{m}^\alpha = \mathbf{F}_e\mathbf{m}_0^\alpha$$

$$\mathbf{n}^\alpha = \mathbf{F}_e^{-T}\mathbf{n}_0^\alpha$$

### 4. Stress-Strain Relation

The 2nd Piola Kirchhoff stress is related to the Cauchy stress by:

$$\tilde{\mathbf{S}} = \det(\mathbf{F}_e)(\mathbf{F}_e^{-1}\mathbf{\sigma}\mathbf{F}_e^{-T})$$

Using a hyperelastic formulation, the stress-strain relationship is:

$$\tilde{\mathbf{S}} = \frac{\partial \psi}{\partial \tilde{\mathbf{E}}_e} = \overline{\overline{\mathbf{C}}} : \tilde{\mathbf{E}}_e$$

where $\overline{\overline{\mathbf{C}}}$ is the 4th order elasticity tensor reflecting the cubic symmetry of copper, expressed in Voigt notation as:

$$\overline{\overline{\mathbf{C}}} = \begin{bmatrix}
C_{11} & C_{12} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{11} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{12} & C_{11} & 0 & 0 & 0 \\
0 & 0 & 0 & C_{44} & 0 & 0 \\
0 & 0 & 0 & 0 & C_{44} & 0 \\
0 & 0 & 0 & 0 & 0 & C_{44}
\end{bmatrix}$$

### 5. Hardening Models

The evolution of slip resistance is defined as:

$$\dot{S_\rho^\alpha} = \sum_{\beta=1}^{N_s} h_{\alpha\beta}|\dot{\gamma}^\beta|$$

where $h_{\alpha\beta}$ is the hardening matrix or instantaneous strain hardening moduli. Three different hardening formulations are implemented:

#### a. Asaro and Needleman

The hardening matrix is defined as:

$$h_{\alpha\beta} = h(\gamma^\alpha)q_{\alpha\beta}$$

where the self-hardening function is:

$$h(\gamma^\alpha) = h_s + (h_0 - h_s)\text{sech}^2\left(\frac{h_0-h_s}{\tau_s-\tau_0}\bar{\gamma}\right)$$

with $\bar{\gamma} = \sum_{\alpha=1}^{N_s}|\gamma^\alpha|$

The latent hardening matrix $q_{\alpha\beta}$ for an FCC crystal is:

$$q_{\alpha\beta} = \begin{bmatrix}
A & qA & qA & qA \\
qA & A & qA & qA \\
qA & qA & A & qA \\
qA & qA & qA & A
\end{bmatrix}$$

where $q=1.4$ for copper and $A$ is a 3×3 matrix of ones.

#### b. Bassani and Wu

The diagonal components of the hardening matrix are:

$$h_{\alpha\alpha} = h(\gamma^\alpha)G$$

with the interactive hardening term:

$$G = 1 + \sum_{\alpha \neq \beta}f_{\alpha\beta}\tanh\left(\frac{\gamma^\beta}{\gamma_0}\right)$$

and a modified self-hardening function:

$$h(\gamma^\alpha) = h_s + (h_0 - h_s)\text{sech}^2\left(\frac{h_0-h_s}{\tau_s-\tau_0}\gamma^\alpha\right)$$

The off-diagonal components are:

$$h_{\alpha\beta} = h_{\alpha\alpha}q_{\alpha\beta}, \quad \alpha \neq \beta \quad \text{(no sum on } \alpha \text{)}$$

For best fits, $q_{\alpha\beta}$ is often set to zero, resulting in a diagonal hardening matrix.

#### c. Luscher et al.

The hardening modulus depends on the evolution of saturation stress:

$$h_{\alpha\beta} = h_0[r + (1-r)\delta_{\alpha\beta}]\left(\frac{S_s^\beta-S_\rho^\beta}{S_s^\beta-S_0^\beta}\right)$$

where the saturation stress is:

$$S_s^\beta = S_{s0}\left|\frac{\dot{\gamma}^\beta}{\dot{\gamma}_0}\right|^{kT/A}$$

with $S_0^\beta$ as the initial slip system resistance, $h_0$ as initial self-hardening, $r$ as the slip system hardening interaction coefficient, $S_{s0}$ as saturation stress at 0K, and $A$ as activation energy.

## Numerical Algorithm

The implementation uses an explicit numerical integration scheme. The key steps are:

1. Update the plastic deformation gradient using the exponential map:
   $$\mathbf{F}_p^{t+\Delta t} = \exp(\tilde{\mathbf{L}}_p^t \Delta t)\mathbf{F}_p^t$$

2. Compute the plastic velocity gradient:
   $$\tilde{L_p}^t = \sum_{\alpha=1}^{N_s}\dot{\gamma}_t^\alpha m_0^\alpha \otimes n_0^\alpha$$

3. Update the slip resistance and shear strain:
   $$S_{p,t+\Delta t}^\alpha = S_{p,t}^\alpha + \dot{S_{p,t}}^\alpha \Delta t$$
   $$\gamma_{t+\Delta t}^\alpha = \gamma_t^\alpha + \dot{\gamma}_t^\alpha \Delta t$$

4. Update the elastic deformation gradient:
   $$F_{t+\Delta t}^e = F_{t+\Delta t}(F_{t+\Delta t}^p)^{-1}$$

This algorithm can be used for any form of the deformation gradient, though the current implementation focuses on uniaxial strain where:

$$\mathbf{F} = \begin{bmatrix}
\lambda(t) & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

with $\lambda$ increasing linearly from 1 to 1.1.

## Convergence Analysis

To ensure accuracy, convergence has been investigated in two ways:

1. **Visual Comparison**: Plotting results with increasing time steps to check for consistency and reduction in oscillations.

2. **Point-to-Point Comparison**: Comparing results at specific points for different numbers of steps.

The minimum number of steps required for convergence with a tolerance of 10^-7 are shown in the table below:

| **Direction** | **Asaro and Needleman** | **Bassani and Wu** | **Luscher et al.** |
|---------------|-------------------------|--------------------|--------------------|
| <100>         | 2758                    | 2750               | 2758               |
| <110>         | 3300                    | 3500               | 2538               |
| <310>         | 2500                    | 2998               | 2578               |
| <580>         | 2542                    | 2972               | 2538               |

Based on these results, 10^4 steps (dt = 10^-12 seconds) is sufficient for all models and orientations.

## Results and Discussion

### B1. Cauchy Stress vs. Logarithmic Strain

The stress-strain curves for all hardening models and crystal orientations are relatively linear and similar to each other. This is because the elastic deformation gradient is relatively the same in all directions at a given time, as the plastic deformation dominates the total deformation. With a linear constitutive equation in terms of the elastic deformation gradient, the stress becomes relatively linear despite slip dislocations.

### B2. Shear Strain Rate Evolution Equation

Examining the shear strain rate of each slip system over time reveals which slip systems are activated. Notably, for any hardening model and any direction, only eight slip systems are activated, with four remaining inactive. The inactive slip systems depend on the crystal orientation:

| **Crystal orientation** | **Inactive slip systems** |
|-------------------------|---------------------------|
| <100>                   | B5, C5, D6, A6            |
| <110>                   | B2, C1, D1, A2            |
| <310>                   | B5, C3, D4, A6            |
| <580>                   | B2, C1, D1, A2            |

This behavior is related to the FCC crystal structure of copper. When copper is under uniaxial strain along certain directions, slip systems with a Schmid factor of zero (where resolved shear stress and shear strain become zero) remain inactive.

### B3. Back Flow Stress

The back stress increases for all hardening models due to the evolution equation:

$$\dot{S_\rho}^\alpha = \sum_{\beta=1}^{N_s} h_{\alpha\beta}|\dot{\gamma}^\beta|$$

Although shear strain rates change, their multiplication with the hardening moduli leads to increasing back stress. The differences between hardening models are evident in the back stress evolution:

- **Asaro-Needleman**: Similar flow stress for all slip systems due to its co-block 12×12 hardening matrix
- **Bassani and Wu**: Flow stress for non-active slip systems remains zero due to zero latent hardening
- **Luscher et al.**: Dynamic hardening dependent on slip rates

### B4. Hardening Models

- **Asaro and Needleman**: The hardening moduli are constant for all 12 slip systems.
- **Bassani and Wu**: Hardening moduli depend on accumulated slip $\gamma^\alpha$. This model is considered more physically motivated and captures the effect of crystal orientation more accurately.
- **Luscher et al.**: Hardening moduli depend on slip rates $\dot{\gamma}^\alpha$. Because slip rates don't change significantly in the simulations, this model's behavior is somewhat similar to the Asaro model.

## Usage

Run the main script to start the simulation:

```python
crystal_plasticity_python.py
```

## Usage

Run the main script to start the simulation. 

The user will be prompted to:
1. Select crystal orientation (or specify a custom one)

```text
   Crystal orientation:
   1- <100>
   2- <110>
   3- <310>
   4- <580>
   5- other
   > 
```

1. Choose the hardening model to use

```text
   Which hardening model:
   1- Asaro and Needleman hardening model
   2- Bassani and Wu hardening model
   3- Luscher et al. hardening model
   > 
``` 

## Output

The program generates four plots:
1. Stress vs. strain curves (normalized by copper yield stress, σ₀ = 117 MPa)
2. Resolved shear stress vs. shear strain 
3. Shear strain rate vs. time
4. Back stress vs. time

Each plot includes data for all 12 slip systems, allowing detailed analysis of slip system activation and interaction.

## References

1. Pierce, D., Asaro, R.J., Needleman, A., "Material rate dependence and localized deformation in crystalline solids," *Acta metallurgica*, Vol. 31, 1951-1976. 1983
2. Asaro, R.J., Needleman, A., "Texture development and strain hardening in rate dependent polycrystals," *Acta metallurgica*, Vol. 33, 923-953. 1985
3. Bassani, J.L., Wu, T., "Latent Hardening in Single Crystals II. Analytical Characterization and Predictions," *Mathematical and Physical Sciences*, Vol. 453, 21-41. 1991
4. Luscher, D.J., Bronkhorst, C.A., Alleman, C.N., Addessio, F.L., "A model for finite-deformation nonlinear thermomechanical response of a single crystal copper under shock conditions," *Journal of Mechanics and Physics of Solids*, Vol. 61, 1877-1894. 2013
5. Wu, P.D., Neale, K.W., Van der Giessen, E., "Simulation of the behavior of FCC polycrystals during reversed torsion," *International Journal of Plasticity*, Vol. 12, 1199-1219. 1996
6. Hill, R., Rice, JR., "Constitutive analysis of elastic-plastic crystals at arbitrary strain," *Journal of Mechanics and Physics of Solids*, Vol. 20, 401-413. 1972
7. Hull, D., Bacon, D.J., *Introduction to Dislocations* Elsevier, Oxford, 2011
8. Belytscho, T. Wing, K.L., Moran, B., and Elkhodary, K.I., *Nonlinear Finite Elements for Continua and Structures*, 2nd ed., Wiley, West Sussex, UK, 2014

## Author
Aiden Azarnoush (mazarnou@purdue.edu)

---

*Dedicated to my mother - Simin Nematpour*
