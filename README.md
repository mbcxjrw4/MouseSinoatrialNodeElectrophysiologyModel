# MouseSinoatrialNodeElectrophysiologyModel
A set of mathematical models simulating the electrophysiological activity of the mouse sinoatrial node (SAN) and surrounding atrial tissue, from single-cell (0D) level to 2D tissue level.
This work supported published research on Sick Sinus Syndrome ([TBX18 overexpression enhances pacemaker function in a rat subsidiary atrial pacemaker model of sick sinus syndrome](https://pubmed.ncbi.nlm.nih.gov/30259525/)).

## 🧪 Overview

- 0D Single-cell model (mouse0D.c): Simulates ion channel dynamics and action potentials in individual SAN cells.

- 1D Tissue model (mouse1D.c): Simulates electrical conduction along a 1D fiber of SAN tissue.

- 2D Tissue models (mouse2D.c, mouse2D_operator_splitting.c): Simulate electrical conduction and pacemaking in a 2D SAN–atrium tissue sheet. Includes an operator-splitting method for efficient stiff ODE integration.

## ⚡ Features
- Modular, extensible C implementation using object-oriented design
  - `cell` base structure with `san_cell` and `atrium_cell` child structures
  
- Parallel simulation using OpenMP

- ODE solving via CVODE (from SUNDIALS)

- Modular, extensible C implementation

- Configurable initial conditions and tissue geometries

## 📁 Repository Structure
```bash
project_root/ 
├── FILE/ 
│   ├── 2D/                 # 2D tissue geometry 
│   └── initial_values/     # initial state values 
│ 
├── Scripts/ 
│   ├── mouse0D.c            # single cell model 
│   ├── mouse1D.c            # 1D tissue model 
│   ├── mouse2D.c             # 2D tissue model 
│   └── mouse2D_operator_splitting.c # 2D model with operator splitting 
│ 
├── LIBRARY/ 
│   ├── constants.h 
│   ├── for_cvode.h 
│   ├── initial_value.h 
│   └── single_cell.h 
│ 
├── Makefile 
└── README.md 
```
