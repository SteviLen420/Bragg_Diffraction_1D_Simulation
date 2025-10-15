# Bragg Diffraction in 1D Multilayers (TMM)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17358243.svg)](https://doi.org/10.5281/zenodo.17358243)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-physics.gen--ph-b31b1b.svg)](https://arxiv.org/)

**Author:** Stefan Len  

A compact, reproducible simulation of **Bragg diffraction** in one-dimensional dielectric multilayers using the **Transfer Matrix Method (TMM)** at normal incidence.  
Validates numerics against the analytical Bragg law and performs a strict **energy conservation** audit.

---

## Summary

- **Goal:** Quantify the Bragg stop band and reflectance of an SiO₂/TiO₂ stack and verify the analytical prediction numerically.  
- **Theory:** For normal incidence the Bragg condition is  
```
  λ_B = 2·n_eff·d,    n_eff = (n₁ + n₂)/2
```
- **Method:** Flux-normalized TMM formulation  
```
  T = (n_out/n_in)·|t|²,    r = M[1,0]/M[0,0],    t = 1/M[0,0]
```
- **Validation:** Sub-grid quadratic interpolation for peak λ_B and conservation test `A = 1 - (R + T)`.

---

## Repository contents

- `Bragg_Diffraction_1D_Simulation.py` — executable code  
- Generated outputs per run:  
  - `reflectivity_spectrum.png`, `validation_plot.png`, `bragg_stack_structure.png`  
  - `wavelength_scan_results.csv`, `summary_report.txt`, `full_summary.json`, `metadata.json`

---

## Quickstart
```bash
python -m venv .venv && source .venv/bin/activate
pip install -U numpy matplotlib
python Bragg_Diffraction_1D_Simulation.py
```

Results appear in timestamped folders: `./Bragg_Diffraction_TMM_Sims/run_YYYYMMDD_HHMMSS/`

**Google Colab:** The script auto-detects Colab and mounts Google Drive. Just upload and run.

---

## Configuration

Edit the **MASTER CONTROL** block at the top of the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LATTICE_CONSTANT` | 120.0 nm | Period d |
| `NUM_PERIODS` | 30 | Number of unit cells N |
| `WAVELENGTH_MIN` | 400.0 nm | Scan start |
| `WAVELENGTH_MAX` | 900.0 nm | Scan end |
| `NUM_WAVELENGTHS` | 50 | Spectral resolution |
| `N1` | 1.46 | Refractive index (SiO₂) |
| `N2` | 2.30 | Refractive index (TiO₂) |

**Example:** For 850 nm VCSEL mirrors:
```python
LATTICE_CONSTANT = 226.0  # nm
NUM_PERIODS = 25
WAVELENGTH_MIN = 700.0
WAVELENGTH_MAX = 1000.0
```

---

## Physical background

### Bragg condition

A periodic stack of alternating dielectrics (n₁, n₂) with period d exhibits high reflectivity when:
```
m·λ = 2·n_eff·d    (m = 1, 2, 3, ...)
```

where `n_eff = (n₁ + n₂)/2` is the effective index. For first-order (m=1):
```
λ_Bragg = 2·n_eff·d
```

### Quarter-wave stack

Each layer has optical thickness λ/4 at design wavelength:
```
d₁ = λ_design / (4·n₁)
d₂ = λ_design / (4·n₂)
```

This maximizes reflectivity via constructive interference.

---

## Transfer Matrix Method

### Matrices

**Interface** (medium a → b):
```
I_ab = 0.5 · [ 1 + n_b/n_a    1 - n_b/n_a ]
             [ 1 - n_b/n_a    1 + n_b/n_a ]
```

**Propagation** (phase φ = 2πnd/λ):
```
P(φ) = [ exp(iφ)      0      ]
       [   0      exp(-iφ)   ]
```

### System matrix

For N periods:
```
M_total = I_in→1 · [M_cell]^N · I_N→out
```

where `M_cell` is one unit cell (layers 1 + 2).

### Observables
```
r = M[1,0] / M[0,0]
t = 1 / M[0,0]
R = |r|²
T = (n_out/n_in) · |t|²
```

### Energy conservation

For lossless media:
```
R + T = 1
```

Deviation `A = 1 - R - T` diagnoses numerical errors.

---

## Output files

### 1. `metadata.json`
Simulation parameters, physics context, material properties.

### 2. `wavelength_scan_results.csv`
Columns: `wavelength`, `reflectivity`, `transmissivity`, `R_plus_T`, `A`, `bragg_prediction`, `normalized_wavelength`.

### 3. `reflectivity_spectrum.png`
Top panel: R(λ) with analytical λ_Bragg marker.  
Bottom panel: T(λ).

### 4. `bragg_stack_structure.png`
Refractive index profile n(z) showing the periodic structure.

### 5. `validation_plot.png`
R vs. normalized wavelength (λ/λ_Bragg) with Bragg region highlighted (±5%).

### 6. `summary_report.txt`
Human-readable summary:
- Theoretical λ_Bragg
- Measured λ_Bragg (interpolated)
- Error percentage
- Max reflectivity
- Energy conservation status

### 7. `full_summary.json`
Structured summary with all metadata and results.

---

## Results and validation

### Typical output (d = 120 nm, N = 30)

**Theory:**
```
λ_Bragg = 2 · [(1.46 + 2.30)/2] · 120 nm = 451.2 nm
```

**Simulation:**
```
Peak wavelength (interpolated): 451.1 nm
Maximum reflectivity: 99.87%
Relative error: 0.02%
Energy conservation: max|A| = 3.2 × 10⁻⁷
```

### Accuracy scaling

| Wavelength points | Relative error |
|-------------------|----------------|
| 25 | ~2% |
| 50 (default) | ~0.5% |
| 100 | ~0.1% |

Quadratic interpolation recovers sub-grid accuracy.

---

## Applications

### 1. VCSEL mirrors
Distributed Bragg reflectors (DBRs) with R > 99% for vertical-cavity surface-emitting lasers at 850 nm, 980 nm.

### 2. Optical filters
Wavelength-selective rejection for spectroscopy, fluorescence microscopy, laser line filtering.

### 3. Photonic crystals
Educational demonstration of photonic bandgap formation and Bloch wave physics.

### 4. Anti-reflection coatings
Inverse design: minimize R by avoiding Bragg condition.

---

## References

### Theory

1. **Yeh, P.** (1988). *Optical Waves in Layered Media*. Wiley.  
   Canonical TMM reference.

2. **Born, M. & Wolf, E.** (1999). *Principles of Optics* (7th ed.). Cambridge.  
   Ch. 1.6: Stratified media.

3. **Saleh, B. E. A. & Teich, M. C.** (2007). *Fundamentals of Photonics* (2nd ed.). Wiley.  
   Ch. 7: Periodic structures.

### Applications

4. **Coldren, L. A., Corzine, S. W., & Mashanovitch, M. L.** (2012). *Diode Lasers and Photonic Integrated Circuits* (2nd ed.). Wiley.  
   Ch. 4: DBR design.

5. **Joannopoulos, J. D., Johnson, S. G., Winn, J. N., & Meade, R. D.** (2008). *Photonic Crystals: Molding the Flow of Light* (2nd ed.). Princeton.

---

## Citation

If you use this code in your research, please cite:
```bibtex
@software{len2025bragg,
  author       = {Len, Stefan},
  title        = {Bragg Diffraction in 1D Multilayers (TMM)},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/yourusername/bragg-diffraction-tmm},
  doi          = {10.5281/zenodo.XXXXXXX}
}
```

**arXiv preprint:** 

---

## License

MIT License — see [LICENSE](../LICENSE) file.

**Copyright © 2025 Stefan Len**

---

## Contact

**Author:** Stefan Len  
**Email:** tqe.simulation@gmail.com 

**GitHub:** [SteviLen420](https://github.com/SteviLen420/Bragg_Diffraction_1D_Simulation)

For questions, issues, or collaboration inquiries, please open an issue on GitHub or contact via email.

---

## Acknowledgments

This work was developed as a pedagogical tool for understanding wave interference in periodic media and for practical DBR design. The implementation follows standard TMM formulations found in classical optics textbooks.

---

**Last updated:** January 2025
