# Numerical Validation of Bragg Diffraction in One-Dimensional Photonic Crystals using Transfer Matrix Method

**Stefan Len**

*Independent Researcher*

**Date:** October 15, 2025

---

## Abstract

I present a rigorous numerical implementation of the Transfer Matrix Method (TMM) for simulating electromagnetic wave propagation in one-dimensional periodic dielectric structures. The simulation validates the Bragg diffraction condition for normal incidence through quantitative comparison with theoretical predictions. Using a representative SiO₂/TiO₂ multilayer stack with 30 periods, equal physical layer thicknesses of 60 nm, and total period d = 120 nm, I achieve excellent agreement between numerical and analytical results, with a relative error of 0.040% in the Bragg wavelength determination. Energy conservation is verified to machine precision (|A| < 10⁻¹⁴), confirming the numerical stability of my implementation. The computed reflectivity spectrum exhibits the characteristic photonic bandgap with R > 99.99% at the Bragg wavelength λ_B = 451.02 nm, along with Fabry-Pérot oscillations outside the stop band. This work provides a validated computational framework for designing distributed Bragg reflectors (DBRs) used in vertical-cavity surface-emitting lasers (VCSELs) and optical filter applications.

**Keywords:** Bragg diffraction, Transfer Matrix Method, photonic crystals, distributed Bragg reflector, numerical simulation, VCSEL

---

## 1. Introduction

### 1.1 Physical Motivation

One-dimensional photonic crystals, consisting of alternating layers of dielectric materials with different refractive indices, exhibit unique optical properties that have enabled transformative technologies in optoelectronics. When the periodicity of such structures approaches the wavelength of light, constructive interference of reflected waves leads to the formation of a photonic bandgap—a range of wavelengths for which propagation through the structure is forbidden [1,2]. This phenomenon, known as Bragg diffraction, is the fundamental principle behind distributed Bragg reflectors (DBRs), which serve as high-reflectivity mirrors in vertical-cavity surface-emitting lasers (VCSELs) [3], narrow-band optical filters [4], and other photonic devices.

### 1.2 Theoretical Background

For a periodic multilayer stack with alternating layers of optical thickness n₁d₁ and n₂d₂, the Bragg condition for constructive interference at normal incidence is:

$$m\lambda = 2(n_1 d_1 + n_2 d_2)$$

where m is the diffraction order (m = 1, 2, 3, ...) and λ is the vacuum wavelength.

For the specific case of **equal physical thickness** (d₁ = d₂ = d/2), where d is the period, this simplifies to:

$$m\lambda = d(n_1 + n_2)$$

For first-order diffraction (m = 1), we can define an effective refractive index:

$$n_{\text{eff}} = \frac{n_1 + n_2}{2}$$

which yields the familiar Bragg wavelength expression:

$$\lambda_B = 2n_{\text{eff}}d$$

**Important note**: This n_eff definition (arithmetic mean) applies specifically to structures with equal physical thicknesses (d₁ = d₂). For quarter-wave stack designs where the optical thicknesses are equal (n₁d₁ = n₂d₂ = λ_B/4), the appropriate effective index is the geometric mean: n_eff = √(n₁·n₂). These represent two distinct design philosophies with different optimization criteria.

### 1.3 Computational Approach

The Transfer Matrix Method (TMM) is a powerful analytical tool for computing the optical response of stratified media [5,6]. By representing the electromagnetic field propagation through each layer and interface as matrix operations, the TMM enables exact computation of reflection and transmission coefficients for multilayer structures. This method is particularly well-suited for Bragg reflectors, as it naturally accounts for multiple internal reflections and interference effects.

### 1.4 Objectives

The primary objectives of this work are:

- To implement a numerically stable TMM algorithm with proper boundary conditions
- To validate the numerical implementation against analytical Bragg law predictions for equal-thickness designs
- To verify energy conservation in the lossless limit
- To characterize the photonic bandgap structure of a representative DBR system
- To provide a reproducible computational framework for DBR design

---

## 2. Theory and Method

### 2.1 Transfer Matrix Formulation

For **normal incidence (θ = 0°)**, the TE and TM polarizations are degenerate, and the following formulation applies to both. The electromagnetic field in a stratified medium can be represented by forward and backward propagating waves. At any position, the field amplitudes are related by:

$$(E^+) = M \begin{pmatrix} E_0^+ \\ E_0^- \end{pmatrix}$$

where M is the transfer matrix connecting the input and output fields.

#### 2.1.1 Interface Matrix

When light encounters an interface between media with refractive indices n_a and n_b, the transfer matrix accounting for Fresnel reflection and transmission is:

$$I_{ab} = \frac{1}{t_{ab}}\begin{pmatrix} 1 & r_{ab} \\ r_{ab} & 1 \end{pmatrix}$$

where the Fresnel coefficients for normal incidence are:

$$r_{ab} = \frac{n_a - n_b}{n_a + n_b}, \quad t_{ab} = \frac{2n_a}{n_a + n_b}$$

This can be written explicitly as:

$$I_{ab} = \frac{1}{2}\begin{pmatrix} 1 + \frac{n_b}{n_a} & 1 - \frac{n_b}{n_a} \\ 1 - \frac{n_b}{n_a} & 1 + \frac{n_b}{n_a} \end{pmatrix}$$

This matrix ensures proper field continuity at the interface.

#### 2.1.2 Propagation Matrix

Propagation through a homogeneous layer of thickness d with refractive index n introduces a phase shift:

$$\phi = \frac{2\pi nd}{\lambda}$$

The propagation matrix is:

$$P(\phi) = \begin{pmatrix} e^{i\phi} & 0 \\ 0 & e^{-i\phi} \end{pmatrix}$$

#### 2.1.3 Unit Cell Matrix

For a periodic structure with two layers per unit cell, the transfer matrix for one complete period is:

$$M_{\text{cell}} = P(\phi_1) \cdot I_{12} \cdot P(\phi_2) \cdot I_{21}$$

where subscripts 1 and 2 denote the two materials.

#### 2.1.4 Total System Matrix

For N periods, the total transfer matrix is:

$$M_{\text{total}} = I_{\text{in}\to 1} \cdot [M_{\text{cell}}]^N \cdot I_{N\to\text{out}}$$

where the entrance and exit interface matrices account for the surrounding medium (typically air).

### 2.2 Reflection and Transmission Coefficients

The complex reflection and transmission amplitudes are obtained from the total matrix elements:

$$r = \frac{M_{10}}{M_{00}}, \quad t = \frac{1}{M_{00}}$$

The power reflectivity R and transmissivity T, accounting for energy flux through interfaces, are:

$$R = |r|^2, \quad T = \frac{n_{\text{out}}}{n_{\text{in}}}|t|^2$$

The factor n_out/n_in ensures proper energy flux normalization when the refractive indices of input and output media differ.

### 2.3 Energy Conservation

For lossless dielectric media with real refractive indices, energy conservation requires:

$$R + T = 1$$

Any deviation from unity indicates numerical errors or unphysical absorption. I define the absorption/error parameter:

$$A = 1 - (R + T)$$

For a numerically stable implementation, I require |A| < 10⁻¹⁰ across the entire spectral range.

### 2.4 Numerical Implementation

#### 2.4.1 Wavelength Scan

I compute R(λ) and T(λ) by scanning wavelength from λ_min to λ_max with N_λ points. At each wavelength:

1. Compute phase shifts φ₁(λ) and φ₂(λ)
2. Construct unit cell matrix M_cell(λ)
3. Compute M_total(λ) = I_in→1 · [M_cell(λ)]^N · I_N→out
4. Extract r(λ) and t(λ) from matrix elements
5. Calculate R(λ), T(λ), and verify A(λ) ≈ 0

#### 2.4.2 Peak Wavelength Determination

To achieve sub-grid accuracy in determining λ_B, I employ quadratic interpolation around the reflectivity maximum. Given three points (λ_i, R_i) near the peak, I fit a parabola:

$$R(\lambda) = a\lambda^2 + b\lambda + c$$

and solve for the vertex position:

$$\lambda_{\text{peak}} = -\frac{b}{2a}$$

This technique typically improves accuracy by an order of magnitude compared to simple maximum finding.

---

## 3. System Configuration

### 3.1 Material Selection

I simulate a SiO₂/TiO₂ multilayer stack, a common material combination for visible-wavelength DBRs:

- **Layer 1:** Silicon dioxide (SiO₂), n₁ = 1.46 (at ~450 nm)
- **Layer 2:** Titanium dioxide (TiO₂), n₂ = 2.30 (at ~450 nm)

These materials offer:

- Large refractive index contrast (Δn = 0.84)
- Excellent optical transparency in the visible range
- Well-established thin-film deposition techniques
- Thermal and chemical stability

**Note on dispersion:** Refractive indices exhibit wavelength dependence n(λ), but we use constant values as a first approximation. For high-precision applications, wavelength-dependent Sellmeier equations should be incorporated.

### 3.2 Geometric Parameters

The structure consists of N = 30 periods with equal physical layer thicknesses:

- **Period:** d = 120 nm
- **Layer 1 thickness:** d₁ = 60 nm
- **Layer 2 thickness:** d₂ = 60 nm
- **Total stack thickness:** L = N·d = 3600 nm = 3.6 µm

This **equal-thickness design** (d₁ = d₂) simplifies fabrication and yields an effective refractive index:

$$n_{\text{eff}} = \frac{n_1 + n_2}{2} = \frac{1.46 + 2.30}{2} = 1.88$$

**Note:** This arithmetic mean definition is valid specifically for equal physical thicknesses. A quarter-wave stack design (n₁d₁ = n₂d₂ = λ_B/4) would use the geometric mean n_eff = √(n₁·n₂) ≈ 1.832 instead, representing a different optimization approach.

### 3.3 Analytical Prediction

From the Bragg condition for equal-thickness layers, I predict first-order maximum reflectivity at:

$$\lambda_B = 2n_{\text{eff}}d = 2 \times 1.88 \times 120\,\text{nm} = 451.2\,\text{nm}$$

This wavelength lies in the blue region of the visible spectrum.

### 3.4 Simulation Parameters

- **Wavelength range:** 400–900 nm
- **Number of wavelength points:** N_λ = 50
- **Spectral resolution:** Δλ ≈ 10.2 nm
- **Input/output medium:** Air (n_in = n_out = 1.0)
- **Incidence angle:** θ = 0° (normal incidence)

---

## 4. Results

### 4.1 Refractive Index Profile

![Refractive Index Profile](figs/bragg_stack_structure.png)

**Figure 1:** Refractive index profile of the 30-period SiO₂/TiO₂ multilayer stack. Blue regions (n = 1.46) represent SiO₂ layers, red regions (n = 2.30) represent TiO₂ layers. The periodic structure extends over 3600 nm with uniform 60 nm layer thicknesses.

Figure 1 shows the one-dimensional refractive index profile n(z) of the simulated structure. The sharp discontinuities at each interface give rise to partial reflections that interfere constructively at the Bragg wavelength. The regularity and uniformity of the structure are critical for achieving high reflectivity over a narrow spectral range.

### 4.2 Reflectivity and Transmissivity Spectra

![Reflectivity and Transmissivity](figs/bragg_reflectivity.png)

**Figure 2:** (Top) Reflectivity spectrum showing the photonic bandgap centered at λ_B = 451.2 nm. The blue curve represents TMM numerical results; the red dashed line indicates the analytical Bragg wavelength prediction. (Bottom) Corresponding transmissivity spectrum demonstrating complementary behavior (R + T = 1).

Figure 2 presents the central result: the wavelength-dependent reflectivity R(λ) and transmissivity T(λ). Several key features are evident:

#### 4.2.1 Photonic Bandgap (Stop Band)

In the range 400–500 nm, the structure exhibits extraordinarily high reflectivity (R > 0.999), with the peak occurring at λ ≈ 451 nm. This "stop band" represents the photonic bandgap where light propagation is forbidden. The numerical results show:

- **Maximum reflectivity:** R_max = 0.9999999999717146 ≈ 1.0 (99.99999999717%)
- **Stop band width (FWHM):** Δλ ≈ 50 nm
- **Spectral selectivity:** Δλ/λ_B ≈ 11%

The near-unity reflectivity with 30 periods demonstrates the efficiency of coherent interference in periodic structures.

#### 4.2.2 Band Edge

At λ ≈ 500 nm, the reflectivity drops sharply from R ≈ 1 to R ≈ 0.2 over ~20 nm. This abrupt transition defines the band edge, beyond which light can propagate through the structure. The steepness of this edge is determined by the refractive index contrast and number of periods.

#### 4.2.3 Fabry-Pérot Oscillations

For λ > 500 nm, the reflectivity exhibits periodic oscillations with gradually decreasing amplitude. These Fabry-Pérot resonances arise from multiple reflections between the front and back interfaces of the finite stack. The oscillation period decreases with wavelength, and the envelope decays as 1/N, where N is the number of periods.

#### 4.2.4 Complementary Transmission

The transmissivity spectrum (Figure 2, bottom) is perfectly complementary to the reflectivity: where R ≈ 1, T ≈ 0, and vice versa. Within the stop band, virtually no light is transmitted (T < 10⁻⁹). Outside the stop band, transmission ranges from 40% to nearly 100%, modulated by the Fabry-Pérot resonances.

### 4.3 Validation Against Analytical Theory

![Validation](figs/bragg_validation.png)

**Figure 3:** Validation of numerical simulation against analytical Bragg law. Reflectivity plotted versus normalized wavelength λ/λ_B. The red dashed line marks the first-order Bragg condition (λ/λ_B = 1), and the pink shaded region indicates ±5% tolerance. TMM results (blue) show excellent agreement with theory.

Figure 3 presents the validation analysis by plotting reflectivity against the normalized wavelength λ/λ_B. This dimensionless representation allows direct comparison with the universal Bragg condition and reveals several important features:

#### 4.3.1 Primary Bragg Peak

The dominant reflectivity maximum occurs precisely at λ/λ_B = 0.9996, corresponding to λ = 451.02 nm. The extremely close alignment with the theoretical prediction (λ/λ_B = 1.0, marked by the red dashed line) validates my numerical implementation. The peak position was refined using quadratic interpolation around the maximum, achieving sub-grid accuracy of ±0.2 nm despite the 10.2 nm sampling interval.

#### 4.3.2 Fabry-Pérot Resonances Beyond the Stop Band

Additional reflectivity maxima are visible at λ/λ_B ≈ 1.3, 1.6, 1.9, etc. These are **NOT higher-order Bragg peaks** (which would occur at λ/λ_B = 1/m < 1 for orders m = 2, 3, ...). Instead, they are **Fabry-Pérot resonances** arising from the finite thickness of the multilayer stack, representing standing wave modes within the structure.

**True higher-order Bragg diffraction** would occur at:
- m=2: λ₂ = λ_B/2 ≈ 225 nm (UV, outside our simulation range)
- m=3: λ₃ = λ_B/3 ≈ 150 nm (deep UV, outside our simulation range)

These shorter-wavelength higher-order peaks are not captured in the 400-900 nm spectral window investigated here.

#### 4.3.3 Bragg Region

The pink shaded region (0.95 < λ/λ_B < 1.05) encompasses wavelengths within 5% of the Bragg condition. The reflectivity remains above 99% throughout this region, indicating the robustness of the photonic bandgap against small wavelength variations.

### 4.4 Quantitative Validation Metrics

**Table 1:** Comparison of theoretical predictions with TMM numerical results.

| Quantity | Theoretical | Numerical | Error |
|----------|-------------|-----------|-------|
| Bragg wavelength λ_B (nm) | 451.20 | 451.02 | 0.18 nm |
| Relative error in λ_B | — | — | 0.040% |
| Maximum reflectivity R_max | 1.0 (ideal) | 0.999999999972 | 2.8×10⁻¹¹ |
| Stop band width Δλ (nm) | ~45 (estimated) | ~50 | ~10% |

The exceptional agreement (0.040% error in λ_B) validates both the physical model and the numerical implementation. The sub-0.1% accuracy is well within the requirements for practical DBR design.

### 4.5 Energy Conservation Verification

A critical test of numerical accuracy is verification of energy conservation. For all wavelength points, energy conservation is satisfied to machine precision:

$$R(\lambda) + T(\lambda) = 1 \pm \varepsilon_{\text{machine}}$$

Key findings:

- **Maximum deviation:** |A|_max = 2.58 × 10⁻¹⁴
- **Mean deviation:** ⟨|A|⟩ = 5.2 × 10⁻¹⁵
- **All points satisfy:** |A| < 10⁻¹³

These values are at the level of machine precision for double-precision floating-point arithmetic (ε ≈ 2.2 × 10⁻¹⁶), confirming the numerical stability of my TMM implementation. No unphysical absorption or numerical artifacts are present.

---

## 5. Discussion

### 5.1 Physical Interpretation

The exceptional agreement between theory and simulation (0.040% error) confirms that the Transfer Matrix Method correctly captures the physics of Bragg diffraction in one-dimensional photonic crystals. The near-unity reflectivity (R > 99.99%) at the Bragg wavelength arises from constructive interference of waves reflected from the 60 interfaces in the 30-period stack. Each interface contributes a small reflection (~4% for the n₁/n₂ contrast), but phase coherence causes these reflections to add constructively at λ_B, yielding total reflection.

The photonic bandgap width (Δλ ≈ 50 nm) is determined primarily by the refractive index contrast. Larger Δn produces wider bandgaps and steeper band edges. The choice of SiO₂/TiO₂ (Δn = 0.84) provides a good balance between bandgap width and material compatibility.

### 5.2 Comparison with Infinite Stack Theory

For an infinite periodic structure (N → ∞), the Bragg condition predicts complete reflection (R = 1) at λ_B with zero bandwidth. My finite stack exhibits:

- Slightly broadened stop band (~50 nm vs. δ-function)
- Fabry-Pérot oscillations outside the gap
- Non-zero transmission at band edges

These deviations are characteristic of finite-size effects and diminish as N increases. With 30 periods, I am in the regime where the structure behaves nearly as an ideal Bragg reflector within the stop band, but finite-size effects are visible outside it.

### 5.3 Equal-Thickness vs. Quarter-Wave Design

My structure uses **equal physical thickness** (d₁ = d₂ = 60 nm), which simplifies fabrication but differs from the optimal **quarter-wave stack** design.

**Quarter-wave stack** (optimal for maximum reflectivity):
- Condition: n₁d₁ = n₂d₂ = λ_B/4
- Layer 1: d₁ = λ_B/(4n₁) = 451.2/(4×1.46) ≈ 77.2 nm
- Layer 2: d₂ = λ_B/(4n₂) = 451.2/(4×2.30) ≈ 49.0 nm
- Effective index: n_eff = √(n₁·n₂) ≈ 1.832
- Bragg wavelength: λ_B = 2·n_eff·(d₁+d₂)/2 ≈ 463 nm

**Equal-thickness design** (my approach):
- Condition: d₁ = d₂ = 60 nm
- Effective index: n_eff = (n₁+n₂)/2 = 1.88
- Bragg wavelength: λ_B = 2·n_eff·d₁ = 451.2 nm
- Advantage: Simpler fabrication (single thickness control)
- Trade-off: ~5-10% lower peak reflectivity than optimal quarter-wave design at longer wavelengths

The two designs target different wavelengths for the same total period and use different effective index definitions. The equal-thickness approach prioritizes manufacturing simplicity, while the quarter-wave approach maximizes reflectivity for a given number of periods.

### 5.4 Practical Implications for DBR Design

The validated simulation framework enables optimization of DBR structures for specific applications:

**VCSEL Mirrors**

For 850 nm VCSELs (common in optical data links), my results suggest:
- Required period: d ≈ 226 nm
- 25–30 periods sufficient for R > 99.9%
- AlGaAs/GaAs material system (Δn ≈ 0.5) requires ~40 periods for equivalent reflectivity

**Optical Filters**

For narrow-band filtering:
- Increase N to sharpen band edges (steepness ∝ N)
- Reduce Δn for narrower stop bands (Δλ/λ_B ∝ Δn)
- Graded interfaces reduce sidelobes

**Broadband Reflectors**

For wide stop bands:
- Maximize Δn (e.g., Si/SiO₂: Δn ≈ 2.6)
- Use chirped structures (varying d along stack)
- Stack multiple Bragg reflectors at different λ_B

### 5.5 Numerical Accuracy and Limitations

My implementation achieves 0.040% error in λ_B determination, limited primarily by:

1. **Wavelength sampling:** Δλ = 10.2 nm discrete grid
   - Mitigated by quadratic interpolation
   - Could be reduced to <0.01% with Δλ = 1 nm

2. **Effective index approximation:** n_eff = (n₁+n₂)/2 for equal-thickness design
   - Exact within the equal-thickness framework
   - Different approximation needed for quarter-wave stacks

3. **Normal incidence assumption:** θ = 0°
   - Angular dependence: λ_B(θ) = λ_B(0)√(1 - sin²θ/n²_eff)
   - Extension to oblique incidence: straightforward

4. **Material dispersion neglected:** constant n₁, n₂
   - Effect: <1% for narrow spectral ranges
   - Sellmeier equations available for high precision

The energy conservation verification (|A| < 10⁻¹⁴) confirms that numerical round-off errors are negligible. Matrix exponentiation ([M_cell]^N) is numerically stable for N ≤ 100, beyond which alternative methods (e.g., Bloch mode decomposition) may be preferable.

### 5.6 Extensions and Future Work

Potential extensions of this work include:

- **Oblique incidence:** Angle-dependent reflectivity R(λ, θ) and polarization-dependent effects (TE vs. TM)
- **Absorption:** Complex refractive indices n = n' + in'' for lossy materials
- **Defect modes:** Engineered defects create transmission resonances within the bandgap
- **Nonlinear effects:** Intensity-dependent n for optical switching applications
- **Chirped structures:** Spatially varying d(z) for broadband response
- **Two-dimensional photonic crystals:** Lateral periodicity for full photonic bandgaps
- **Thermal tuning:** Temperature-dependent n(T) for tunable filters

---

## 6. Conclusions

I have presented a comprehensive numerical study of Bragg diffraction in one-dimensional photonic crystals using the Transfer Matrix Method. The key findings are:

1. **Excellent validation:** The numerical simulation reproduces the analytical Bragg wavelength for equal-thickness designs with 0.040% error, demonstrating the accuracy of the TMM implementation.

2. **Energy conservation:** All numerical results satisfy R + T = 1 to within 10⁻¹⁴, confirming the physical consistency and numerical stability of the algorithm.

3. **Photonic bandgap characterization:** The SiO₂/TiO₂ stack with N = 30 periods achieves 99.99999999717% reflectivity at λ_B = 451.02 nm, with a stop band width of approximately 50 nm.

4. **Finite-size effects:** The simulation correctly captures Fabry-Pérot oscillations outside the stop band, validating the method's ability to describe complex interference phenomena in finite structures.

5. **Design framework:** The distinction between equal-thickness and quarter-wave stack designs is clarified, providing guidance for application-specific optimization.

6. **Sub-grid accuracy:** Quadratic interpolation enables determination of the Bragg wavelength to 0.2 nm precision despite 10 nm spectral sampling.

7. **Practical applicability:** The validated framework provides a reliable tool for designing DBRs for VCSELs, optical filters, and other photonic devices.

This work demonstrates that the Transfer Matrix Method, when carefully implemented with proper boundary conditions and flux normalization, provides an exact and efficient approach for simulating wave propagation in stratified media. The complete agreement with analytical theory, combined with rigorous energy conservation, establishes this as a trustworthy computational framework for photonic crystal design and analysis.

The methodology and open-source implementation provided here enable researchers and engineers to:

- Design custom DBR structures for specific wavelengths and applications
- Optimize layer thicknesses and material combinations
- Predict optical performance before fabrication
- Understand the interplay between structure and optical response

Future applications of this framework to more complex geometries, oblique incidence, and nonlinear effects will further extend its utility in the rapidly advancing field of photonics.

---

## Acknowledgments

The author thanks the open-source scientific Python community for providing the numerical tools (NumPy, Matplotlib) that enabled this work. This research was conducted independently without external funding.

---

## References

[1] Yeh, P. (1988). *Optical Waves in Layered Media*. John Wiley & Sons, New York.

[2] Joannopoulos, J. D., Johnson, S. G., Winn, J. N., & Meade, R. D. (2008). *Photonic Crystals: Molding the Flow of Light* (2nd ed.). Princeton University Press.

[3] Coldren, L. A., Corzine, S. W., & Mashanovitch, M. L. (2012). *Diode Lasers and Photonic Integrated Circuits* (2nd ed.). John Wiley & Sons.

[4] Macleod, H. A. (2010). *Thin-Film Optical Filters* (4th ed.). CRC Press.

[5] Born, M., & Wolf, E. (1999). *Principles of Optics* (7th ed.). Cambridge University Press.

[6] Saleh, B. E. A., & Teich, M. C. (2007). *Fundamentals of Photonics* (2nd ed.). John Wiley & Sons.

---

## Appendix A: Computational Details

### A.1 Software Implementation

The simulation was implemented in Python 3.7+ using:

- NumPy 1.21+ for numerical linear algebra
- Matplotlib 3.4+ for visualization
- Standard library modules for I/O and data management

The complete source code is available at: [GitHub repository URL]

### A.2 Algorithm Pseudocode
```python
for each wavelength λ in [λ_min, λ_max]:
    # Compute phase shifts
    φ₁ = 2π·n₁·d₁/λ
    φ₂ = 2π·n₂·d₂/λ
    
    # Construct matrices
    I₁₂ = interface_matrix(n₁, n₂)
    I₂₁ = interface_matrix(n₂, n₁)
    P₁ = propagation_matrix(φ₁)
    P₂ = propagation_matrix(φ₂)
    
    # Unit cell
    M_cell = P₁ · I₁₂ · P₂ · I₂₁
    
    # Total stack
    M_stack = [M_cell]^N
    M_total = I_in→1 · M_stack · I_N→out
    
    # Reflection and transmission
    r = M_total[1,0] / M_total[0,0]
    t = 1 / M_total[0,0]
    R[λ] = |r|²
    T[λ] = (n_out/n_in)·|t|²
    
    # Energy check
    A[λ] = 1 - (R[λ] + T[λ])
```

### A.3 Computational Performance

- **Single wavelength point:** ~0.5 ms (Intel Core i7, 3.6 GHz)
- **Full 50-point scan:** ~25 ms
- **Memory usage:** <10 MB
- **Matrix exponentiation:** O(log N) using repeated squaring

The algorithm scales efficiently: doubling N increases computation time by ~10%, doubling N_λ doubles computation time linearly.

---

## Appendix B: Data Availability

All simulation data, including:

- Raw reflectivity and transmissivity values (CSV format)
- Metadata and parameters (JSON format)
- High-resolution figures (PNG, 300 DPI)
- Complete simulation summary (TXT format)

are archived with DOI: [Zenodo DOI](https://doi.org/10.5281/zenodo.17358243) and available at the associated GitHub repository.

---

**Manuscript Version:** 1.0  
**Word Count:** ~5,200  
**Figures:** 3  
**Tables:** 1  
**Code Availability:** [Bragg_Diffraction_1D_Simulation](https://github.com/SteviLen420/Bragg_Diffraction_1D_Simulation)  
**Data Availability:** [Zenodo DOI](https://doi.org/10.5281/zenodo.17358243) 

---

*Correspondence:* Stefan Len, tqe.simulation@gmail.com, 
[GitHub: @SteviLen420](https://github.com/SteviLen420/Bragg_Diffraction_1D_Simulation)
 
*Date:* October 15, 2025
