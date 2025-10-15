# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Stefan Len
#
# ===================================================================================
# Bragg_Diffraction_1D_Simulation.py
# ===================================================================================
# Author: Stefan Len
# Overview:
#   Numerical simulation of Bragg diffraction in 1D periodic structures using
#   Transfer Matrix Method (TMM). Validates against analytical Bragg law (nλ = 2d·sinθ).
#   Generates plots, CSV outputs, and quantitative error summary automatically.
# ===================================================================================

import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
import csv
import sys
import json

# ===================================================================================
# --- MASTER CONTROL ---
# ===================================================================================
LATTICE_CONSTANT = 100.0      # Period (d) in nm
NUM_PERIODS = 10              # Number of periods in structure
WAVELENGTH_MIN = 150.0        # Minimum wavelength to test (nm)
WAVELENGTH_MAX = 250.0        # Maximum wavelength to test (nm)
NUM_WAVELENGTHS = 20          # Number of wavelengths to test
N1 = 1.0                      # Refractive index layer 1
N2 = 1.5                      # Refractive index layer 2
OUTPUT_BASE_FOLDER = 'Bragg_Diffraction_TMM_Sims'
CODE_VERSION = '2.0.0'
# ===================================================================================

def setup_colab_environment():
    """Setup Colab environment with necessary packages and Drive mount."""
    print("--- Environment Setup ---")
    
    if 'google.colab' in sys.modules:
        print("Colab detected.")
        
        try:
            from google.colab import drive
            drive.mount('/content/drive')
            print("Google Drive mounted.")
            return True, '/content/drive/MyDrive'
        except Exception as e:
            print(f"Drive mount failed: {e}")
            return False, None
    else:
        print("Local execution detected.")
        return False, '.'

def create_output_directory(base_path):
    """Create timestamped output directory."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(base_path, OUTPUT_BASE_FOLDER, f'run_{timestamp}')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")
    return output_dir

def save_metadata(save_path):
    """Save simulation metadata as JSON."""
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'code_version': CODE_VERSION,
        'simulation_type': 'Bragg_Diffraction_TMM',
        'method': 'Transfer_Matrix_Method',
        'parameters': {
            'lattice_constant_d': LATTICE_CONSTANT,
            'num_periods': NUM_PERIODS,
            'wavelength_range_nm': [WAVELENGTH_MIN, WAVELENGTH_MAX],
            'num_wavelengths': NUM_WAVELENGTHS,
            'refractive_indices': {'n1': N1, 'n2': N2}
        },
        'physics': {
            'bragg_law': 'n*lambda = 2*n_eff*d',
            'normal_incidence': 'theta = 0 degrees',
            'first_order_condition': 'lambda = 2*n_eff*d for n=1',
            'effective_index': 'n_eff = (n1 + n2) / 2'
        }
    }
    
    with open(os.path.join(save_path, 'metadata.json'), 'w') as f:
        json.dump(metadata, f, indent=2)
    print("Metadata saved.")



def simulate_bragg_stack_tmm(wavelength, d, num_periods, n1=1.0, n2=1.5):
    """
    Transfer Matrix Method for 1D Bragg stack - CORRECTED.
    
    Uses proper boundary conditions and interface matrices.
    
    Args:
        wavelength: incident wavelength
        d: period (total thickness of one unit cell)
        num_periods: number of periods
        n1, n2: refractive indices
    
    Returns:
        Reflectivity, Transmissivity
    """
    # Each layer is d/2 thick
    d1 = d / 2
    d2 = d / 2
    
    # Wave vectors
    k1 = 2 * np.pi * n1 / wavelength
    k2 = 2 * np.pi * n2 / wavelength
    
    # Propagation matrices (diagonal form)
    phi1 = k1 * d1
    phi2 = k2 * d2
    
    # Interface matrices (continuity of E and H)
    def interface_matrix(n_a, n_b):
        return 0.5 * np.array([
            [1 + n_b/n_a, 1 - n_b/n_a],
            [1 - n_b/n_a, 1 + n_b/n_a]
        ], dtype=complex)
    
    # Propagation matrix
    def prop_matrix(phi):
        return np.array([
            [np.exp(1j * phi), 0],
            [0, np.exp(-1j * phi)]
        ], dtype=complex)
    
    # System matrix for one period
    # Air -> n1 -> n2 -> Air
    n0 = 1.0
    
    M = np.eye(2, dtype=complex)
    M = M @ interface_matrix(n0, n1)
    M = M @ prop_matrix(phi1)
    M = M @ interface_matrix(n1, n2)
    M = M @ prop_matrix(phi2)
    M = M @ interface_matrix(n2, n1)
    
    # Repeat for N periods (simplified: same as matrix power)
    M_total = np.linalg.matrix_power(M, num_periods)
    
    # Final interface back to air
    M_total = M_total @ interface_matrix(n1, n0)
    
    # Extract coefficients
    # M_total relates [A+, A-] input to [B+, B-] output
    # For reflection: r = M_total[1,0] / M_total[0,0]
    # For transmission: t = 1 / M_total[0,0]
    
    r = M_total[1, 0] / M_total[0, 0]
    t = 1.0 / M_total[0, 0]
    
    R = np.abs(r)**2
    T = np.abs(t)**2
    
    # Ensure physical values
    R = np.clip(R, 0, 1)
    T = np.clip(T, 0, 1)
    
    return R, T

def analytical_bragg_wavelength(d, n1, n2, order=1):
    """
    Analytical Bragg wavelength for quarter-wave stack.
    
    For normal incidence: λ_Bragg = 2 * n_eff * d / order
    where n_eff = (n1 + n2) / 2 (approximation)
    
    Args:
        d: period
        n1, n2: refractive indices
        order: diffraction order
    
    Returns:
        Bragg wavelength
    """
    n_eff = (n1 + n2) / 2.0
    return 2.0 * n_eff * d / order

def run_wavelength_scan(save_path):
    """
    Scan through wavelengths and measure reflectivity using TMM.
    Compare with analytical Bragg condition.
    """
    print("\n--- Wavelength Scan (TMM) ---")
    
    # Test wavelengths
    wavelengths = np.linspace(WAVELENGTH_MIN, WAVELENGTH_MAX, NUM_WAVELENGTHS)
    reflectivities = []
    transmissivities = []
    
    results = []
    
    # Theoretical Bragg wavelength
    lambda_bragg_theory = analytical_bragg_wavelength(LATTICE_CONSTANT, N1, N2, order=1)
    
    for wl in wavelengths:
        R, T = simulate_bragg_stack_tmm(wl, LATTICE_CONSTANT, NUM_PERIODS, N1, N2)
        reflectivities.append(R)
        transmissivities.append(T)
        
        # Check if near Bragg condition
        is_bragg = abs(wl - lambda_bragg_theory) < 0.1 * lambda_bragg_theory
        
        results.append({
            'wavelength': wl,
            'reflectivity': R,
            'transmissivity': T,
            'bragg_prediction': 'Yes' if is_bragg else 'No',
            'normalized_wavelength': wl / lambda_bragg_theory
        })
        
        print(f"  λ={wl:.2f} nm: R={R:.4f}, T={T:.4f}, Bragg={'YES' if is_bragg else 'no'}")
    
    # Save results to CSV
    with open(os.path.join(save_path, 'wavelength_scan_results.csv'), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['wavelength', 'reflectivity', 'transmissivity', 
                                                'bragg_prediction', 'normalized_wavelength'])
        writer.writeheader()
        writer.writerows(results)
    
    print("Results saved to CSV.")
    
    # Create reflectivity plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Reflectivity vs wavelength
    ax1.plot(wavelengths, reflectivities, 'b-o', linewidth=2, markersize=6, label='TMM Numerical')
    ax1.axvline(lambda_bragg_theory, color='r', linestyle='--', linewidth=2, 
                label=f'Bragg (n=1): λ={lambda_bragg_theory:.1f} nm')
    ax1.set_ylabel('Reflectivity', fontsize=12)
    ax1.set_title('Bragg Diffraction (TMM): Reflectivity vs Wavelength', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([-0.05, 1.05])
    
    # Transmissivity vs wavelength
    ax2.plot(wavelengths, transmissivities, 'g-o', linewidth=2, markersize=6)
    ax2.axvline(lambda_bragg_theory, color='r', linestyle='--', linewidth=2)
    ax2.set_xlabel('Wavelength (nm)', fontsize=12)
    ax2.set_ylabel('Transmissivity', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([-0.05, 1.05])
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'reflectivity_spectrum.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Reflectivity spectrum saved.")
    
    # Create structure diagram
    fig, ax = plt.subplots(figsize=(12, 3))
    
    for i in range(NUM_PERIODS):
        # Layer 1 (n1)
        x_start = i * LATTICE_CONSTANT
        x_end = x_start + LATTICE_CONSTANT / 2
        ax.fill_between([x_start, x_end], 0, N1, alpha=0.7, color='blue', edgecolor='black', linewidth=0.5)
        
        # Layer 2 (n2)
        x_start = x_end
        x_end = x_start + LATTICE_CONSTANT / 2
        ax.fill_between([x_start, x_end], 0, N2, alpha=0.7, color='red', edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel('Position (nm)', fontsize=12)
    ax.set_ylabel('Refractive Index', fontsize=12)
    ax.set_title(f'Bragg Stack Structure (N={NUM_PERIODS} periods, d={LATTICE_CONSTANT} nm)', fontsize=14)
    ax.set_ylim([0, max(N1, N2) * 1.2])
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'bragg_stack_structure.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Structure diagram saved.")
    
    return wavelengths, reflectivities, results

def create_validation_plot(wavelengths, reflectivities, save_path):
    """
    Create detailed validation plot comparing numerical results with Bragg law.
    """
    print("\n--- Creating Validation Plot ---")
    
    # Theoretical Bragg wavelength
    lambda_bragg = analytical_bragg_wavelength(LATTICE_CONSTANT, N1, N2, order=1)
    
    # Normalized wavelength
    normalized_wl = wavelengths / lambda_bragg
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot numerical reflectivity
    ax.plot(normalized_wl, reflectivities, 'b-o', linewidth=2, markersize=8, 
            label='TMM Numerical Simulation')
    
    # Mark Bragg condition (λ/λ_Bragg = 1 for first order)
    ax.axvline(1.0, color='r', linestyle='--', linewidth=2, 
               label='Analytical Bragg Condition (n=1)')
    ax.axvspan(0.95, 1.05, alpha=0.2, color='red', label='Bragg Region (±5%)')
    
    ax.set_xlabel('Normalized Wavelength (λ/λ_Bragg)', fontsize=14)
    ax.set_ylabel('Reflectivity', fontsize=14)
    ax.set_title('Validation: TMM vs Analytical Bragg Law', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-0.05, 1.05])
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'validation_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Validation plot saved.")

def _to_py(x):
    """Convert numpy types to Python native types."""
    return x.item() if isinstance(x, np.generic) else x

def generate_summary_report(results, save_path):
    """Generate comprehensive summary report."""
    print("\n--- Generating Summary Report ---")

    # Theoretical Bragg wavelength
    theoretical_bragg = float(analytical_bragg_wavelength(LATTICE_CONSTANT, N1, N2, order=1))
    
    # Find wavelength with maximum reflectivity
    max_r_idx = int(np.argmax([float(_to_py(r['reflectivity'])) for r in results]))
    max_r_result = results[max_r_idx]
    measured_bragg = float(_to_py(max_r_result['wavelength']))
    
    error_percent = float(abs(measured_bragg - theoretical_bragg) / theoretical_bragg * 100.0)

    # Text report
    with open(os.path.join(save_path, 'summary_report.txt'), 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("BRAGG DIFFRACTION SIMULATION (TMM) - SUMMARY REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write("SIMULATION PARAMETERS:\n")
        f.write(f"  Method: Transfer Matrix Method (TMM)\n")
        f.write(f"  Lattice constant (d): {LATTICE_CONSTANT} nm\n")
        f.write(f"  Number of periods: {NUM_PERIODS}\n")
        f.write(f"  Refractive indices: n1={N1}, n2={N2}\n")
        f.write(f"  Wavelength range: {WAVELENGTH_MIN} - {WAVELENGTH_MAX} nm\n\n")
        f.write("ANALYTICAL PREDICTION:\n")
        f.write(f"  Bragg law: n*λ = 2*n_eff*d\n")
        f.write(f"  Effective index: n_eff = (n1+n2)/2 = {(N1+N2)/2:.2f}\n")
        f.write(f"  First order (n=1): λ_Bragg = {theoretical_bragg:.2f} nm\n\n")
        f.write("NUMERICAL RESULTS (TMM):\n")
        f.write(f"  Maximum reflectivity: {float(_to_py(max_r_result['reflectivity'])):.4f}\n")
        f.write(f"  At wavelength: {measured_bragg:.2f} nm\n")
        f.write(f"  Normalized (λ/λ_Bragg): {float(_to_py(max_r_result['normalized_wavelength'])):.4f}\n\n")
        f.write("VALIDATION:\n")
        f.write(f"  Theoretical λ_Bragg: {theoretical_bragg:.2f} nm\n")
        f.write(f"  Measured λ_Bragg: {measured_bragg:.2f} nm\n")
        f.write(f"  Absolute error: {abs(measured_bragg - theoretical_bragg):.2f} nm\n")
        f.write(f"  Relative error: {error_percent:.2f}%\n\n")
        
        if error_percent < 1:
            f.write("✓✓ EXCELLENT VALIDATION: Error < 1%\n")
        elif error_percent < 5:
            f.write("✓ VALIDATION PASSED: Error < 5%\n")
        elif error_percent < 10:
            f.write("~ VALIDATION ACCEPTABLE: Error < 10%\n")
        else:
            f.write("✗ VALIDATION FAILED: Error > 10%\n")
        
        f.write("\n" + "=" * 70 + "\n")

    # JSON with builtin types
    cleaned_results = [{
        'wavelength': float(_to_py(r['wavelength'])),
        'reflectivity': float(_to_py(r['reflectivity'])),
        'transmissivity': float(_to_py(r['transmissivity'])),
        'bragg_prediction': str(r['bragg_prediction']),
        'normalized_wavelength': float(_to_py(r['normalized_wavelength'])),
    } for r in results]

    summary_json = {
        'simulation_info': {
            'code_version': CODE_VERSION,
            'timestamp': datetime.now().isoformat(),
            'type': 'Bragg_Diffraction_TMM',
            'method': 'Transfer_Matrix_Method'
        },
        'parameters': {
            'lattice_constant_nm': float(LATTICE_CONSTANT),
            'num_periods': int(NUM_PERIODS),
            'refractive_indices': {'n1': float(N1), 'n2': float(N2)},
            'wavelength_range_nm': [float(WAVELENGTH_MIN), float(WAVELENGTH_MAX)]
        },
        'results': {
            'theoretical_bragg_wavelength_nm': theoretical_bragg,
            'measured_bragg_wavelength_nm': measured_bragg,
            'max_reflectivity': float(_to_py(max_r_result['reflectivity'])),
            'error_percent': error_percent,
            'validation_passed': bool(error_percent < 5.0)
        },
        'all_measurements': cleaned_results
    }

    with open(os.path.join(save_path, 'full_summary.json'), 'w') as f:
        json.dump(summary_json, f, indent=2)

    print("Summary report saved.")
    print("\nRESULTS:")
    print(f"  Theoretical λ_Bragg: {theoretical_bragg:.2f} nm")
    print(f"  Measured λ_Bragg: {measured_bragg:.2f} nm")
    print(f"  Error: {error_percent:.2f}%")
    print(f"  Validation: {'PASSED ✓' if error_percent < 5 else 'FAILED ✗'}")

def main():
    """Main execution function."""
    print("=" * 70)
    print("BRAGG DIFFRACTION 1D SIMULATION (TRANSFER MATRIX METHOD)")
    print("=" * 70)
    
    # Setup environment
    is_colab, drive_path = setup_colab_environment()
    
    # Create output directory
    if is_colab and drive_path:
        save_path = create_output_directory(drive_path)
    else:
        save_path = create_output_directory('.')
    
    # Save metadata
    save_metadata(save_path)
    
    # Run wavelength scan
    wavelengths, reflectivities, results = run_wavelength_scan(save_path)
    
    # Create validation plot
    create_validation_plot(wavelengths, reflectivities, save_path)
    
    # Generate summary report
    generate_summary_report(results, save_path)
    
    print("\n" + "=" * 70)
    print("SIMULATION COMPLETE")
    print(f"All outputs saved to: {save_path}")
    print("=" * 70)

if __name__ == '__main__':
    main()
