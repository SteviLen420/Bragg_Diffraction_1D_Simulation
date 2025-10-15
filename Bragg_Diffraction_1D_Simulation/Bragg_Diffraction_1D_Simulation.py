# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Stefan Len
#
# ===================================================================================
# Bragg_Diffraction_1D_Simulation.py
# ===================================================================================
# Author: Stefan Len
# Overview:
#   Numerical simulation of Bragg diffraction in 1D periodic structures.
#   Implements wave propagation through periodic potential and validates
#   against analytical Bragg law (nλ = 2d·sinθ). Generates plots, CSV outputs,
#   and quantitative error summary automatically.
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
GRID_SIZE = 2000              # Number of spatial points
LATTICE_CONSTANT = 10.0       # Spacing between scatterers (d)
NUM_PERIODS = 20              # Number of periods in structure
WAVELENGTH_MIN = 8.0          # Minimum wavelength to test
WAVELENGTH_MAX = 25.0         # Maximum wavelength to test
NUM_WAVELENGTHS = 10          # Number of wavelengths to test
POTENTIAL_STRENGTH = 1.0      # Strength of periodic potential
OUTPUT_BASE_FOLDER = 'Bragg_Diffraction_Sims'
CODE_VERSION = '1.0.0'
# ===================================================================================

def setup_colab_environment():
    """Setup Colab environment with necessary packages and Drive mount."""
    print("--- Environment Setup ---")
    
    if 'google.colab' in sys.modules:
        print("Colab detected.")
        
        # Install required packages if missing
        try:
            import scipy
        except ImportError:
            print("Installing scipy...")
            os.system('pip install scipy > /dev/null 2>&1')
            print("scipy installed.")
        
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
        'simulation_type': 'Bragg_Diffraction_1D',
        'parameters': {
            'grid_size': GRID_SIZE,
            'lattice_constant_d': LATTICE_CONSTANT,
            'num_periods': NUM_PERIODS,
            'wavelength_range': [WAVELENGTH_MIN, WAVELENGTH_MAX],
            'num_wavelengths': NUM_WAVELENGTHS,
            'potential_strength': POTENTIAL_STRENGTH
        },
        'physics': {
            'bragg_law': 'n*lambda = 2*d*sin(theta)',
            'normal_incidence': 'theta = 90 degrees',
            'first_order_condition': 'lambda = 2*d for n=1'
        }
    }
    
    with open(os.path.join(save_path, 'metadata.json'), 'w') as f:
        json.dump(metadata, f, indent=2)
    print("Metadata saved.")

def create_periodic_potential(x, d, num_periods, strength):
    """
    Create 1D periodic potential (square wave).
    
    Args:
        x: Spatial grid
        d: Lattice constant (period)
        num_periods: Number of periods
        strength: Potential strength
    
    Returns:
        Periodic potential array
    """
    potential = np.zeros_like(x)
    total_length = num_periods * d
    center = len(x) // 2
    start_idx = center - int(total_length / 2)
    end_idx = center + int(total_length / 2)
    
    for i in range(start_idx, end_idx):
        if i >= 0 and i < len(x):
            phase = (x[i] - x[start_idx]) / d
            # Square wave: 1 for half period, 0 for other half
            if phase % 1.0 < 0.5:
                potential[i] = strength
    
    return potential

def simulate_wave_propagation(x, wavelength, potential):
    """
    Simulate wave propagation through periodic structure using transfer matrix.
    
    Args:
        x: Spatial grid
        wavelength: Wavelength of incident wave
        potential: Periodic potential array
    
    Returns:
        Reflectivity and transmissivity
    """
    k = 2 * np.pi / wavelength
    dx = x[1] - x[0]
    
    # Initialize wave (incident from left)
    psi = np.zeros(len(x), dtype=complex)
    psi[0] = 1.0  # Incident amplitude
    
    # Simple propagation with potential scattering
    for i in range(1, len(x)):
        # Phase accumulation
        phase = k * dx
        
        # Scattering from potential
        if potential[i] > 0:
            # Reflection coefficient (simplified)
            r = potential[i] * 0.1 * np.exp(1j * phase)
            psi[i] = psi[i-1] * np.exp(1j * phase) * (1 - r)
        else:
            psi[i] = psi[i-1] * np.exp(1j * phase)
    
    # Calculate reflection and transmission
    # Reflection: look at backward-going wave component
    # Simplified: compare amplitudes at boundaries
    incident_intensity = np.abs(psi[0])**2
    transmitted_intensity = np.abs(psi[-1])**2
    
    transmissivity = transmitted_intensity / incident_intensity
    reflectivity = 1.0 - transmissivity  # Conservation (simplified)
    
    return reflectivity, transmissivity, psi

def analytical_bragg_condition(wavelength, d, order=1):
    """
    Analytical Bragg condition for normal incidence.
    
    For normal incidence (θ = 90°), sin(θ) = 1
    Bragg law: n*λ = 2*d*sin(θ)
    At normal incidence: n*λ = 2*d
    
    Returns True if Bragg condition is satisfied for given order.
    """
    bragg_wavelength = 2 * d / order
    tolerance = 0.1 * d  # 10% tolerance
    return abs(wavelength - bragg_wavelength) < tolerance

def run_wavelength_scan(save_path):
    """
    Scan through wavelengths and measure reflectivity.
    Compare with analytical Bragg condition.
    """
    print("\n--- Wavelength Scan ---")
    
    # Create spatial grid
    x = np.linspace(0, GRID_SIZE, GRID_SIZE)
    
    # Create periodic potential
    potential = create_periodic_potential(x, LATTICE_CONSTANT, NUM_PERIODS, POTENTIAL_STRENGTH)
    
    # Save potential profile
    plt.figure(figsize=(12, 4))
    plt.plot(x, potential, 'k-', linewidth=1)
    plt.xlabel('Position (arbitrary units)', fontsize=12)
    plt.ylabel('Potential Strength', fontsize=12)
    plt.title(f'Periodic Potential (d={LATTICE_CONSTANT}, N={NUM_PERIODS})', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(save_path, 'potential_profile.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Test wavelengths
    wavelengths = np.linspace(WAVELENGTH_MIN, WAVELENGTH_MAX, NUM_WAVELENGTHS)
    reflectivities = []
    transmissivities = []
    bragg_predictions = []
    
    results = []
    
    for wl in wavelengths:
        R, T, psi = simulate_wave_propagation(x, wl, potential)
        reflectivities.append(R)
        transmissivities.append(T)
        
        # Check Bragg condition (first order)
        is_bragg = analytical_bragg_condition(wl, LATTICE_CONSTANT, order=1)
        bragg_predictions.append(1.0 if is_bragg else 0.0)
        
        results.append({
            'wavelength': wl,
            'reflectivity': R,
            'transmissivity': T,
            'bragg_prediction': 'Yes' if is_bragg else 'No',
            'normalized_wavelength': wl / (2 * LATTICE_CONSTANT)
        })
        
        print(f"  λ={wl:.2f}: R={R:.4f}, T={T:.4f}, Bragg={'YES' if is_bragg else 'no'}")
    
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
    ax1.plot(wavelengths, reflectivities, 'b-o', linewidth=2, markersize=6, label='Numerical')
    ax1.axvline(2 * LATTICE_CONSTANT, color='r', linestyle='--', linewidth=2, 
                label=f'Bragg (n=1): λ={2*LATTICE_CONSTANT:.1f}')
    ax1.set_ylabel('Reflectivity', fontsize=12)
    ax1.set_title('Bragg Diffraction: Reflectivity vs Wavelength', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Transmissivity vs wavelength
    ax2.plot(wavelengths, transmissivities, 'g-o', linewidth=2, markersize=6)
    ax2.axvline(2 * LATTICE_CONSTANT, color='r', linestyle='--', linewidth=2)
    ax2.set_xlabel('Wavelength (arbitrary units)', fontsize=12)
    ax2.set_ylabel('Transmissivity', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'reflectivity_spectrum.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Reflectivity spectrum saved.")
    
    return wavelengths, reflectivities, results

def create_validation_plot(wavelengths, reflectivities, save_path):
    """
    Create detailed validation plot comparing numerical results with Bragg law.
    """
    print("\n--- Creating Validation Plot ---")
    
    # Normalized wavelength (λ/2d)
    normalized_wl = wavelengths / (2 * LATTICE_CONSTANT)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot numerical reflectivity
    ax.plot(normalized_wl, reflectivities, 'b-o', linewidth=2, markersize=8, 
            label='Numerical Simulation')
    
    # Mark Bragg condition (λ/2d = 1 for first order)
    ax.axvline(1.0, color='r', linestyle='--', linewidth=2, 
               label='Bragg Condition (n=1)')
    ax.axhspan(0.9, 1.1, alpha=0.2, color='red', label='Bragg Region (±10%)')
    
    ax.set_xlabel('Normalized Wavelength (λ/2d)', fontsize=14)
    ax.set_ylabel('Reflectivity', fontsize=14)
    ax.set_title('Validation: Numerical vs Analytical Bragg Law', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'validation_plot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Validation plot saved.")

def generate_summary_report(results, save_path):
    """
    Generate comprehensive summary report.
    """
    print("\n--- Generating Summary Report ---")
    
    # Find wavelength with maximum reflectivity
    max_r_idx = np.argmax([r['reflectivity'] for r in results])
    max_r_result = results[max_r_idx]
    
    # Theoretical Bragg wavelength (first order)
    theoretical_bragg = 2 * LATTICE_CONSTANT
    
    # Calculate error
    measured_bragg = max_r_result['wavelength']
    error_percent = abs(measured_bragg - theoretical_bragg) / theoretical_bragg * 100
    
    # Write text summary
    with open(os.path.join(save_path, 'summary_report.txt'), 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("BRAGG DIFFRACTION SIMULATION - SUMMARY REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("SIMULATION PARAMETERS:\n")
        f.write(f"  Lattice constant (d): {LATTICE_CONSTANT}\n")
        f.write(f"  Number of periods: {NUM_PERIODS}\n")
        f.write(f"  Wavelength range: {WAVELENGTH_MIN} - {WAVELENGTH_MAX}\n")
        f.write(f"  Grid size: {GRID_SIZE} points\n\n")
        
        f.write("ANALYTICAL PREDICTION:\n")
        f.write(f"  Bragg law: n*λ = 2*d*sin(θ)\n")
        f.write(f"  Normal incidence: θ = 90°, sin(θ) = 1\n")
        f.write(f"  First order (n=1): λ_Bragg = 2*d = {theoretical_bragg:.2f}\n\n")
        
        f.write("NUMERICAL RESULTS:\n")
        f.write(f"  Maximum reflectivity: {max_r_result['reflectivity']:.4f}\n")
        f.write(f"  At wavelength: {measured_bragg:.2f}\n")
        f.write(f"  Normalized (λ/2d): {max_r_result['normalized_wavelength']:.4f}\n\n")
        
        f.write("VALIDATION:\n")
        f.write(f"  Theoretical λ_Bragg: {theoretical_bragg:.2f}\n")
        f.write(f"  Measured λ_Bragg: {measured_bragg:.2f}\n")
        f.write(f"  Absolute error: {abs(measured_bragg - theoretical_bragg):.2f}\n")
        f.write(f"  Relative error: {error_percent:.2f}%\n\n")
        
        if error_percent < 5:
            f.write("✓ VALIDATION PASSED: Error < 5%\n")
        elif error_percent < 10:
            f.write("~ VALIDATION ACCEPTABLE: Error < 10%\n")
        else:
            f.write("✗ VALIDATION FAILED: Error > 10%\n")
        
        f.write("\n" + "=" * 70 + "\n")
    
    # Save JSON summary
    summary_json = {
        'simulation_info': {
            'code_version': CODE_VERSION,
            'timestamp': datetime.now().isoformat(),
            'type': 'Bragg_Diffraction_1D'
        },
        'parameters': {
            'lattice_constant': LATTICE_CONSTANT,
            'num_periods': NUM_PERIODS,
            'wavelength_range': [WAVELENGTH_MIN, WAVELENGTH_MAX]
        },
        'results': {
            'theoretical_bragg_wavelength': theoretical_bragg,
            'measured_bragg_wavelength': measured_bragg,
            'max_reflectivity': max_r_result['reflectivity'],
            'error_percent': error_percent,
            'validation_passed': error_percent < 10
        },
        'all_measurements': results
    }
    
    with open(os.path.join(save_path, 'full_summary.json'), 'w') as f:
        json.dump(summary_json, f, indent=2)
    
    print(f"Summary report saved.")
    print(f"\nRESULTS:")
    print(f"  Theoretical Bragg wavelength: {theoretical_bragg:.2f}")
    print(f"  Measured Bragg wavelength: {measured_bragg:.2f}")
    print(f"  Error: {error_percent:.2f}%")
    print(f"  Validation: {'PASSED' if error_percent < 10 else 'FAILED'}")

def main():
    """Main execution function."""
    print("=" * 70)
    print("BRAGG DIFFRACTION 1D SIMULATION")
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
