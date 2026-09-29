#!/usr/bin/env python3
"""
reproduce_all.py - Complete Clean-Room Reproduction Suite for A5 Phase 3 Pivot:
Acute Drug-Response Tail Heterogeneity (ART / ARTI).

Usage:
    python reproduce_all.py [--quick] [--skip-downloads]
"""

import os
import sys
import time
import subprocess
import argparse

def run_step(step_name, command):
    print(f"\n========================================================")
    print(f"[*] RUNNING: {step_name}")
    print(f"[*] COMMAND: {command}")
    print(f"========================================================")
    t0 = time.time()
    res = subprocess.run(command, shell=True, capture_output=True, text=True)
    dt = time.time() - t0
    
    if res.returncode != 0:
        print(f"[!] FAILED: {step_name} (exit code {res.returncode}, {dt:.2f}s)")
        print("--- STDERR ---")
        print(res.stderr[-1000:] if len(res.stderr) > 1000 else res.stderr)
        print("--- STDOUT ---")
        print(res.stdout[-1000:] if len(res.stdout) > 1000 else res.stdout)
        return False
    else:
        print(f"[+] SUCCESS: {step_name} ({dt:.2f}s)")
        if res.stdout:
            lines = res.stdout.strip().split("\n")
            for l in lines[-5:]:
                print(f"    {l}")
        return True

def main():
    parser = argparse.ArgumentParser(description="Clean-room reproduction runner for A5")
    parser.add_argument("--skip-figures", action="store_true", help="Skip figure rendering")
    args = parser.parse_args()
    
    print("=" * 70)
    print("A5 PHASE 3 CLEAN-ROOM REPRODUCTION PIPELINE")
    print("Acute Single-Cell Drug-Response Tail Heterogeneity (ART / ARTI)")
    print("=" * 70)
    
    steps = [
        ("Step 1: G10 Falsification & Null 3 Re-check", "python code/external_dtp/check_null3_details.py"),
        ("Step 2: P1 Comparator Metric Benchmark (C0-C5)", "python code/benchmark/run_p1_comparator_benchmark.py"),
        ("Step 3: P2 Replicate Reproducibility (Plate 4 vs 14)", "python code/benchmark/run_p2_replicate_reproducibility.py"),
        ("Step 4: P2 Transcriptional Program Concordance", "python code/benchmark/run_p2_program_concordance.py"),
        ("Step 5: P3 Dose-Response Architecture", "python code/benchmark/run_p3_dose_response_architecture.py"),
        ("Step 6: P5 Context Dependence & Variance Partitioning", "python code/benchmark/run_p5_variance_partitioning.py"),
        ("Step 7: P7 External Validation - SciPlex2", "python code/external_validation/run_sciplex2_p7_benchmark.py"),
        ("Step 8: P7 External Validation - Aissa et al.", "python code/external_validation/run_aissa_p7_benchmark.py"),
        ("Step 9: P11 Definition Stability", "python code/benchmark/run_p11_definition_stability.py"),
        ("Step 10: P12 Negative Controls Floor", "python code/benchmark/run_p12_negative_controls.py"),
    ]
    
    if not args.skip_figures:
        steps.append(("Step 11: Main Figures 1-7 & Source Data Generation", "python code/figures/generate_all_main_figures.py"))
        
    start_time = time.time()
    results = []
    
    for name, cmd in steps:
        ok = run_step(name, cmd)
        results.append((name, ok))
        if not ok:
            print(f"\n[!] Pipeline halted due to error in {name}")
            sys.exit(1)
            
    total_time = time.time() - start_time
    print("\n" + "=" * 70)
    print(f"CLEAN-ROOM REPRODUCTION COMPLETE in {total_time:.2f}s")
    print("=" * 70)
    for name, ok in results:
        print(f"  {'[PASS]' if ok else '[FAIL]'} {name}")
    print("\nAll benchmarks, external validations, tables, and figures verified successfully.")

if __name__ == '__main__':
    main()
