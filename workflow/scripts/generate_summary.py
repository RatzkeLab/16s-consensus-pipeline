#!/usr/bin/env python3
"""
Generate pipeline summary report showing sample attrition.
"""

import sys
from pathlib import Path


def read_checkpoint_summary(path):
    """Read a checkpoint summary TSV and return dict of sample -> status."""
    results = {}
    with open(path) as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                sample, reads, threshold, status = parts[0], parts[1], parts[2], parts[3]
                results[sample] = {'reads': int(reads), 'threshold': int(threshold), 'status': status}
    return results


def main(initial_check, filtered_check, output_path):
    """Generate summary report."""
    
    # Read checkpoint data
    initial = read_checkpoint_summary(initial_check)
    filtered = read_checkpoint_summary(filtered_check)
    
    all_samples = sorted(set(initial.keys()) | set(filtered.keys()))
    
    # Write summary
    with open(output_path, 'w') as f:
        f.write("# 16S Consensus Pipeline Summary Report\n\n")
        
        # Overall statistics
        initial_total = len(initial)
        initial_pass = sum(1 for s in initial.values() if s['status'] == 'PASS')
        filtered_total = len(filtered)
        filtered_pass = sum(1 for s in filtered.values() if s['status'] == 'PASS')
        
        # Get thresholds from first available sample in each dict
        initial_threshold = next(iter(initial.values()))['threshold'] if initial else 'N/A'
        filtered_threshold = next(iter(filtered.values()))['threshold'] if filtered else 'N/A'
        
        f.write("## Pipeline Summary\n\n")
        f.write(f"Total samples processed: {initial_total}\n")
        f.write(f"Passed initial QC (≥{initial_threshold} reads): {initial_pass}/{initial_total}\n")
        f.write(f"Passed filtered QC (≥{filtered_threshold} reads): {filtered_pass}/{filtered_total}\n")
        f.write(f"Final consensus sequences: {filtered_pass}\n\n")
        
        # Sample attrition table
        f.write("## Sample Processing Details\n\n")
        f.write("| Sample | Initial Reads | Initial Status | Filtered Reads | Filtered Status | Final |\n")
        f.write("|--------|--------------|----------------|----------------|-----------------|-------|\n")
        
        for sample in all_samples:
            init_data = initial.get(sample, {})
            filt_data = filtered.get(sample, {})
            
            init_reads = init_data.get('reads', 'N/A')
            init_status = init_data.get('status', 'N/A')
            filt_reads = filt_data.get('reads', 'N/A') if init_status == 'PASS' else '-'
            filt_status = filt_data.get('status', 'N/A') if init_status == 'PASS' else '-'
            
            final = '✓' if filt_status == 'PASS' else '✗'
            
            f.write(f"| {sample} | {init_reads} | {init_status} | {filt_reads} | {filt_status} | {final} |\n")
        
        # Failure reasons
        f.write("\n## Samples Excluded\n\n")
        
        failed_initial = [s for s in all_samples if initial.get(s, {}).get('status') == 'FAIL']
        failed_filtered = [s for s in all_samples if initial.get(s, {}).get('status') == 'PASS' 
                          and filtered.get(s, {}).get('status') == 'FAIL']
        
        if failed_initial:
            f.write(f"### Failed Initial QC ({len(failed_initial)} samples)\n")
            f.write(f"Samples with <{initial_threshold} reads: {', '.join(failed_initial)}\n\n")
        
        if failed_filtered:
            f.write(f"### Failed Post-Filter QC ({len(failed_filtered)} samples)\n")
            f.write(f"Samples with <{filtered_threshold} reads after filtering: {', '.join(failed_filtered)}\n\n")


if __name__ == "__main__":
    # Snakemake script interface
    main(
        snakemake.input.read_summary,
        snakemake.input.filter_summary,
        snakemake.output.report
    )
