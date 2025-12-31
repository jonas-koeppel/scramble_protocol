#!/usr/bin/env python3
"""
Count insertion sequences in amplicon sequencing data
"""

import argparse
import re
import pandas as pd
import gzip
import glob
from pathlib import Path
from collections import Counter

def read_fastq(file_path):
    """Read sequences from FASTQ file"""
    open_func = gzip.open if str(file_path).endswith('.gz') else open
    sequences = []
    
    with open_func(file_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence lines
                sequences.append(line.strip())
    return sequences

def count_insertions(sequences, insertion_seq):
    """Count sequences with and without insertion""" 
    pattern = r"TTCTACCAGAGGTACA(.*?)AGG.GG"
    
    # Extract matching sequences
    matches = [re.search(pattern, seq) for seq in sequences]
    guides = [m.group() for m in matches if m]
    
    # Calculate expected lengths based on insertion sequence
    wildtype_length = 22  # Base length without any insertion
    insertion_length = wildtype_length + len(insertion_seq)
    
    # Count by length
    lengths = Counter(len(g) for g in guides)
    
    n_wildtype = lengths.get(wildtype_length, 0)
    n_insertion = lengths.get(insertion_length, 0)
    n_other = sum(count for length, count in lengths.items() 
                 if length not in [wildtype_length, insertion_length])
    
    total_analyzed = n_wildtype + n_insertion + n_other
    perc_insertion = (n_insertion / (n_wildtype + n_insertion) * 100 
                     if (n_wildtype + n_insertion) > 0 else 0)
    
    return {
        'total_reads': len(sequences),
        'guide_matches': len(guides),
        'n_wildtype': n_wildtype,
        'n_insertion': n_insertion,
        'n_other': n_other,
        'total_analyzed': total_analyzed,
        'perc_insertion': round(perc_insertion, 2)
    }

def main():
    parser = argparse.ArgumentParser(description='Count insertion sequences in amplicon data')
    parser.add_argument('-f', '--fastq', required=True, 
                       help='Input FASTQ file(s) - can use wildcards (e.g., folder/*.fastq.gz)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file (e.g., insertion_count.csv)')
    parser.add_argument('-i', '--insertion', 
                       default='TAACTTCGTATAATGTACATTATACGAAGTTATG',
                       help='Insertion sequence (default: 34bp LINE1 sequence)')
    
    args = parser.parse_args()
    
    # Find all matching FASTQ files
    if '*' in args.fastq or '?' in args.fastq:
        fastq_files = glob.glob(args.fastq)
        if not fastq_files:
            print(f"No files found matching: {args.fastq}")
            return
    else:
        fastq_files = [args.fastq]
    
    print(f"Found {len(fastq_files)} FASTQ file(s) to process")
    print(f"Insertion sequence: {args.insertion} ({len(args.insertion)} bp)")
    print(f"Looking for: wildtype={22}bp, insertion={22 + len(args.insertion)}bp")
    
    # Process all files
    all_results = []
    for fastq_file in sorted(fastq_files):
        print(f"\nProcessing: {Path(fastq_file).name}")
        
        # Read sequences
        sequences = read_fastq(fastq_file)
        print(f"  {len(sequences)} sequences")
        
        # Count insertions
        results = count_insertions(sequences, args.insertion)
        
        # Create output record
        sample_name = Path(fastq_file).stem.split('.')[0]
        output_data = {'sample': sample_name, **results}
        all_results.append(output_data)
        
        # Print summary for this file
        print(f"  Wildtype: {results['n_wildtype']}, Insertions: {results['n_insertion']}, Rate: {results['perc_insertion']}%")
    
    # Save combined results
    df = pd.DataFrame(all_results)
    df.to_csv(args.output, index=False)
    
    # Print overall summary
    total_reads = df['total_reads'].sum()
    total_insertions = df['n_insertion'].sum()
    total_wildtype = df['n_wildtype'].sum()
    overall_rate = total_insertions / (total_insertions + total_wildtype) * 100 if (total_insertions + total_wildtype) > 0 else 0
    
    print(f"\n=== OVERALL SUMMARY ===")
    print(f"Files processed: {len(all_results)}")
    print(f"Total reads: {total_reads:,}")
    print(f"Total insertions: {total_insertions:,}")
    print(f"Overall insertion rate: {overall_rate:.2f}%")
    print(f"Results saved to: {args.output}")

if __name__ == "__main__":
    main() 