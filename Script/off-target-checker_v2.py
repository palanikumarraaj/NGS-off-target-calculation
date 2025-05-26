#!/usr/bin/env python3
"""
NGS Target Region Analysis Pipeline - Optimized Version

This script analyzes BAM files to count reads and bases in different regions:
1. Reads in target region
2. Reads in extended target region (target + read length buffer)
3. Reads completely outside target regions
4. Base counts for target and untargeted regions

Requirements:
- pysam
- pandas
- numpy

Install with: pip install pysam pandas numpy
"""

import pysam
import pandas as pd
import numpy as np
import argparse
import sys
from collections import defaultdict
import logging
from intervaltree import IntervalTree
import time

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NGSTargetAnalyzer:
    def __init__(self, bam_file, bed_file, read_length):
        """
        Initialize the analyzer with input files and parameters
        
        Args:
            bam_file (str): Path to BAM file
            bed_file (str): Path to BED file with target regions (4 columns, no header)
            read_length (int): Read length for buffer calculations
        """
        self.bam_file = bam_file
        self.bed_file = bed_file
        self.read_length = int(read_length)
        self.target_trees = {}  # IntervalTree for each chromosome
        self.extended_trees = {}  # IntervalTree for extended regions
        self.results = {}
        
    def load_target_regions(self):
        """Load target regions from BED file and create interval trees for fast lookup"""
        try:
            logger.info(f"Loading target regions from {self.bed_file}")
            
            # Read BED file (4 columns: chr, start, end, name)
            bed_df = pd.read_csv(self.bed_file, sep='\t', header=None, 
                               names=['chr', 'start', 'end', 'name'])
            
            logger.info(f"Found {len(bed_df)} regions in BED file")
            
            # Create interval trees for fast overlap detection
            for _, row in bed_df.iterrows():
                chr_name = str(row['chr'])
                start = int(row['start'])
                end = int(row['end'])
                
                # Initialize trees if not exists
                if chr_name not in self.target_trees:
                    self.target_trees[chr_name] = IntervalTree()
                    self.extended_trees[chr_name] = IntervalTree()
                
                # Add target region
                self.target_trees[chr_name][start:end] = 'target'
                
                # Add extended region (with read_length buffer)
                ext_start = max(0, start - self.read_length)
                ext_end = end + self.read_length
                self.extended_trees[chr_name][ext_start:ext_end] = 'extended'
            
            total_regions = sum(len(tree) for tree in self.target_trees.values())
            logger.info(f"Created interval trees for {len(self.target_trees)} chromosomes with {total_regions} total regions")
            
        except Exception as e:
            logger.error(f"Error loading BED file: {e}")
            sys.exit(1)
    
    def calculate_overlap_length(self, chr_name, start, end, trees_dict):
        """Calculate overlap length using interval trees"""
        if chr_name not in trees_dict:
            return 0
        
        overlaps = trees_dict[chr_name][start:end]
        if not overlaps:
            return 0
        
        # Calculate total overlap length
        total_overlap = 0
        for interval in overlaps:
            overlap_start = max(start, interval.begin)
            overlap_end = min(end, interval.end)
            total_overlap += overlap_end - overlap_start
        
        return total_overlap
    
    def analyze_bam(self):
        """Analyze BAM file and count reads in different regions"""
        logger.info(f"Analyzing BAM file: {self.bam_file}")
        
        try:
            # Open BAM file
            bamfile = pysam.AlignmentFile(self.bam_file, "rb")
            
            # Get total number of mapped reads for progress tracking
            try:
                total_mapped = bamfile.mapped
                logger.info(f"Total mapped reads in BAM: {total_mapped:,}")
            except:
                total_mapped = None
                logger.info("Could not determine total read count, progress will be estimated")
            
            # Initialize counters
            stats = {
                'total_reads': 0,
                'mapped_reads': 0,
                'reads_in_target': 0,
                'reads_in_extended_target': 0,
                'reads_completely_untargeted': 0,
                'bases_in_target': 0,
                'bases_in_extended_target': 0,
                'bases_completely_untargeted': 0,
                'total_bases_mapped': 0
            }
            
            # Per-chromosome stats
            chr_stats = defaultdict(lambda: defaultdict(int))
            
            start_time = time.time()
            logger.info("Processing reads...")
            
            for read_count, read in enumerate(bamfile.fetch()):
                # Progress reporting
                if read_count % 250000 == 0 and read_count > 0:
                    elapsed = time.time() - start_time
                    rate = read_count / elapsed
                    if total_mapped:
                        progress = (read_count / total_mapped) * 100
                        eta = (total_mapped - read_count) / rate if rate > 0 else 0
                        logger.info(f"Processed {read_count:,} reads ({progress:.1f}%) - Rate: {rate:.0f} reads/sec - ETA: {eta/60:.1f} min")
                    else:
                        logger.info(f"Processed {read_count:,} reads - Rate: {rate:.0f} reads/sec - Elapsed: {elapsed/60:.1f} min")
                
                stats['total_reads'] += 1
                
                # Skip unmapped reads
                if read.is_unmapped:
                    continue
                    
                stats['mapped_reads'] += 1
                
                # Get read coordinates
                read_chr = str(read.reference_name)
                read_start = read.reference_start
                read_end = read.reference_end
                read_length_actual = read_end - read_start
                
                stats['total_bases_mapped'] += read_length_actual
                chr_stats[read_chr]['total_bases'] += read_length_actual
                chr_stats[read_chr]['total_reads'] += 1
                
                # Check overlaps using interval trees
                target_overlap = read_chr in self.target_trees and bool(self.target_trees[read_chr][read_start:read_end])
                extended_overlap = read_chr in self.extended_trees and bool(self.extended_trees[read_chr][read_start:read_end])
                
                if target_overlap:
                    stats['reads_in_target'] += 1
                    overlap_bases = self.calculate_overlap_length(read_chr, read_start, read_end, self.target_trees)
                    stats['bases_in_target'] += overlap_bases
                    chr_stats[read_chr]['reads_in_target'] += 1
                    chr_stats[read_chr]['bases_in_target'] += overlap_bases
                
                if extended_overlap:
                    stats['reads_in_extended_target'] += 1
                    overlap_bases = self.calculate_overlap_length(read_chr, read_start, read_end, self.extended_trees)
                    stats['bases_in_extended_target'] += overlap_bases
                    chr_stats[read_chr]['reads_in_extended'] += 1
                    chr_stats[read_chr]['bases_in_extended'] += overlap_bases
                
                if not extended_overlap:  # Completely outside extended regions
                    stats['reads_completely_untargeted'] += 1
                    stats['bases_completely_untargeted'] += read_length_actual
                    chr_stats[read_chr]['reads_untargeted'] += 1
                    chr_stats[read_chr]['bases_untargeted'] += read_length_actual
            
            bamfile.close()
            
            elapsed = time.time() - start_time
            logger.info(f"Completed processing {stats['total_reads']:,} total reads ({stats['mapped_reads']:,} mapped) in {elapsed/60:.2f} minutes")
            
            # Calculate additional metrics
            stats['reads_in_extended_only'] = (stats['reads_in_extended_target'] - 
                                             stats['reads_in_target'])
            stats['bases_in_extended_only'] = (stats['bases_in_extended_target'] - 
                                             stats['bases_in_target'])
            
            # Calculate percentages
            if stats['mapped_reads'] > 0:
                stats['pct_reads_in_target'] = (stats['reads_in_target'] / 
                                              stats['mapped_reads'] * 100)
                stats['pct_reads_in_extended'] = (stats['reads_in_extended_target'] / 
                                                stats['mapped_reads'] * 100)
                stats['pct_reads_untargeted'] = (stats['reads_completely_untargeted'] / 
                                               stats['mapped_reads'] * 100)
            
            if stats['total_bases_mapped'] > 0:
                stats['pct_bases_in_target'] = (stats['bases_in_target'] / 
                                              stats['total_bases_mapped'] * 100)
                stats['pct_bases_in_extended'] = (stats['bases_in_extended_target'] / 
                                                stats['total_bases_mapped'] * 100)
                stats['pct_bases_untargeted'] = (stats['bases_completely_untargeted'] / 
                                               stats['total_bases_mapped'] * 100)
            
            self.results = {'global_stats': stats, 'chr_stats': dict(chr_stats)}
            logger.info("Analysis completed successfully")
            
        except Exception as e:
            logger.error(f"Error analyzing BAM file: {e}")
            sys.exit(1)
    
    def calculate_target_region_stats(self):
        """Calculate total target region size and coverage"""
        total_target_size = 0
        total_extended_size = 0
        
        # Calculate total size from interval trees
        for chr_name, tree in self.target_trees.items():
            for interval in tree:
                total_target_size += interval.end - interval.begin
                
        for chr_name, tree in self.extended_trees.items():
            for interval in tree:
                total_extended_size += interval.end - interval.begin
        
        self.results['target_region_size'] = total_target_size
        self.results['extended_region_size'] = total_extended_size
        
        # Calculate coverage
        if total_target_size > 0:
            self.results['target_coverage'] = (self.results['global_stats']['bases_in_target'] / 
                                             total_target_size)
        if total_extended_size > 0:
            self.results['extended_coverage'] = (self.results['global_stats']['bases_in_extended_target'] / 
                                               total_extended_size)
    
    def generate_report(self, output_file=None):
        """Generate a comprehensive report"""
        report = []
        stats = self.results['global_stats']
        
        report.append("=" * 80)
        report.append("NGS TARGET REGION ANALYSIS REPORT")
        report.append("=" * 80)
        report.append(f"BAM File: {self.bam_file}")
        report.append(f"BED File: {self.bed_file}")
        report.append(f"Read Length (buffer): {self.read_length}")
        report.append(f"Number of chromosomes with targets: {len(self.target_trees)}")
        report.append(f"Total target regions: {sum(len(tree) for tree in self.target_trees.values())}")
        report.append(f"Total target region size: {self.results['target_region_size']:,} bp")
        report.append(f"Total extended region size: {self.results['extended_region_size']:,} bp")
        report.append("")
        
        report.append("GLOBAL STATISTICS")
        report.append("-" * 40)
        report.append(f"Total reads: {stats['total_reads']:,}")
        report.append(f"Mapped reads: {stats['mapped_reads']:,}")
        report.append(f"Total mapped bases: {stats['total_bases_mapped']:,}")
        report.append("")
        
        report.append("READ COUNTS BY REGION")
        report.append("-" * 40)
        report.append(f"Reads in target region: {stats['reads_in_target']:,} "
                     f"({stats.get('pct_reads_in_target', 0):.2f}%)")
        report.append(f"Reads in extended target region: {stats['reads_in_extended_target']:,} "
                     f"({stats.get('pct_reads_in_extended', 0):.2f}%)")
        report.append(f"Reads in buffer only (extended - target): {stats['reads_in_extended_only']:,}")
        report.append(f"Reads completely outside target: {stats['reads_completely_untargeted']:,} "
                     f"({stats.get('pct_reads_untargeted', 0):.2f}%)")
        report.append("")
        
        report.append("BASE COUNTS BY REGION")
        report.append("-" * 40)
        report.append(f"Bases in target region: {stats['bases_in_target']:,} "
                     f"({stats.get('pct_bases_in_target', 0):.2f}%)")
        report.append(f"Bases in extended target region: {stats['bases_in_extended_target']:,} "
                     f"({stats.get('pct_bases_in_extended', 0):.2f}%)")
        report.append(f"Bases in buffer only (extended - target): {stats['bases_in_extended_only']:,}")
        report.append(f"Bases completely outside target: {stats['bases_completely_untargeted']:,} "
                     f"({stats.get('pct_bases_untargeted', 0):.2f}%)")
        report.append("")
        
        report.append("COVERAGE STATISTICS")
        report.append("-" * 40)
        report.append(f"Average coverage in target regions: {self.results.get('target_coverage', 0):.2f}x")
        report.append(f"Average coverage in extended regions: {self.results.get('extended_coverage', 0):.2f}x")
        report.append("")
        
        # Per-chromosome breakdown (top 10 chromosomes by read count)
        report.append("TOP CHROMOSOMES BY READ COUNT")
        report.append("-" * 40)
        chr_stats = self.results['chr_stats']
        sorted_chrs = sorted(chr_stats.items(), 
                           key=lambda x: x[1].get('total_reads', 0), 
                           reverse=True)[:10]
        
        for chr_name, chr_data in sorted_chrs:
            report.append(f"{chr_name}:")
            report.append(f"  Total reads: {chr_data.get('total_reads', 0):,}")
            report.append(f"  Reads in target: {chr_data.get('reads_in_target', 0):,}")
            report.append(f"  Reads in extended: {chr_data.get('reads_in_extended', 0):,}")
            report.append(f"  Reads untargeted: {chr_data.get('reads_untargeted', 0):,}")
            report.append(f"  Bases in target: {chr_data.get('bases_in_target', 0):,}")
            report.append(f"  Total bases: {chr_data.get('total_bases', 0):,}")
        
        report_text = "\n".join(report)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_text)
            logger.info(f"Report saved to {output_file}")
        else:
            print(report_text)
        
        return report_text
    
    def save_detailed_results(self, output_prefix):
        """Save detailed results to CSV files"""
        # Global statistics
        global_df = pd.DataFrame([self.results['global_stats']])
        global_df.to_csv(f"{output_prefix}_global_stats.csv", index=False)
        
        # Per-chromosome statistics
        chr_data = []
        for chr_name, stats in self.results['chr_stats'].items():
            row = {'chromosome': chr_name}
            row.update(stats)
            chr_data.append(row)
        
        if chr_data:
            chr_df = pd.DataFrame(chr_data)
            chr_df.to_csv(f"{output_prefix}_chromosome_stats.csv", index=False)
        
        logger.info(f"Detailed results saved with prefix: {output_prefix}")
    
    def run_analysis(self, output_file=None, output_prefix=None):
        """Run the complete analysis pipeline"""
        logger.info("Starting NGS target region analysis...")
        
        # Load target regions
        self.load_target_regions()
        
        # Analyze BAM file
        self.analyze_bam()
        
        # Calculate additional statistics
        self.calculate_target_region_stats()
        
        # Generate report
        self.generate_report(output_file)
        
        # Save detailed results if requested
        if output_prefix:
            self.save_detailed_results(output_prefix)
        
        logger.info("Analysis completed successfully!")
        
        return self.results

def main():
    """Main function to run the pipeline from command line"""
    parser = argparse.ArgumentParser(
        description="Analyze NGS data to count reads and bases in target regions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python ngs_pipeline.py -b sample.bam -t targets.bed -r 150
  python ngs_pipeline.py -b sample.bam -t targets.bed -r 150 -o report.txt
  python ngs_pipeline.py -b sample.bam -t targets.bed -r 150 --output-prefix results
  
BED file format (4 columns, tab-separated, no header):
  chr1    1000    2000    target1
  chr1    3000    4000    target2
        """
    )
    
    parser.add_argument('-b', '--bam', required=True,
                        help='Path to BAM file')
    parser.add_argument('-t', '--target', required=True,
                        help='Path to BED file with target regions (4 columns: chr, start, end, name)')
    parser.add_argument('-r', '--read-length', required=True, type=int,
                        help='Read length for buffer calculations')
    parser.add_argument('-o', '--output', 
                        help='Output file for report (default: print to stdout)')
    parser.add_argument('--output-prefix',
                        help='Prefix for detailed CSV output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate input files
    try:
        with open(args.bam, 'rb'):
            pass
    except FileNotFoundError:
        logger.error(f"BAM file not found: {args.bam}")
        sys.exit(1)
    
    try:
        with open(args.target, 'r'):
            pass
    except FileNotFoundError:
        logger.error(f"BED file not found: {args.target}")
        sys.exit(1)
    
    # Create analyzer and run analysis
    analyzer = NGSTargetAnalyzer(args.bam, args.target, args.read_length)
    results = analyzer.run_analysis(args.output, args.output_prefix)
    
    return results

if __name__ == "__main__":
    main()
