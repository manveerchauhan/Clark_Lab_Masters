# Script to count the frequency of cell barcodes in FASTQ files
# Takes a FASTQ file (gzipped or uncompressed) as input and outputs a CSV file with barcode counts
# To use: python returnBarcodeReadInfo.py input.fastq.gz output.csv

import sys
import gzip
from collections import Counter
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_fastq(filename):
    """Process FASTQ file and count barcodes."""
    barcode_counts = Counter()
    total_reads = 0
    invalid_barcodes = 0
    
    # Function to check if a barcode is valid
    def is_valid_barcode(seq):
        return len(seq) == 16 and all(c in 'ACGT' for c in seq)
    
    # Open file (handles both gzipped and regular files)
    opener = gzip.open if filename.endswith('.gz') else open
    with opener(filename, 'rt') as f:
        # Process 4 lines at a time (FASTQ format)
        while True:
            # Read header line
            header = f.readline().strip()
            if not header:  # End of file
                break
                
            # Skip other three lines of FASTQ entry
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            
            total_reads += 1
            
            # Process header to extract barcode
            if header.startswith('@'):
                barcode = header[1:17]  # Take first 16 chars after @
                if is_valid_barcode(barcode):
                    barcode_counts[barcode] += 1
                else:
                    invalid_barcodes += 1
            else:
                invalid_barcodes += 1
            
            # Log progress every million reads
            if total_reads % 1_000_000 == 0:
                logger.info(f"Processed {total_reads:,} reads...")
    
    logger.info(f"Total reads processed: {total_reads:,}")
    logger.info(f"Valid barcodes found: {sum(barcode_counts.values()):,}")
    logger.info(f"Invalid barcodes: {invalid_barcodes:,}")
    
    return barcode_counts

def main(fastq_file, output_csv):
    logger.info(f"Processing {fastq_file}...")
    
    # Count barcodes
    barcode_counts = process_fastq(fastq_file)
    
    logger.info("Creating DataFrame...")
    # Convert to DataFrame and save
    df = pd.DataFrame(barcode_counts.items(), columns=['barcode', 'count'])
    df = df.sort_values('count', ascending=False)
    
    logger.info(f"Saving to {output_csv}...")
    df.to_csv(output_csv, index=False)
    logger.info("Done!")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fastq[.gz] output.csv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
