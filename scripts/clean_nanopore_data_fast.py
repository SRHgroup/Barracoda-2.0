import os
import gzip
import regex
import argparse
import zipfile
import shutil
from Bio import SeqIO

def detect_file_type(filepath):
    """Detects file type using bytes. Returns 'zip', 'gz', or 'other'."""
    if os.path.isdir(filepath): # find directory 
        return "dir"
    with open(filepath, "rb") as f:
        signature = f.read(4)
    if signature.startswith(b"\x50\x4B\x03\x04"):  # find .zip
        return "zip"
    elif signature.startswith(b"\x1F\x8B"):  # find .gz 
        return "gz"
    else:
        return "other" 

def reverse_complement(seq):
    complement = str.maketrans('ATCGN', 'TAGCN')
    return seq.translate(complement)[::-1]

def fuzzyMapFind(subSeq, sequence, subs=2):
    pattern = f"(?e)({subSeq}){{s<={subs}}}"
    match = regex.search(pattern, sequence)
    return -1 if match is None else match.start()

def get_quality_score(quality):
    phred_scores = [ord(char) - 33 for char in quality]
    return sum(phred_scores) / len(phred_scores) if phred_scores else 0

def process_fastq(input_path, output_file, junkfile, primerA, primerB,
                  min_length=132, max_length=136, quality_threshold=10, n_subs=2,
                  extract_dir=None):
    """
    Handles:
      - a single fastq.gz file
      - a directory containing fastq.gz files
      - a zip file containing fastq.gz files (unzipped to extract_dir or outfile dir)
    """
    total_reads, passed_reads = 0, 0

    filetype = detect_file_type(input_path)

    # If zip, extract it first
    if filetype == "zip":
        # if extract_dir not provided, use outfile directory
        if extract_dir is None:
            extract_dir = os.path.dirname(os.path.abspath(output_file))
        os.makedirs(extract_dir, exist_ok=True)
        with zipfile.ZipFile(input_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        input_path = extract_dir  # Now treat as a directory

    with open(output_file, 'a') as final_out, open(junkfile, 'a') as final_junk:
        if os.path.isdir(input_path):
            for root, _, files in os.walk(input_path):
                for fname in files:
                    if detect_file_type(os.path.join(root, fname)) == "gz":
                        filepath = os.path.join(root, fname)
                        print(f"Processing {filepath}...")
                        tr, pr = trim_fastq(filepath, final_out, final_junk,
                                            primerA, primerB,
                                            min_length, max_length,
                                            quality_threshold, n_subs)
                        total_reads += tr
                        passed_reads += pr
        elif filetype == "gz":
            print(f"Processing {input_path}...")
            total_reads, passed_reads = trim_fastq(input_path, final_out, final_junk,
                                                   primerA, primerB,
                                                   min_length, max_length,
                                                   quality_threshold, n_subs)
        else:
            raise ValueError("Invalid input: Must be a .fastq.gz file, .zip file, or a directory containing fastq.gz files.")
    
    percent_passed = (passed_reads / total_reads * 100) if total_reads > 0 else 0
    print(f"Total reads processed: {total_reads}")
    print(f"Reads passed: {passed_reads} ({percent_passed:.2f}%)")

def trim_fastq(input_file, final_out, final_junk, primerA, primerB,
               min_length=132, max_length=136, quality_threshold=10, n_subs=2):
    total_reads, passed_reads = 0, 0
    primerA_reverse = reverse_complement(primerA)
    primerB_reverse = reverse_complement(primerB)
    
    infile = gzip.open(input_file, 'rt') if detect_file_type(input_file) == "gz" else open(input_file, 'r')
    with infile:
        while True:
            header = infile.readline().strip()
            sequence = infile.readline().strip()
            plus_line = infile.readline().strip()
            quality = infile.readline().strip()
            
            if not header:
                break
            
            total_reads += 1
            primerA_index = fuzzyMapFind(primerA, sequence, subs=n_subs)
            primerB_reverse_index = fuzzyMapFind(primerB_reverse, sequence, subs=n_subs)
            primerB_index = fuzzyMapFind(primerB, sequence, subs=n_subs)
            primerA_reverse_index = fuzzyMapFind(primerA_reverse, sequence, subs=n_subs)
            
            if (primerA_index != -1 and primerB_reverse_index != -1) and (primerA_index < primerB_reverse_index):
                trimmed_sequence = sequence[primerA_index-10:primerB_reverse_index+len(primerB)]
                trimmed_quality = quality[primerA_index-10:primerB_reverse_index+len(primerB)]
                
                if (min_length <= len(trimmed_sequence) <= max_length) and (get_quality_score(trimmed_quality[:10]) > quality_threshold):
                    final_out.write(f"{header}\n{trimmed_sequence}\n{plus_line}\n{trimmed_quality}\n")
                    passed_reads += 1
                else:
                    final_junk.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")

            elif (primerB_index != -1 and primerA_reverse_index != -1) and (primerB_index < primerA_reverse_index):
                trimmed_sequence = sequence[primerB_index:primerA_reverse_index+len(primerA_reverse)+10]
                trimmed_quality = quality[primerB_index:primerA_reverse_index+len(primerA_reverse)+10]
                reversed_trimmed_sequence = reverse_complement(trimmed_sequence)
                reversed_trimmed_quality = trimmed_quality[::-1]
                
                if (min_length <= len(trimmed_sequence) <= max_length) and (get_quality_score(reversed_trimmed_quality[:10]) > quality_threshold):
                    final_out.write(f"{header}\n{reversed_trimmed_sequence}\n{plus_line}\n{reversed_trimmed_quality}\n")
                    passed_reads += 1
                else:
                    final_junk.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
            else:
                final_junk.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
    
    return total_reads, passed_reads

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Cleans Nanopore data for Barracoda.")
    parser.add_argument("--infile", "-i", type=str, required=True, help="Input file (.fastq.gz, .zip, or directory)")
    parser.add_argument("--outfile", "-o", type=str, required=True, help="Final output file with trimmed reads")
    parser.add_argument("--junkfile", "-j", type=str, required=True, help="Final junk file with discarded reads")
    parser.add_argument("--primerA", "-a", type=str, default='GAAGTTCCAGCCAGCGTCACAGTTT', help="Forward primer A from Barcoding")
    parser.add_argument("--primerB", "-b", type=str, default='CAATCTTGAGCGTGACTTAAG', help="Forward primer B from Barcoding")
    parser.add_argument("--n_subs", "-n", type=int, default=2, help="Number of substitutions allowed in the forward primer A and B")
    parser.add_argument("--quality", "-q", type=int, default=10, help="Quality threshold for the A-key found within the read")
    parser.add_argument("--min_length", "-min", type=int, default=132, help="Minimum length of reads to keep (default: 132)")
    parser.add_argument("--max_length", "-max", type=int, default=136, help="Maximum length of reads to keep (default: 136)")
    parser.add_argument("--extract_dir", "-x", type=str, help="Directory to extract .zip input into")
    
    args = parser.parse_args()
    
    process_fastq(args.infile, args.outfile, args.junkfile, args.primerA, args.primerB,
                  min_length=args.min_length, max_length=args.max_length,
                  quality_threshold=args.quality, n_subs=args.n_subs,
                  extract_dir=args.extract_dir)
    
    print(f"Processing complete. Results saved to {args.outfile}")

