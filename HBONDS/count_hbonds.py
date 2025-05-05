import argparse
from collections import Counter
import math

def analyze_hbond_counts(input_file, output_file, frame_step):
    # Read the data and parse each line by the first nine fields and frame number (field 10)
    with open(input_file, 'r') as f:
        # Skip header line
        header = f.readline().strip()
        data = []
        frame_numbers = []

        # Parse each line for the hydrogen bond fields and frame number
        for line in f:
            fields = line.strip().split()
            hbond_data = tuple(fields[:9])  # First nine fields
            frame_number = int(fields[9])   # Tenth field as frame number
            data.append(hbond_data)
            frame_numbers.append(frame_number)

    # Count occurrences of each unique hydrogen bond by the first nine fields
    hbond_counter = Counter(data)
    total_frames = math.ceil((max(frame_numbers) + 1) / frame_step)  # Total frames, rounded up

    # Prepare results sorted by frequency (high to low)
    results = hbond_counter.most_common()

    # Write results to output file, including the percentage of frames for each hydrogen bond
    with open(output_file, 'w') as f:
        f.write("Donor_resid Donor_resname Donor_name Hydrogen_resid Hydrogen_resname Hydrogen_name "
                "Acceptor_resid Acceptor_resname Acceptor_name Count Percent_of_Frames\n")
        for fields, count in results:
            percent_of_frames = (count / total_frames) * 100
            f.write(" ".join(fields) + f" {count} {percent_of_frames:.2f}%\n")

def main():
    parser = argparse.ArgumentParser(description='Count occurrences of hydrogen bonds.')
    parser.add_argument('input_file', help='Path to the input file with hydrogen bond data.')
    parser.add_argument('output_file', help='Path to the output file with results.')
    parser.add_argument('--frame_step', type=int, default=1, help='Step size between frames.')
    args = parser.parse_args()

    analyze_hbond_counts(args.input_file, args.output_file, args.frame_step)

if __name__ == '__main__':
    main()

