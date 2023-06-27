import argparse
from collections import defaultdict

def filter_hbond_counts(input_file, output_file, min_percentage, max_percentage):
    filtered_lines = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            values = line.split()

            # Skip lines that don't have the expected number of values
            if len(values) != 4:
                continue

            donor, hydrogen, acceptor, percentage = values
            percentage = float(percentage)

            if min_percentage <= percentage <= max_percentage:
                filtered_lines.append(line)

    with open(output_file, 'w') as f:
        f.write('\n'.join(filtered_lines))

def filter_hbonds(input_file, output_file, hbonds_counts_file):
    accepted_hbonds = set()
    frame_counts = defaultdict(int)

    with open(hbond_counts_file, 'r') as f:
        for line in f:
            line = line.strip()
            donor, hydrogen, acceptor, _ = line.split()
            accepted_hbonds.add((int(float(donor)), int(float(hydrogen)), int(float(acceptor))))

    filtered_lines = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            frame, acceptor, hydrogen, donor, *_ = line.split()

            if (int(float(donor)), int(float(hydrogen)), int(float(acceptor))) in accepted_hbonds:
                filtered_lines.append(f"{int(float(frame))} {int(float(acceptor))} {int(float(hydrogen))} {int(float(donor))}")
                frame_counts[int(float(frame))] += 1

    with open(output_file, 'w') as f:
        f.write('\n'.join(filtered_lines))

    with open(f"{output_file.split('.')[0]}_counts_by_time.txt", 'w') as f:
        for frame, count in sorted(frame_counts.items()):
            f.write(f"{frame} {count}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter hbond counts and hbonds files.')
    parser.add_argument('prefix', type=str, help='Prefix for the input and output file names')

    args = parser.parse_args()

    prefix = args.prefix
    hbond_counts_file = f"{prefix}_hbond_counts_by_ids.txt"
    input_file = f"{prefix}_hbonds.txt"
    output_hbond_counts_file = f"{prefix}_hbond_filtered_counts_by_ids.txt"
    output_hbonds_file = f"{prefix}_hbonds_filtered.txt"

    min_percentage = 50
    max_percentage = 100

    filter_hbond_counts(hbond_counts_file, output_hbond_counts_file, min_percentage, max_percentage)
    filter_hbonds(input_file, output_hbonds_file, output_hbond_counts_file)

