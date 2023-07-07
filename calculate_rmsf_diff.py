import sys

def process_pdb(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        alpha_carbons = [line for line in lines if line[13:15].strip() == 'CA']
        beta_values = [float(line[60:66]) for line in alpha_carbons]
        return beta_values, lines

def write_pdb(filename, lines, new_betas):
    with open(filename, 'w') as file:
        beta_index = 0
        for line in lines:
            if line[13:15].strip() == 'CA':
                file.write(line[:60] + "{:6.2f}".format(new_betas[beta_index]) + line[66:])
                beta_index += 1
            else:
                file.write(line[:60] + "{:6.2f}".format(0.0) + line[66:])

def subtract_pdb_beta(pdb1, pdb2, output):
    beta_values1, lines1 = process_pdb(pdb1)
    beta_values2, _ = process_pdb(pdb2)
    new_betas = [bv1 - bv2 for bv1, bv2 in zip(beta_values1, beta_values2)]
    write_pdb(output, lines1, new_betas)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py pdb1 pdb2 output")
        sys.exit(1)

    subtract_pdb_beta(sys.argv[1], sys.argv[2], sys.argv[3])

