import os
import shutil

def write_model_1_file(input_file, output_file):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        lines = []
        for line in fin:
            split_line = line.split()
            if len(split_line) >= 3 and split_line[1] == '1':
                lines.append(line)
        lines.sort(key=lambda x: float(x.split()[0]))
        fout.write(''.join(lines))
    return len([line for line in lines if float(line.split()[0]) <= -9.0])

def modify_and_copy_pdb(src_file, dest_file):
    with open(src_file, 'r') as fin, open(dest_file, 'w') as fout:
        copy_line = False
        for line in fin:
            if line.strip() == 'MODEL 1':
                copy_line = True
            elif line.strip() == 'ENDMDL':
                copy_line = False
            if copy_line:
                fout.write(line)

def modify_and_copy_pdb_rdkit(src_file, dest_file):
    with open(src_file, 'r') as fin, open(dest_file, 'w') as fout:
        copy_line = False
        for line in fin:
            if line.strip() == 'MODEL 1':
                copy_line = True
            elif line.strip() == 'ENDMDL':
                copy_line = False
            if copy_line:
                if line.endswith('A \n'):
                    line = line.rsplit('A ', 1)[0] + 'C \n'
                elif line.endswith('NA\n'):
                    line = line.rsplit('NA', 1)[0] + 'N \n'
                elif line.endswith('OA\n'):
                    line = line.rsplit('OA', 1)[0] + 'O \n'
                elif line.endswith('SA\n'):
                    line = line.rsplit('SA', 1)[0] + 'S \n'
                fout.write(line)

def copy_top_hits(input_file, src_directory, dest_directory, num_hits):
    if not os.path.exists(dest_directory):
        os.mkdir(dest_directory)
    with open(input_file, 'r') as fin:
        for i, line in enumerate(fin):
            if i >= num_hits:
                break
            zinc_number = line.split()[2].split('log_zinc')[-1]
            src_file = os.path.join(src_directory, 'zinc_out_{}.pdb'.format(zinc_number))
            dest_file = os.path.join(dest_directory, 'zinc_hitnumber_{}_{}.pdb'.format(i+1, zinc_number))
            dest_file_rdkit = os.path.join(dest_directory, 'zinc_hitnumber_{}_{}_rdkit.pdb'.format(i+1, zinc_number))
            modify_and_copy_pdb(src_file, dest_file)
            modify_and_copy_pdb_rdkit(src_file, dest_file_rdkit)

def count_lines(filename):
    with open(filename, 'r') as f:
        return len(f.readlines())

num_hits = write_model_1_file('best_hits.dat', 'best_hits_model1.dat')

# Now count the lines in the 'best_hits_model1.dat' file
num_hits = count_lines('best_hits_model1.dat')

copy_top_hits('best_hits_model1.dat', 'RESULTS', 'top_hits100', 100)
copy_top_hits('best_hits_model1.dat', 'RESULTS', 'top_hits1000', 1000)
copy_top_hits('best_hits_model1.dat', 'RESULTS', 'top_hits10', 10)
copy_top_hits('best_hits_model1.dat', 'RESULTS', f'top_hits{num_hits}', num_hits)


