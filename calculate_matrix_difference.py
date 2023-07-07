import numpy as np
import sys
import seaborn as sns
import matplotlib.pyplot as plt

def subtract_matrices(file1, file2, output_prefix):
    # Load the matrices from the text files
    matrix1 = np.loadtxt(file1)
    matrix2 = np.loadtxt(file2)

    # Check if both matrices have the same shape
    if matrix1.shape != matrix2.shape:
        print("Error: The two matrices don't have the same shape.")
        return

    # Subtract the matrices
    result = matrix1 - matrix2

    # Save the result to the output file
    np.savetxt(output_prefix + '.txt', result)

    # Plot the result as a heatmap
    sns.heatmap(result, cmap='coolwarm',center=0)
    plt.savefig(output_prefix + '.png')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1> <file2> <output_prefix>")
    else:
        subtract_matrices(sys.argv[1], sys.argv[2], sys.argv[3])

