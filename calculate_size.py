import math
import argparse

def calculate_size(tolerance, cutoffdistance, box_length):
    alpha = (1.0/cutoffdistance)*math.sqrt(-math.log(2.0*tolerance))
    size = int(math.ceil(2*alpha*box_length/(3*pow(tolerance, 0.2))))
    return size

def main():
    parser = argparse.ArgumentParser(description='Calculate size.')
    parser.add_argument('-t', '--tolerance', type=float, required=True, help='Input tolerance.')
    parser.add_argument('-c', '--cutoffdistance', type=float, required=True, help='Input cutoff distance.')
    parser.add_argument('-b', '--boxlength', type=float, required=True, help='Input box length.')
    args = parser.parse_args()

    size = calculate_size(args.tolerance, args.cutoffdistance, args.boxlength)
    print(f"The size is: {size}")

if __name__ == "__main__":
    main()

