import random
import sys
import argparse
import csv
import string


def integer_csv_header(snps, samples, delimiter):
    random.seed(42)
    generators = []
    char_set = (string.ascii_letters + string.digits + '"' + "'" + "#&* \t")

    for i in range(0,snps):
        generators.append(lambda: "X")    

    generators.append(lambda: "Y")

    writer = csv.writer(sys.stdout, delimiter=delimiter)
    writer.writerow([g() for g in generators])


def integer_csv(snps, samples, delimiter, target, randSeed):
    random.seed(randSeed)
    generators = []
    char_set = (string.ascii_letters + string.digits + '"' + "'" + "#&* \t")

    for i in range(0,snps):
        generators.append(lambda: random.randint(0, 2))    
    
    generators.append(lambda: target)

    writer = csv.writer(sys.stdout, delimiter=delimiter)
    for x in xrange(samples):
        writer.writerow([g() for g in generators])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates CSV file representing a case-control dataset.')
    parser.add_argument('snps', type=int, help='number of SNP columns to generate')
    parser.add_argument('samples', type=int, help='number of samples to generate (half cases/controls)')
    parser.add_argument('--delimiter', type=str, default=',', required=False, help='the CSV delimiter')

    args = parser.parse_args()
    integer_csv_header(args.snps, args.samples, args.delimiter)  # The header 
    integer_csv(args.snps, (args.samples / 2), args.delimiter, 1, 19)  # CASES
    integer_csv(args.snps, (args.samples / 2), args.delimiter, 0, 13)  # CONTROLS


