#!/usr/bin/env python

# standard import
import csv
import sys
import argparse

def main(args: None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser("h2g2")
    subparsers = parser.add_subparsers(dest='subparser_name', title='subcommands', help='additional help')
    
    subparsers_prepare = subparsers.add_parser("prepare", help='Prepare pipeline run of h2g2')
    subparsers_prepare.add_argument("-r", "--reference", help="path to the reference file you want use")
    subparsers_prepare.add_argument("-a", "--assemblies", nargs="*", required=True, help="list of assemblies you want use")
    
    subparsers_split = subparsers.add_parser("split", help='Generate a new vcf with overlapping variant combine')
    subparsers_split.add_argument("-v", "--vcf-files", nargs="*", required=True, help="Path to vcf")
    subparsers_split.add_argument("-r", "--references", required=True, help="Path to references")
    subparsers_split.add_argument("-o", "--output", required=True, help="Path where output was write")

    subparsers_extract = subparsers.add_parser("extract", help="Create a set of sequence of overlapping variant region")
    subparsers_extract.add_argument("-v", "--vcf-files", nargs="*", required=True, help="Path to vcf")
    subparsers_extract.add_argument("-r", "--references", required=True, help="Path to references")
    subparsers_extract.add_argument("-k", "--kmer-size", type=int, required=True, help="Size of k")
    subparsers_extract.add_argument("-o", "--output", required=True, help="Prefix where output was write")
    
    args = parser.parse_args()
    
    if args.subparser_name == "prepare":
        from .prepare import prepare # import only when it's required
        return prepare(args)
    elif args.subparser_name == "split":
        from .split import split # import only when it's required
        return split(args)
    elif args.subparser_name == "extract":
        from .extract import extract # import only when it's required
        return extract(args)
    else:
        print("You need to use a subcommands", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
