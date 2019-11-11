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
    subparsers_prepare.add_argument("-a", "--assembly", nargs="*", required=True, help="list of assembly you want use")
    
    subparsers_split_variant = subparsers.add_parser("split", help='Extract complicate region from vcf')
    subparsers_split_variant.add_argument("-v", "--vcf-files", nargs="*", required=True, help="Path to vcf")
    subparsers_split_variant.add_argument("-r", "--references", required=True, help="Path to references")
    
    subparsers_split_variant.add_argument("-o", "--output", required=True, help="Path where output was write")
    args = parser.parse_args()

    
    if args.subparser_name == "prepare":
        from .prepare import prepare # import only when it's required
        return prepare(args)
    elif args.subparser_name == "split":
        from .split import split # import only when it's required
        return split(args)
    else:
        print("You need to use a subcommands", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
