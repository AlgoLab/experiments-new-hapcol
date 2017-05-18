#!/usr/bin/env python

import sys
import os
import argparse
import logging

def main():
    parser = argparse.ArgumentParser(prog = "checkFileNames",
                                     description = "Check File Name Structure.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inck', help = "Increase k HapChat dir.",
                        required = True, dest = 'inck_dir')
    parser.add_argument('-w', '--whatshap', help = "WhatsHap dir.",
                        required = False, dest = 'wh_dir')
    parser.add_argument('-v', '--verbose',
                        help='increase output verbosity',
                        action='count', default=0)
    args = parser.parse_args()

    if args.verbose == 0:
        log_level = logging.INFO
    elif args.verbose == 1:
        log_level = logging.DEBUG
    else:
        log_level = logging.DEBUG
    
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt="%y%m%d %H%M%S")

    logging.info("Program Started")
    logging.info("Parsing Inck files")
    for df in os.listdir(args.inck_dir):
        tokens = df.rstrip().split(".")
        if(len(tokens) == 9):
            nn = ".".join(tokens[:-2]) + ".no_merged.no_downs." + ".".join([tokens[-2], tokens[-1]])
            os.rename(args.inck_dir + df.rstrip(), args.inck_dir + nn)
            print(df.rstrip())
        if(len(tokens) == 10):
            nn = ".".join(tokens[:-3]) + ".no_merged." + ".".join([tokens[-3], tokens[-2], tokens[-1]])
            os.rename(args.inck_dir + df.rstrip(), args.inck_dir + nn)
            print(df.rstrip())
    for df in os.listdir(args.inck_dir):
        tokens = df.rstrip().split(".")
        if(len(tokens) != 11):
            logging.error("Error: %s", df.rstrip())
            exit()
    logging.info("Program Finshed")

if __name__ == "__main__":
    main()

