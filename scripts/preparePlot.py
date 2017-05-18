#!/usr/bin/env python

import sys
import os
import argparse
import logging

def getWhInfo(wh_log_file):
    info = {}
    with open(wh_log_file, 'r') as whlog:
        for line in whlog.readlines():
            if line.startswith("MEC cost:"):
                info["score"] = int(line.strip().split(':')[1])
            if line.startswith("Time spent phasing:"):
                info["time"] = float(line.strip().replace(" ", "").split(':')[1][:-1])
        return info

def getScore(log_file):
    with open(log_file, 'r') as log:
        for line in log.readlines():
            line = line.strip()
            if line.startswith("* ") and " | " in line:
                cleanline = line.split(" | ")[1]
                if cleanline.startswith("OPTIMUM:"):
                    return int(cleanline.split(":")[1])

def getSER(diff_file):
    with open(diff_file, 'r') as d:
        nl = 0
        for l in d.readlines():
            nl += 1
            if nl == 23:
                return(l.rstrip().replace(" ", "").split(":")[1]).replace("%", "")

def getTimeMem(log_file):
    info = {}
    with open(log_file, 'r') as lf:
        for line in lf.readlines():
            line = line.strip()
            if line.startswith("Elapsed (wall clock) time"):
                wallv = line.split(": ")[1].split(":")
                hours = 0 if len(wallv) == 2 else int(wallv[-3])
                mins = int(wallv[-2])
                secs = float(wallv[-1])
                info["time"] = ((hours*60) + mins)*60+secs
            if line.startswith("Maximum resident set size"):
                info["mem"] = int(line.split(":")[1])/(1024)#*1024)
        return info

def main():
    parser = argparse.ArgumentParser(prog = "preparePlot",
                                     description = "Prepare Hplotyping data.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inck', help = "Increase k HapChat dir.",
                        required = True, dest = 'inck_dir')
    parser.add_argument('-p', '--pre-input', help = "Input WIF file dir.",
                        required = True, dest = 'wif_dir')
    parser.add_argument('-w', '--whatshap', help = "WhatsHap dir.",
                        required = False, dest = 'wh_dir')
    parser.add_argument('-o', '--out-file', help = "CSV output file.",
                        required = True, dest = 'out_file')
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
    with open(args.out_file, 'w') as out:
        out.write("Tool,Data,Technology,Individual,Chr,MeanCov,RawRealigned,WhDowns")
        out.write(",Merge,MergeE,MergeM,MergeT,MergeN,MergeTimeSec,MergeMaxMemMB")
        out.write(",RndDowns,RndDownsSeed,RndDownsMaxCov,RndDownsTimeSec,RndDownsMaxMemMB")
        out.write(",Epsilon,Alpha")
        out.write(",SwErrRatePerc,MecScore,PhasTimeSec,PhasMaxMemMB")
        out.write("\n")
        # Hapchat
        num_inck = 0
        logging.info("Parsing Inck files")
        for df in os.listdir(args.inck_dir):
            if df.endswith(".diff"):
                num_inck += 1
                out.write("HapChat,")
                ds = df.rstrip().split(".")[:-1]
                dataset = ".".join(ds)
                ser = getSER(args.inck_dir + df)
                tfile = args.inck_dir + dataset + ".time"
                if not os.path.isfile(tfile):
                    logging.error("File not found: %s", tfile)
                    exit()
                info = getTimeMem(tfile)
                lfile = args.inck_dir + dataset + ".log"
                if not os.path.isfile(lfile):
                    logging.error("File not found: %s", lfile)
                    exit()
                score = getScore(lfile)
                # print(dataset + " " + str(score) +
                #       ser + " " +
                #       str(info["time"]) + " " +
                #       str(info["mem"]))
                merge = "no,NA,NA,NA,NA,NA,NA"
                if(ds[7] != "no_merged"):
                    mfile = args.wif_dir + ".".join(ds[0:8]) + ".wif.time"
                    if not os.path.isfile(mfile):
                        logging.error("File not found: %s", mfile)
                        exit()
                    merge_info = getTimeMem(mfile)
                    mfields = ds[7].split('_')
                    merge = "yes" + "," + \
                            "0." + mfields[1][1:] + "," + \
                            "0." + mfields[2][1:] + "," + \
                            "1e+" + mfields[3][1:] + "," + \
                            "1e+" + mfields[4][1:] + "," + \
                            str(merge_info["time"]) + "," + \
                            str(merge_info["mem"])
                downs = "no,NA,NA,NA,NA"
                if(ds[8] != "no_downs"):
                    dfile = args.wif_dir + ".".join(ds[0:7])
                    if(ds[7] != "no_merged"):
                        dfile += "." + ds[7]
                    dfile += ".wif.sample_" + "_".join(ds[8].split('_')[1:]) + ".time"
                    if not os.path.isfile(dfile):
                        logging.error("File not found: %s", dfile)
                        exit()
                    downs_info = getTimeMem(dfile)
                    dfields = ds[8].split('_')
                    downs = "yes" + "," + \
                            dfields[1][1:] + "," + \
                            dfields[2][1:] + "," + \
                            str(downs_info["time"]) + "," + \
                            str(downs_info["mem"])
                out.write(",".join(ds[0:4]) + "," +
                          ds[4].replace("cov", "") + "," +
                          ds[5] + "," +
                          ds[6].replace("h", "") + "," +
                          merge + "," +
                          downs + "," +
                          "0." + ds[9].split("_")[0] + "," +
                          "0." + ds[9].split("_")[1] + "," +
                          ser + "," +
                          str(score) + "," +
                          str(info["time"]) + "," +
                          str(info["mem"]))
                out.write("\n")
        logging.info("Parsed %d Inck files.", num_inck)
        # WhatsHap
        num_wh = 0
        logging.info("Parsing WhatsHap files")
        for df in os.listdir(args.wh_dir):
            if df.endswith(".diff"):
                num_wh += 1
                out.write("WhatsHap,")
                ds = df.rstrip().split(".")[:-1]
                dataset = ".".join(ds)
                ser = getSER(args.wh_dir + df)
                tfile = args.wh_dir + dataset + ".time"
                if not os.path.isfile(tfile):
                    logging.error("File not found: %s", tfile)
                    exit()
                info = getTimeMem(tfile)
                whlfile = args.wh_dir + dataset + ".log"
                if not os.path.isfile(whlfile):
                    logging.error("File not found: %s", whlfile)
                    exit()
                wh_info = getWhInfo(whlfile)
                out.write(",".join(ds[0:4]) + "," +
                          ds[4].replace("cov", "") + "," +
                          ds[5] + "," +
                          ds[6].replace("h", "") + "," +
                          "no,NA,NA,NA,NA,NA,NA" + "," +
                          "no,NA,NA,NA,NA" + "," +
                          "NA,NA" + "," +
                          ser + "," +
                          str(wh_info["score"]) + "," +
                          str(wh_info["time"]) + "," +
                          str(info["mem"]))
                out.write("\n")
                
        logging.info("Parsed %d WhatsHap files.", num_wh)
        
    logging.info("Program Finshed")

if __name__ == "__main__":
    main()
