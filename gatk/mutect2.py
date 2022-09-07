#!/usr/bin/env python

## java -jar /home/exacloud/clinical/installedTest/GATK/GenomeAnalysisTK-nightly-2016-03-06.jar -R /opt/installed/galaxy_genomes/hg19/Homo_sapiens_assembly19.fasta -T MuTect2 -stand_call_conf 10 -stand_emit_conf 10 -mbq 0 -PON /home/exacloud/clinical/Galaxy/database/files/034/dataset_34785.dat -I:tumor test_4253.bam -o 4253_10_all.vcf -L solidtumor.interval_list

VERSION = '0.5.7'

import argparse
import subprocess
import shutil
import pysam
import logging
import os

def build_cmd(args):
    """
    Put together the command that will run MuTect2.
    """

    cmd = ["/home/groups/clinical/installedTest/jdk1.8.0_65/bin/java", "-jar", args.mutect2, "-R", args.ref_genome, "-T", "MuTect2", "-I:tumor", "input.bam", "-o", args.output_vcf, "--disable_auto_index_creation_and_locking_when_reading_rods"]

    ### Arguments specifying input resource files.
    if args.intervals:
        cmd.extend(["-L", "intervals.interval_list"])
    if args.interval_str:
        cmd.extend(["-L", args.interval_str])
    if args.panel_of_normals:
        cmd.extend(["-PON", args.panel_of_normals])
    if args.dbsnp:
        cmd.extend(["--dbsnp", args.dbsnp])
    if args.cosmic:
        cmd.extend(["--cosmic", args.cosmic])

    ### Arguments specifying integer parameter values.

    if args.interval_padding:
        cmd.extend(["--interval_padding", args.interval_padding])
    # if args.stand_call_conf:
    #     cmd.extend(["-stand_call_conf", args.stand_call_conf])
    # if args.stand_emit_conf:
    #     cmd.extend(["-stand_emit_conf", args.stand_emit_conf])
    if args.min_base_quality_score:
        cmd.extend(["-mbq", args.min_base_quality_score])
    if args.initial_tumor_lod:
        cmd.extend(["--initial_tumor_lod", args.initial_tumor_lod])
    if args.initial_normal_lod:
        cmd.extend(["--initial_normal_lod", args.initial_normal_lod])
    if args.tumor_lod:
        cmd.extend(["--tumor_lod", args.tumor_lod])
    if args.normal_lod:
        cmd.extend(["--normal_lod", args.normal_lod])
    if args.dbsnp_normal_lod:
        cmd.extend(["--dbsnp_normal_lod", args.dbsnp_normal_lod])
    if args.active_region_in:
        cmd.extend(["--activeRegionIn", "intervals.interval_list"])
    ### Arguments not unique to MuTect2.

    if args.min_pruning:
        cmd.extend(["--minPruning", args.min_pruning])
    if args.min_dangling_branch_length:
        cmd.extend(["--minDanglingBranchLength", args.min_dangling_branch_length])

    ### JHL Added 12/03/16
    if args.kmerSize:
        cmd.extend(["--kmerSize", args.kmerSize])
    if args.downsampling_type != 'NO_SELECT':
        cmd.extend(["--downsampling_type", args.downsampling_type])
    if args.max_alternate_alleles:
        cmd.extend(["--max_alternate_alleles", args.max_alternate_alleles])
    if args.dontUseSoftClippedBases:
        cmd.extend(["--dontUseSoftClippedBases"])

    print(' '.join(cmd))
    return cmd


def cmd_caller(cmd):
    """
    Runs the command that will invoke this tool.
    """

    logging.info("RUNNING: %s" % (' '.join(cmd)))
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if len(stderr):
        print(stderr)
    return p.returncode


def main():

    parser = argparse.ArgumentParser(description='Run MuTect2.')    

    ### Required
    parser.add_argument('mutect2', help='Path to MuTect2.')
    parser.add_argument('ref_genome', help='Path to the reference genome.')
    parser.add_argument('input_bam', help='Input BAM to process.')
    parser.add_argument('output_vcf', help='Output VCF.')
    parser.add_argument('log_file', help='Logging goes here.')

    ### Arguments specifying input resource files.
    parser.add_argument('--intervals', help='File of genomic coordinate intervals to call variants against of form <chrom:start-stop>.')
    parser.add_argument('--interval_str', help='Genomic coordinate intervals to call variants against of form <chrom:start-stop>.')
    parser.add_argument('--panel_of_normals', help='Panel of normals to compare tumor calls against.')
    parser.add_argument('--dbsnp', help='dbSNP file')
    parser.add_argument('--cosmic', help='VCF file of COSMIC sites')

    ### Arguments specifying integer parameter values.
    parser.add_argument('--interval_padding', help='Amount of padding in bp to add to each interval')
    # parser.add_argument('--stand_call_conf', help='The minimum phred-scaled confidence threshold at which variants should be called')
    # parser.add_argument('--stand_emit_conf', help='The minimum phred-scaled confidence threshold at which variants should be emitted')
    parser.add_argument('--min_base_quality_score', help='Minimum base quality required to consider a base for calling')
    parser.add_argument('--initial_tumor_lod', help='Initial LOD threshold for calling tumor variant')
    parser.add_argument('--initial_normal_lod', help='Initial LOD threshold for calling normal variant')
    parser.add_argument('--tumor_lod', help='LOD threshold for calling tumor variant')
    parser.add_argument('--normal_lod', help='LOD threshold for calling normal non-germline')
    parser.add_argument('--dbsnp_normal_lod', help='LOD threshold for calling normal non-variant at dbsnp sites')
    parser.add_argument('--active_region_in', help='Use this interval list file as the active regions to process')
    ### Arguments not unique to MuTect2.
    parser.add_argument('--min_pruning', help='Minimum support to not prune paths in the graph')
    parser.add_argument('--min_dangling_branch_length', help='Minimum length of a dangling branch to attempt recovery')
    ## JHL Added 12/03/16
    parser.add_argument('--kmerSize', help='Kmer size to use in the read threading assembler')
    parser.add_argument('--downsampling_type', help='Type of read downsampling to employ at a given locus [NONE, ALL_READS, BY_SAMPLE]')
    parser.add_argument('--max_alternate_alleles', help='Maximum number of alternate alleles to genotype')
    parser.add_argument('--dontUseSoftClippedBases', action='store_true', help='If specified, we will not analyze soft clipped bases in the reads')

    args = parser.parse_args()

    logging.basicConfig(filename=args.log_file,level=logging.INFO)
    work_bam = "input.bam"
    work_bai = "input.bam.bai"

    logging.info("Symlink " + args.input_bam + " with " + work_bam)

    if os.path.isfile(work_bam):
        os.remove(work_bam)
    if os.path.isfile(work_bai):
        os.remove(work_bai)

    os.symlink(args.input_bam, work_bam)
    logging.info("Indexing " + work_bam)
    pysam.index(work_bam)

    if args.intervals:
        logging.info("Including interval list in working directory.")
        shutil.copyfile(args.intervals, "intervals.interval_list")

    logging.info("Building command in preparation for invocation.")
    cmd = build_cmd(args)
    proc = cmd_caller(cmd)


if __name__ == "__main__":
    main()

