# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import datetime

# using class
from vprimer.param import Param
from vprimer.ini_file import IniFile
from vprimer.conf import Conf
from vprimer.ref_fasta import RefFasta
from vprimer.outlist import OutList

# ---------------------------------------- #
pick_mode_list = ['all', 'indel', 'snp']

ini_section = "[vprimer]"
ini_version = "20210201"

##########################################
# for region
GENOME_ALL = 'all'

##########################################
# allele status
AL_HOMO     = 'homo'
AL_HETERO   = 'hetero'

##########################################
# pick_mode
MODE_ALL    = 'all'
MODE_INDEL  = 'indel'   # indel
MODE_SNPMNV = 'snpmnv'  # snp nmv mind ???
MODE_SNP    = 'snp'     # snp

##########################################
# var_type (variation_type)
OutOfRange  = "oor"
# Code for indel allele, includes substitutions of unequal length
INDEL       = "indel"
# Code for single nucleotide variant allele
SNP         = "snp"
# Code for a multi nucleotide variant allele
MNV         = "mnv"
# mini_indel 1 =< diff_len min_indel_len
MIND        = "mind"


# long_smpl_side
SAME_LENGTH = -1

# analysis SKIP as variant
SKIP_DONT_SKIP      = -1
SKIP_SAME_HOMO      = 1
SKIP_SAME_HETERO    = 2
SKIP_None           = 3
SKIP_DIFF_INGROUP   = 10

# formtxt COMMENT
COMMENT_nop = '-'
COMMENT_dup = 'dup'

# segregation_pattern
segr_ptn_NOP                        = 'nop'
# ./.
segr_ptn_NOT_EXIST_ALLELE           = 'not_exist_allele'
# AA,AA BB,BB CC,CC DD,DD
segr_ptn_SAME_HOMO                  = 'same_homo'
# AB,AB AC,AC AD,AD BC,BC BD,BD CD,CD
segr_ptn_SAME_HETERO                = 'same_hetero'
# AA,BB AA,CC AA,DD BB,CC BB,DD CC,DD
# BB,AA CC,AA DD,AA CC,BB DD,BB DD,CC
segr_ptn_HOMO_HOMO                  = 'hoho_1'
# AA,AB AA,AC AA,AD BB,AB BB,BC BB,BD CC,AC CC,BC CC,CD DD,AD
# DD,BD DD,CD
# AB,AA AC,AA AD,AA AB,BB BC,BB BD,BB AC,CC BC,CC CD,CC AD,DD
# BD,DD CD,DD
segr_ptn_HOMO_HETERO_SHARE          = 'hohe_s1'
# AA,BC AA,BD AA,CD BB,AC BB,AD BB,CD CC,AB CC,AD CC,BD DD,AB
# DD,AC DD,BC
# BC,AA BD,AA CD,AA AC,BB AD,BB CD,BB AB,CC AD,CC BD,CC AB,DD
# AC,DD BC,DD
segr_ptn_HOMO_HETERO_NOT_SHARE      = 'hohe_n2'
# AB,AC AB,AD AB,BC AB,BD AC,AD AC,BC AC,CD AD,BD AD,CD BC,BD
# BC,CD BD,CD
segr_ptn_HETERO_HETERO_SHARE        = 'hehe_s3'
# AB,CD AC,BD AD,BC
segr_ptn_HETERO_HETERO_NOT_SHARE    = 'hehe_n4'

# 使用する制限酵素の最長のrecognition site length
AROUND_SEQ_LEN = 20


def init(prog_name):

    global program_name
    program_name = prog_name

    global now_epochtime, now_datetime_str, now_datetime_form
    now_epochtime, now_datetime_str, now_datetime_form = get_now_time()

    global cwd
    cwd = os.getcwd()

    # local variable
    param = Param()
    ini = IniFile()

    global conf
    conf = Conf()

    global ref
    ref = RefFasta()

    global outlist
    outlist = OutList()

    #--------------------------------
    # get command line parameter
    param = param.get_args()

    # if specified, read ini file
    if param.p['ini_file'] is not None:
        ini.read_ini_file(param.p['ini_file'])

    # collect variables from param, ini, default
    conf.collect_param_ini(param, ini)

    # start logging
    conf.out_dir_logging_start()

    param.open_log()
    conf.open_log_disting()
    conf.open_log_vcffile()
    conf.open_log_currset()
    conf.open_log_enzyme()
    ref.open_log()
    outlist.open_log()


def get_now_time():

    now_datetime = datetime.datetime.now()      # 2021-02-09 10:56:57.209040
    now_epochtime = now_datetime.timestamp()    # 1612835817.20904

    now_dt_ms = now_datetime.strftime("%f")     # 209040
    now_dt_ms2 = now_dt_ms[0:2]                 # 20
    now_dt_hr = now_datetime.strftime("%Y%m%d%H%M%S")   # 20210209102610
    now_datetime_form = now_datetime.strftime("%Y-%m-%d %H:%M:%S")

    now_datetime_str = "{}{}".format(
        now_dt_hr, now_dt_ms2)                          # 2021020910261020

    return now_epochtime, now_datetime_str, now_datetime_form

