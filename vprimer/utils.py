# -*- coding: utf-8 -*-

import sys
import os
import re
import glob
import copy
import pprint

#
import pandas as pd
#import time
import datetime

import logging
import logging.config

# global configuration
import vprimer.glv as glv

import subprocess as sbp
from subprocess import PIPE

#=========================================================
def open_log():
    '''
    '''
    logging.config.dictConfig(glv.conf.log.config)
    global log
    log = logging.getLogger(__name__)
    log.info("logging start {}".format(__name__))


def elapsed_time(now_time, start_time):
    '''
    '''
    elapsed_time = now_time - start_time
    #3:46:11.931354
    return "elapsed_time {}".format(
        datetime.timedelta(seconds=elapsed_time))


def strip_hash_comment(line):
    '''
    '''
    return line.split('#')[0].strip()


def full_path(file_name):
    '''
    after checking either file_name is absolute one or relative one,
    make an appropriate absolute path
    '''

    full_path = ''

    if file_name == '' or file_name is None:
        pass

    else:
        if file_name.startswith('/'):
            full_path = file_name
        elif file_name == '':
            pass
        else:
            full_path = "{}/{}".format(glv.cwd, file_name)

    return full_path


def prelog(line, name):
    '''
    '''
    print("prelog: {} {}".format(name, line), file=sys.stderr)


def try_exec(cmd, can_log=True):
    '''
    '''
    try:
        if can_log == True:
            log.info("do {}".format(cmd))
        else:
            #print("do {}".format(cmd), file=sys.stderr)
            prelog("do {}".format(cmd), __name__)

        sbp.run(cmd,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
            shell=True,
            check=True)

    except sbp.CalledProcessError as e:
        if can_log == True:
            log.error("{}.".format(e.stderr))
        else:
            prelog("{}.".format(e.stderr), cmd)

        sys.exit(1)


def try_exec_error(cmd):
    '''
    '''
    # https://qiita.com/HidKamiya/items/e192a55371a2961ca8a4

    err_str = ''

    try:
        log.info("do {}".format(cmd))

        sbp.run(cmd,
            stdout=PIPE,
            stderr=PIPE,
            text=True,
            shell=True,
            check=True)

    except sbp.CalledProcessError as e:
        err_str = e.stderr

    return err_str


def save_to_tmpfile(file_path, can_log=True, copy_mod=False):
    '''
    '''

    ret = False
    if os.path.isfile(file_path):
        # /a/b/c.txt
        # /a/b/bak/c.txt
        dirname_file = os.path.dirname(file_path)
        basename_file = os.path.basename(file_path)

        file_bak_path = "{}/{}".format(
            glv.conf.out_bak_dir_path, basename_file)

        ret = True
        new_file_path = "{}.{}.bak".format(
            file_bak_path, glv.now_datetime_str)

        mode = "mv"
        if copy_mod == True:
            mode = "cp"
            cp_f(file_path, new_file_path)

        else:
            mv_f(file_path, new_file_path)
        
        if can_log:
            log.info("found {}. {} to {}".format(
                file_path, mode, new_file_path))

    return ret


def ln_s(source_file, symbolic_link_file):
    '''
    '''
    cmd1 = "{} {} {}".format(
        'ln -s', source_file, symbolic_link_file)
    try_exec(cmd1)


def rm_f(remove_file):
    '''
    '''
    cmd1 = "{} {}".format(
        'rm -f', remove_file)
    try_exec(cmd1)


def makedirs(dir_name):
    '''
    '''
    cmd1 = "{} {}".format(
        'mkdir -p', dir_name)
    try_exec(cmd1, False)


def mv_f(source_file, renamed_file):
    '''
    '''
    cmd1 = "{} {} {}".format(
        'mv', source_file, renamed_file)
    try_exec(cmd1, False)


def cp_f(source_file, copyed_file):
    '''
    '''
    cmd1 = "{} {} {}".format(
        'cp', source_file, copyed_file)
    try_exec(cmd1, False)


def tabix(vcf_file):
    '''
    '''
    #-f -p vcf
    cmd1 = "{} {} {}".format(
        'tabix',
        '-f -p vcf',
        vcf_file)

    try_exec(cmd1)


#=========================================================
def write_df_to_csv(dataframe, out_file, index=True, force=True):
    '''
    '''
    file_name = out_file
    moved = False
    if force == True:
        moved = save_to_tmpfile(out_file)

    dataframe.to_csv(out_file, sep='\t', index=index)


#def progress_check(now_progress):
#    '''
#    '''
#    stat = False    # False if don't do this progress
#    param_progress = glv.conf.progress

#    log.info("now_progress={} param_progress={}".format(
#        now_progress, param_progress))

    #log.debug("now_progress={} {} param_progress={} {}".format(
    #    now_progress,
    #    now_progress_no,
    #    param_progress,
    #    param_progress_no))

#    if param_progress == 'all':
#        stat = True

#    else:
#        now_progress_no = glv.outlist.outf_prefix[now_progress]['no']
#        param_progress_no = glv.outlist.outf_prefix[param_progress]['no']
#        if now_progress_no >= param_progress_no:
#            stat = True

#    return stat

def print_distin_info(proc_name, distin_dict, proc_cnt, simple=False):

    reg = distin_dict['region']
    reg_dict = glv.conf.regions_dict[reg]

    log.info("proc:        {}, {} / {}".format(
        proc_name, proc_cnt, len(glv.outlist.distin_files)))

    if simple == True:
        log.info("disting_str  {}\n".format(distin_dict['disting_str']))
        return
    else:
        log.info("disting_str  {}".format(distin_dict['disting_str']))

    log.info("all regions  {}".format(
        ', '.join(distin_dict['regions_list'])))
    log.info("groups       {} / {}".format(
        distin_dict[0], distin_dict[1]))
    log.info("region_name  {}".format(distin_dict['region']))
    log.info("region_str   {}".format(reg_dict['reg']))
    log.info("pick_mode    {}".format(distin_dict['pick_mode']))
    log.info("indel_size   {}".format(distin_dict['indel_size']))
    log.info("members      {} : {}".format(
        distin_dict[0],
        ', '.join(glv.conf.group_members_dict[distin_dict[0]])))
    log.info("members      {} : {}\n".format(
        distin_dict[1],
        ', '.join(glv.conf.group_members_dict[distin_dict[1]])))


def decide_action_stop(now_progress):
    '''
    '''

    # Current location number
    now_progress_no = glv.outlist.outf_prefix[now_progress]['no']

#    print("now_progress={}, now_progress_no={}".format(
#        now_progress, now_progress_no))

    # User-specified start point and number
    progress_name = glv.conf.progress
    progress_no = 100

    # User-specified stop point and number
    stop_name = glv.conf.stop
    stop_no = 100

    # If progress_name is "all", progress_no is 0
    if progress_name != "all":
        progress_no = glv.outlist.outf_prefix[progress_name]['no']

    # If stop_name is "no", stop_no is 100
    if stop_name != "no":
        stop_no = glv.outlist.outf_prefix[stop_name]['no']

    ret_status = "stop"

#    print("progress_name={}, progress_no={}".format(
#        progress_name, progress_no))
#    print("stop_name={}, stop_no={}".format(
#        stop_name, stop_no))

    # decide
    if progress_no == 100:
        ret_status = "action"
    elif now_progress_no < progress_no:
        ret_status = "gothrough"
    else:
        ret_status = "action"

    if stop_no < now_progress_no:
        ret_status = "stop"

#    print("ret_status={}".format(ret_status))

    return ret_status


#def stop(now_progress):
#    '''
#    '''
#    if glv.conf.stop == 'no':
#        return

#    now_progress_no = glv.outlist.outf_prefix[now_progress]['no']
#    param_stop_no = glv.outlist.outf_prefix[glv.conf.stop]['no']
#    if now_progress_no >= param_stop_no:
#        log.info("stop {}".format(glv.conf.stop))
#        sys.exit(1)


def is_my_pick_mode(var_type, pick_mode):
    '''
    '''

    ret = False

    if var_type == glv.OutOfRange:
        # Exceeds the specified indel length
        pass

    elif pick_mode == glv.MODE_ALL:
        # return True if pick_mode is ALL, no matter what var_type is
        ret = True

    elif pick_mode == glv.MODE_INDEL:
        if var_type == glv.INDEL:
            # return True if pick_mode is INDEL and var_type is INDEL
            ret = True

    elif pick_mode != glv.MODE_INDEL:
        # return True if pick_mode is not INDEL and var_type is not INDEL
        if var_type != glv.INDEL:
            ret = True

    return ret


def check_for_files(filepath):
    '''
    '''
    # filepath is pattern
    fobj_list = list()

    for filepath_object in glob.glob(filepath):

        if os.path.isfile(filepath_object):
            fobj_list.append(filepath_object)

    return sorted(fobj_list)


def is_same_gt(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
    same_gt = False

    # same homo: AA,AA
    if is_same_homo(s0_0, s0_1, s1_0, s1_1):
        same_gt = True

    # same hetero: AB,AB
    elif is_same_hetero(s0_0, s0_1, s1_0, s1_1):
        same_gt = True

    return same_gt


def is_None(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
    return s0_0 is None or s0_1 is None or \
           s1_0 is None or s1_1 is None


def is_homo_homo(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
           #  1 == 2           3 == 4
    return s0_0 == s0_1 and s1_0 == s1_1


def is_homo_hetero(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
           #  1 == 2           3 != 4
    return s0_0 == s0_1 and s1_0 != s1_1 or \
           s0_0 != s0_1 and s1_0 == s1_1
           #  1 != 2           3 == 4


def is_hetero_hetero(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
           #  1 != 2           3 != 4
    return s0_0 != s0_1 and s1_0 != s1_1

def is_same_homo(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
    return is_homo_homo(s0_0, s0_1, s1_0, s1_1) and \
           s0_0 == s1_0
           #  1 == 3

def is_same_hetero(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
    return is_hetero_hetero(s0_0, s0_1, s1_0, s1_1) and \
           s0_0 == s1_0 and s0_1 == s1_1
           #  1 == 3           2 == 4

def is_share(s0_0, s0_1, s1_0, s1_1):
    '''
    '''
           #  1 == 3          1 == 4
    return s0_0 == s1_0 or s0_0 == s1_1 or \
           s0_1 == s1_0 or s0_1 == s1_1
           #  2 == 3          2 == 4


def is_Not_NoneAndNull(value):
    '''
    '''

    ret = False

    is_none = False
    is_null = False

    if value is None:
        is_none = True

    if value == "":
        is_null = True

    if is_none != True and is_null != True:
        ret = True

    return ret


def sort_file(
        proc, distin_dict, out_file_name,
        nm_chrom, nm_pos, nm_order, n):
    '''
    '''

    # sort command option
    if n == 'number':
        n = 'n'
    else:
        n = ''


    # cmd.sort
    hdr_dict = distin_dict[proc]['hdr_dict']
    out_txt_file = distin_dict[proc]['out_path']
    sorted_file = "{}.sorted".format(out_txt_file)

    col_chrom = hdr_dict[nm_chrom]
    col_pos = hdr_dict[nm_pos]
    col_order = hdr_dict[nm_order]

    cmd_sort = "{} -k {},{} -k {},{}n -k {},{}{} {} > {}".format(
        'sort',
        col_chrom, col_chrom,
        col_pos, col_pos,
        col_order, col_order,
        n,
        out_txt_file,
        sorted_file)

    try_exec(cmd_sort)

    # rm file
    rm_f(out_txt_file)

    # make header.txt
    out_txt_header = "{}.header_txt".format(out_txt_file)
    with open(out_txt_header, mode='w') as f:
        f.write("{}\n".format(distin_dict[proc]['hdr_text']))


    # cat header file
    cmd_cat = "{} {} {} > {}".format(
        'cat',
        out_txt_header,
        sorted_file,
        out_txt_file)

    try_exec(cmd_cat)

    # rm header sorted
    rm_f(out_txt_header)
    rm_f(sorted_file)


def get_ordered_sample_list(group_list):
    '''
    '''

    sample_nickname_ordered_list = list()
    sample_fullname_ordered_list = list()

    # glv.conf.vcf_sample_nickname_list
    group0 = group_list[0]
    group1 = group_list[1]

    # Deep copy is required here
    # sample_nickname_ordered_list
    sample_nickname_ordered_list = \
        get_sample_list_from_groupname(group_list, "nickname")
    sample_fullname_ordered_list = \
        get_sample_list_from_groupname(group_list, "fullname")

    nogrp_nickname_list = list()
    nogrp_fullname_list = list()

    # not grouped nickname and fullname
    for nickname in glv.conf.vcf_sample_nickname_list:
        fullname = get_fullname(nickname)

        # The rest of the members behind
        if nickname not in sample_nickname_ordered_list:
            nogrp_nickname_list += [nickname]
            nogrp_fullname_list += [fullname]

    sample_nickname_ordered_list += nogrp_nickname_list
    sample_fullname_ordered_list += nogrp_fullname_list

    return \
        sample_nickname_ordered_list, \
        sample_fullname_ordered_list


def get_sample_list_from_groupname(group_list, need_name):
    '''
    '''

    all_sample_list = list()

    for group_name in group_list:

        if group_name == "ref":
            continue

        sample_list = glv.conf.group_members_dict[group_name]

        for sn in sample_list:

            if need_name == "nickname":
                sample_name = get_nickname(sn)
            elif need_name == "basename":
                sample_name = get_basename(sn)
            elif need_name == "fullname":
                sample_name = get_fullname(sn)
            else:
                sample_name = get_nickname(sn)

            all_sample_list.append(sample_name)

    return all_sample_list


#-----------------------------------------------------
def is_sample_name(sample_name):
    '''
    '''
    ret = False

    if sample_name.lower() == "ref":
        ret = True

    elif sample_name in glv.conf.vcf_sample_nickname_list:
        ret = True

    elif sample_name in glv.conf.vcf_sample_basename_list:
        ret = True

    elif sample_name in glv.conf.vcf_sample_fullname_list:
        ret = True

    return ret


def get_nickname(sample_name):
    '''
    '''
    nickname = ""

    if sample_name.lower() == "ref":
        nickname = "ref"

    elif sample_name in glv.conf.vcf_sample_nickname_list:
        nickname = sample_name

    elif sample_name in glv.conf.vcf_sample_basename_list:
        nickname = glv.conf.vcf_sample_basename_dict[sample_name]['nickname']

    elif sample_name in glv.conf.vcf_sample_fullname_list:
        nickname = glv.conf.vcf_sample_fullname_dict[sample_name]['nickname']

    return nickname


def get_basename(sample_name):
    '''
    '''
    basename = ""

    if sample_name.lower() == "ref":
        basename = "ref"

    elif sample_name in glv.conf.vcf_sample_nickname_list:
        basename = glv.conf.vcf_sample_nickname_dict[sample_name]['basename']

    elif sample_name in glv.conf.vcf_sample_basename_list:
        basename = sample_name

    elif sample_name in glv.conf.vcf_sample_fullname_list:
        basename = glv.conf.vcf_sample_fullname_dict[sample_name]['basename']

    return basename


def get_fullname(sample_name):
    '''
    '''
    fullname = ""

    if sample_name.lower() == "ref":
        fullname = "ref"

    elif sample_name in glv.conf.vcf_sample_nickname_list:
        fullname = glv.conf.vcf_sample_nickname_dict[sample_name]['fullname']

    elif sample_name in glv.conf.vcf_sample_basename_list:
        fullname = glv.conf.vcf_sample_basename_dict[sample_name]['fullname']

    elif sample_name in glv.conf.vcf_sample_fullname_list:
        fullname = sample_name

    return fullname


# primer.py: prepare_from_marker_file,
# formtxt.py: _prepare_from_primer_file
def get_basic_primer_info(df_row, hdr_dict):

    marker_id = ""
    if 'marker_id' in hdr_dict.keys():
        marker_id = str(df_row[hdr_dict['marker_id']])

    chrom = str(df_row[hdr_dict['chrom']])
    pos = int(df_row[hdr_dict['pos']])

    targ_grp = str(df_row[hdr_dict['targ_grp']])
    g0_name, g1_name = targ_grp.split(',')

    targ_ano = str(df_row[hdr_dict['targ_ano']])
    g0_ano, g1_ano = map(int, targ_ano.split(','))

    vseq_gno_str = df_row[hdr_dict['vseq_gno_str']]

    gts_segr_lens = df_row[hdr_dict['gts_segr_lens']]
    var_type = str(df_row[hdr_dict['var_type']])

#    set_enz_cnt = str(df_row[hdr_dict['set_enz_cnt']])
    set_enz_cnt = ""
    if 'set_enz_cnt' in hdr_dict:
        set_enz_cnt = str(df_row[hdr_dict['set_enz_cnt']])
    else:
        set_enz_cnt = str(df_row[hdr_dict['set_n']])

    marker_info = ""
    target_gno = -1
    target_len = 0
    enzyme_name = '-'
    digest_pattern = '-'

    if 'marker_info' in hdr_dict:
        marker_info = str(df_row[hdr_dict['marker_info']])
        if var_type == glv.INDEL:
            longer_group, \
            longer_length, \
            shorter_length, \
            diff_length, \
            digested_pos = \
                map(int, marker_info.split('.'))
            target_gno = longer_group
            target_len = longer_length
        
        else:
            enzyme_name, \
            digested_gno, \
            found_pos, \
            digest_pattern, \
            digested_pos = \
                marker_info.split('.')

            target_gno = int(digested_gno)
            target_len = int(found_pos)

    vseq_lens_ano_str = ""
    if 'vseq_lens_ano_str' in hdr_dict:
        vseq_lens_ano_str = str(df_row[hdr_dict['vseq_lens_ano_str']])

    return \
        marker_id, \
        chrom, \
        pos, \
        targ_grp, \
        g0_name, \
        g1_name, \
        targ_ano, \
        g0_ano, \
        g1_ano, \
        vseq_gno_str, \
        gts_segr_lens, \
        var_type, \
        set_enz_cnt, \
        marker_info, \
        vseq_lens_ano_str, \
        enzyme_name, \
        digest_pattern, \
        target_gno, \
        target_len


def flip_gno(gno):

    flip_gno = 0
    if gno == 0:
        flip_gno = 1

    return flip_gno


