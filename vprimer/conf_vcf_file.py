# -*- coding: utf-8 -*-

import sys
import os
import errno
import re

import pprint
import vcfpy
import pandas as pd

# global variants
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf

class ConfVcfFile(object):

    def open_log_vcffile(self):

        # in glv
        global log
        log = LogConf.open_log(__name__)


    def prepare_vcf(self):
        """
        """

        # self.user_vcf_file_path: user asigned vcf in param or ini file
        if os.path.isfile(self.user_vcf_file_path):
            log.info("found {}.".format(self.user_vcf_file_path))
        else:
            log.info("not found {}. exit.".format(
                self.user_vcf_file_path))
            sys.exit(1)

        # symbolic link user_vcf_file_path to vcf_file
        if os.path.isfile(self.vcf_file_slink_system):
            log.info("found {}.".format(self.vcf_file_slink_system))
        else:
            utl.ln_s(
                self.user_vcf_file_path, self.vcf_file_slink_system)

        # vcf_file
        if os.path.isfile(self.vcf_file_path):
            log.info("found {}.".format(self.vcf_file_path))
        else:
            cmd1 = "{} annotate --threads {} -O z -x ^FMT/GT -o {} {}".format(
                'bcftools',
                self.parallel_full_thread,
                self.vcf_file_path,
                self.vcf_file_slink_system)

            # if bcftools 1.10.2, it have --force option.
            # if 1.9.0, if we have an error
            # No matching tag in -x ^FMT/GT
            # we will not worry, continue

            err_str = utl.try_exec_error(cmd1)

            if err_str == '':
                pass

            elif err_str.startswith('No matching tag in'):

                log.error("we will go ahead if bcftools says, >{}<.".format(
                    err_str))

                utl.rm_f(self.vcf_file)
                utl.ln_s(
                    self.vcf_file_slink_system, self.vcf_file_path)

            else:
                log.error(">{}<".format(err_str))
                sys.exit(1)

            # make tbi
            utl.tabix(self.vcf_file_path)


    def save_vcf_sample_name_txt(self):
        '''
        '''

        # exist or not, vcf_sample_name_file
        if os.path.isfile(self.vcf_sample_name_file):
            log.info("found.\n{}".format(self.vcf_sample_name_file))
            # Make a backup of vcf_sample_name_file
            # as it may have been edited by the user
            utl.save_to_tmpfile(self.vcf_sample_name_file, True, True)

        else:
            sample_name_list = list()
            # if not, read vcf and pick sample_name
            log.info("not found {}.".format(self.vcf_sample_name_file))
            sample_name_list += ["#{}\t{}\t{}\t{}\t{}".format(
                'no', 'group', 'nickname', 'basename', 'fullname')]

            sample_name_list += self._pick_vcf_sample_list(
                self.vcf_file_path)

            # backup
            utl.save_to_tmpfile(self.vcf_sample_name_file)

            # write to vcf_sample_name_file
            with open(self.vcf_sample_name_file, mode='w') as f:
                f.write("{}\n".format("\n".join(sample_name_list)))

            log.info("save.\n{}".format(self.vcf_sample_name_file))


    def make_vcf_sample_variable(self):

        # a simple list. Used for existence confirmation, etc.
        vcf_sample_nickname_list = list()
        vcf_sample_basename_list = list()
        vcf_sample_fullname_list = list()

        # _is_sample_name(sample)
        # _get_nickname(sample)
        # _get_basename(sample)
        # _get_fullname(sample)

        # The key is nickname, which gives you the basename and fullname.
        vcf_sample_nickname_dict = dict()
        # The key is basename, which gives you the nickname and fullname.
        vcf_sample_basename_dict = dict()
        # The key is fullname, which gives you the nickname and basename.
        vcf_sample_fullname_dict = dict()

        group_members_vcf_str = ""

        group_dict = dict()
        p_list = list()

        with open(self.vcf_sample_name_file, mode='r') as r:
            # 'no', 'group', 'nickname', 'basename', 'fullname')]
            for liner in r:
                r_line = liner.strip()

                # for print to STDOUT
                if r_line == '':
                    continue

                # print except ''
                p_list.append(r_line)
                #print("{}".format(r_line), file=sys.stdout)

                # comment line
                if r_line.startswith('#'):
                    continue

                # group may be separated by comma
                no, groups, nickname, basename, fullname = \
                    r_line.split('\t')

                # group_dict, group_members_vcf_str
                sep_groups = groups.split(',')
                for sep_group in sep_groups:
                    if not sep_group in group_dict:
                        group_dict[sep_group] = list()
                    group_dict[sep_group].append(nickname)

                sample_dict = dict()
                sample_dict = {
                    'no': no,
                    'nickname': nickname,
                    'basename': basename,
                    'fullname':fullname
                }

                # if nickname duplicated, it's ok if contents is same
                if nickname in vcf_sample_nickname_dict:
                    if basename == \
                        vcf_sample_nickname_dict[nickname]['basename'] and \
                        fullname == \
                        vcf_sample_nickname_dict[nickname]['fullname']:
                        # It's ok
                        pass
                    else:
                        log.info("{} >{}< {}".format(
                            "same nickname", nickname,
                            "have different contents. exit."))
                        sys.exit(1)

                vcf_sample_nickname_list.append(nickname)
                vcf_sample_basename_list.append(basename)
                vcf_sample_fullname_list.append(fullname)

                vcf_sample_nickname_dict[nickname] = sample_dict
                vcf_sample_basename_dict[basename] = sample_dict


        log.info("\n{}\n{}\n".format(
            self.vcf_sample_name_file, "\n".join(p_list)))

        for group in sorted(group_dict):
            group_members_vcf_str += "{}:{},".format(
                group, ",".join(group_dict[group]))

        group_members_vcf_str = re.sub(r",$", "", group_members_vcf_str)

        return vcf_sample_nickname_list, \
            vcf_sample_basename_list, \
            vcf_sample_fullname_list, \
            vcf_sample_nickname_dict, \
            vcf_sample_basename_dict, \
            vcf_sample_fullname_dict, \
            group_members_vcf_str


    def _pick_vcf_sample_list(self, vcf_file_path):
        ''' open vcf file and pick sample information as list
        '''

        sample_list = list()
        reader = vcfpy.Reader.from_path(vcf_file_path)

        for (sample_no, vcf_sample_name) in enumerate(
            reader.header.samples.names, 1):

            sample_basename_bam = os.path.basename(vcf_sample_name)
            sample_basename = re.sub(r"\.bam$", "", sample_basename_bam)
            group_name = "-"

            sample_list.append("{}\t{}\t{}\t{}\t{}".format(
                sample_no, group_name, sample_basename,
                sample_basename_bam, vcf_sample_name))

        return sample_list


    def get_vcf_pos_info(self):

        vcf_pos_info_list = list()
        chrom_list = list()
        chrom_dict = dict()
        dict_cnt = 0

        reader = vcfpy.Reader.from_path(self.vcf_file_path)
        for record in reader:
            chrom = record.CHROM
            pos = record.POS

            if chrom not in chrom_list:
                #print(chrom)
                chrom_list.append(chrom)

                if dict_cnt != 0:
                    vcf_pos_info_list.append(chrom_dict)
                dict_cnt += 1
                #pprint.pprint(vcf_pos_info_list)

                chrom_dict = dict()

                chrom_dict = {
                    'chrom': chrom,
                    'start': pos,
                    'end': pos,
                }

            chrom_dict['end'] = pos    
            
        vcf_pos_info_list.append(chrom_dict)

        #pprint.pprint(chrom_list)
        #pprint.pprint(vcf_pos_info_list)

        pd_json = pd.json_normalize(vcf_pos_info_list)
        log.info("\n{}\n".format(pd_json))





