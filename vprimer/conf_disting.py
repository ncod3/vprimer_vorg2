# -*- coding: utf-8 -*-

import sys
import os
import errno
import re

import pprint

# global variants
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf

class ConfDistinG(object):


    def open_log_disting(self):

        global log
        log = LogConf.open_log(__name__)

    def set_regions_str(self):
        ''' regions 
        Variable origin priority and standardize the string format
            1. command parameter        <target>
            2. command parameter        <regions>
            3. ini file                 <regions>

        Format:
            region_name : chrom_name : start_pos-end_pos
        '''

        regions_str = ""

        # 1. self.conf_dict['target']['easy']
        #log.debug("1. self.conf_dict['target']['param']={}".format(
        #    self.conf_dict['target']['param']))
        #log.debug("2. self.conf_dict['regions']['param']={}".format(
        #    self.conf_dict['regions']['param']))
        #log.debug("3. self.conf_dict['regions']['ini']={}".format(
        #    self.conf_dict['regions']['ini']))

        # "" or None
        if utl.is_Not_NoneAndNull(self.conf_dict['target']['param']):
            regions_str = re.sub(r";", ",",
                self.conf_dict['target']['param'])
            # EASY_MODE will format variables later
            regions_str = "<EASY_MODE>{}".format(regions_str)

        elif utl.is_Not_NoneAndNull(self.conf_dict['regions']['param']):
            regions_str = re.sub(r";", ",",
                self.conf_dict['regions']['param'])

        elif utl.is_Not_NoneAndNull(self.conf_dict['regions']['ini']):
            regions_str = re.sub(r";", ",",
                self.conf_dict['regions']['ini'])

        # ini use '=' as separator, so substitute to ','
        regions_str = re.sub(r"=", ",", regions_str)
        #log.debug("regions_str={}".format(regions_str))

        return regions_str


    def set_group_members_str(self):
        ''' group_members
        Variable origin priority
        1. command parameter        <a_sample,b_sample>
        2. sample_name in_vcf_file  <vcf>
        3. command parameter        <group_members>
        4. ini file                 <group_members>

        '''

        group0 = "ref"
        group1 = "ref"

        # group:members
        group_members_str = ""

        #log.debug("1. self.conf_dict['a_sample']['param']={}".format(
        #    self.conf_dict['a_sample']['param']))
        #log.debug("1. self.conf_dict['b_sample']['param']={}".format(
        #    self.conf_dict['b_sample']['param']))

        #log.debug("2. self.group_members_vcf_str={}".format(
        #    self.group_members_vcf_str))

        #log.debug("3. self.conf_dict['group_members']['param']={}".format(
        #    self.conf_dict['group_members']['param']))
        #log.debug("4. self.conf_dict['group_members']['ini']={}".format(
        #    self.conf_dict['group_members']['ini']))

        # It's easy mode if a_sample or b_sample is included

        if self._is_easy_mode() == True:

            a_members = ""
            b_members = ""

            # for group
            if self.conf_dict['a_sample']['param'] != "":
                group0 = "a"
                a_members = "{},".format(
                    self.conf_dict['a_sample']['param'])

            if self.conf_dict['b_sample']['param'] != "":
                group1 = "b"
                b_members = "{},".format(
                    self.conf_dict['b_sample']['param'])

            group_members_str = "{}:{}{}:{}".format(
                group0, a_members, group1, b_members)
            group_members_str = re.sub(r",$", "", group_members_str)

        else:
            group_members_str = self._value_choice('group_members')

            # If neither param nor ini is specified, refer to the vcf file
            if group_members_str == "":
                group_members_str = self.group_members_vcf_str

        #log.debug("group_members_str={}".format(group_members_str))

        return group_members_str


    def set_distinguish_groups_str(self):
        ''' distinguish_groups
        Variable origin priority
        1. command parameter        <a_sample, b_sample>
        2. command parameter        <distinguish_groups>
        3. ini file                 <distinguish_groups>
        '''

        group0 = "ref"
        group1 = "ref"
        region_name = ""
        pick_mode = ""
        indel_size = ""

        # group0/group1:region_name[+region_name...]:pick_mode:indel_size
        distinguish_groups_str = ""

        #log.debug("1. self.conf_dict['a_sample']['param']={}".format(
        #    self.conf_dict['a_sample']['param']))

        #log.debug("1. self.conf_dict['b_sample']['param']={}".format(
        #    self.conf_dict['b_sample']['param']))

        #log.debug(
        #    "2. self.conf_dict['distinguish_groups']['param']={}".format(
        #    self.conf_dict['distinguish_groups']['param']))

        #log.debug(
        #    "3. self.conf_dict['distinguish_groups']['ini']={}".format(
        #    self.conf_dict['distinguish_groups']['ini']))

        # It's easy mode if a_sample or b_sample is included
        if self._is_easy_mode() == True:

            # for group

            if utl.is_Not_NoneAndNull(
                self.conf_dict['a_sample']['param']):
                group0 = "a"

            if utl.is_Not_NoneAndNull(
                self.conf_dict['b_sample']['param']):
                group1 = "b"

            # EASY_MODE will format variables later
            region_name = "<EASY_MODE>"

            # pick_mode
            pick_mode = self.pick_mode

            # indel_size
            indel_size = self.indel_size

            distinguish_groups_str = "{}/{}:{}:{}:{}".format(
                group0, group1, region_name, pick_mode, indel_size)

        else:
            distinguish_groups_str = self._value_choice('distinguish_groups')

        #log.debug("distinguish_groups_str={}".format(distinguish_groups_str))

        return distinguish_groups_str


    def set_regions_dict(self, regions_str):
        ''' make regions dict from regions_str
        '''
        rectified_regions_str = ""
        regions_dict = dict()
        region_name_list = list()

        log.info("original regions_str={}".format(regions_str))

        # Revive the omitted field
        for region in regions_str.split(','):

            field = region.split(':')
            field_cnt = len(field)
            rectified_region = ""

            err_mes = ""

            if field_cnt == 1:
                err_mes = "You cannot define just only the region name."

            else:
                if field_cnt == 2:
                    # Blank has the same meaning as 'whole'
                    field.append("whole")
                    field_cnt += 1

                if field_cnt == 3:
                    # Check the member with the group name with vcf
                    region_name = field[0]
                    chrom_name = field[1]
                    rg_str = field[2]

                    #print("{}, {}, {}".format(
                    #    region_name, chrom_name, rg_str))

                    # 1) chrom_name
                    if not self._is_chrom_name(chrom_name):
                        #print("{} is bad".format(chrom_name))
                        err_mes = "This chrom_name ({}) ".format(chrom_name)
                        err_mes += "is incorrect."

                    # 2) rg_str
                    elif rg_str == 'whole':
                        region_str, start, end, length = \
                            self._get_chrom_info(chrom_name)
                        rg_str = "{}-{}".format(start, end)

                    elif not self._is_valid_int_range(rg_str):
                        err_mes = "This range ({}) ".format(rg_str)
                        err_mes += "is incorrect."

                    # 2) min, max on chromosome
                    elif not self._is_valid_chrom_range(
                        "{}:{}".format(chrom_name, rg_str)):
                        err_mes = "This range ({}) ".format(rg_str)
                        err_mes += "is beyond the range of the chromosome."

                    rectified_region = "{}:{}:{}".format(
                        region_name, chrom_name, rg_str)

                    #print("ok {} {} {}".format(field[0], field[1], field[2]))

                else:
                    err_mes = "The number of fields exceeds "
                    err_mes += "three({}).".format(field_cnt)

            if err_mes != "":
                log.error("{} exit.".format(err_mes))
                log.error("region=\'{}\'.".format(region))
                log.error("regions_str=\'{}\'.".format(regions_str))
                sys.exit(1)
                
            rectified_regions_str += "{},".format(rectified_region)

        rectified_regions_str = re.sub(r",$", "", rectified_regions_str)
        log.info("rectified_regions_str={}".format(rectified_regions_str))


        # set into dict
        for region in rectified_regions_str.split(','):

            region_name, chrom_name, r_range = region.split(':')

            start_pos = 0
            end_pos = 0

            region_str = chrom_name
            if '-' in r_range:
                start_pos, end_pos = r_range.split('-')
                region_str = "{}:{}".format(chrom_name, r_range)

            # ---------------------------------------------
            regions_dict[region_name] = \
                {   'chr': chrom_name,
                    'start': start_pos,
                    'end': end_pos,
                    'reg': region_str
                }

            #log.info("{}".format(regions_dict))

        # region_name_list is not include chrom_name
        region_name_list = list(regions_dict.keys())

        # regions_dict include chrom_name
        for chrom_name in self.ref_fasta_chrom_list:

            region_str, start_pos, end_pos, length = \
                self._get_chrom_info(chrom_name)

            regions_dict[chrom_name] = \
                {   'chr': chrom_name,
                    'start': start_pos,
                    'end': end_pos,
                    'reg': region_str
                }

        #print("\nregion_name_list={}\n".format(region_name_list))
        return \
            rectified_regions_str, \
            regions_dict, \
            region_name_list


    def set_group_members_dict(self, group_members_str):
        '''
        '''
        # easy a:s1,s2,s3,b:s4,s5,s6
        #      g0:s1,s4,s7,g1:s3,s5,s9

        group_members_dict = dict()

        # separate
        group_name = ""

        for member in group_members_str.split(","):
            if ":" in member:
                group_name, member = member.split(":")

                err_mes = ""

                # group name check
                if "ref" == group_name.lower():
                    err_mes = "ref cannot be used as a group name."

                elif self._is_chrom_name(group_name):
                    err_mes = "chrom name is a reserved word, so it cannot "
                    err_mes += "be used as a group name."

                if err_mes != "":
                    log.error("{} exit.".format(err_mes))
                    log.error("group_members_str=\'{}\'.".format(
                    group_members_str))
                    sys.exit(1)

                group_members_dict[group_name] = list()

            # for easy mode
            # If only A or B, it is a comparison with ref.
            if member == "None":
                member = "ref"

            # check members
            if not utl.is_sample_name(member):
                err_mes = "The sample name \'{}\' does not ".format(member)
                err_mes += "match either the nickname, basename, "
                err_mes += "or fullname in the vcf file."

                log.error("{} exit.".format(err_mes))
                log.error("group_members_str=\'{}\'.".format(
                    group_members_str))
                sys.exit(1)

            group_members_dict[group_name].append(member)

        group_name_list = list(group_members_dict.keys())

        # at last, append 'ref'
        group_members_dict['ref'] = list()
        group_members_dict['ref'].append('ref')

        return group_members_dict, group_name_list


    def set_distinguish_groups_list(self, distinguish_groups_str):
        '''
        '''

        # group0/group1
        # region_name[+region_name...]
        # pick_mode
        # indel_size

        rectified_distinguish_groups_str = ""

        # distinguish_groups_str is the data
        # that does not fill in the missing parts in the 4 fields
        for distinguish_group in distinguish_groups_str.split(','):

            #print("distinguish_group={}".format(distinguish_group))

            field = distinguish_group.split(':')
            field_cnt = len(field)

            group_pair = ""
            region_names = ""
            pick_mode = self.pick_mode
            indel_size = self.indel_size

            err_mes = ""
            # 1. check field_cnt (4)
            if field_cnt == 1:
                group_pair = distinguish_group

            elif field_cnt == 2:
                group_pair, region_names = distinguish_group.split(':')

            elif field_cnt == 3:
                group_pair, region_names, pick_mode = \
                    distinguish_group.split(':')

            elif field_cnt == 4:
                group_pair, region_names, pick_mode, indel_size = \
                    distinguish_group.split(':')

            else:
                err_mes = "The number of fields exceeds "
                err_mes += "four({}).".format(field_cnt)

            #print("field_cnt={}".format(field_cnt))
            #print("{} {} {} {}".format(group_pair, region_names,
            #    pick_mode, indel_size))

            # 2. check content of group_pair
            if err_mes == "":
                if group_pair == "":
                    err_mes = "Group pairs must always be specified."

                elif "/" not in group_pair:
                    err_mes = "Group pairs must be separated by slash." 

            if err_mes != "":
                log.error("{} exit.".format(err_mes))
                log.error("distinguish_groups_str=\'{}\'.".format(
                    distinguish_groups_str))
                sys.exit(1)


            err_mes = ""
            # 3. check group0, group1
            group0, group1 = group_pair.split('/')

            group0_ok = 0
            group1_ok = 0

            # if group0 and group1 have the same name.
            if group0 == group1:
                err_mes = "Two groups have the same name."
            
            else:
                # according to each content
                if group0.lower() == 'ref':
                    group0_ok = 1

                if group1.lower() == 'ref':
                    group1_ok = 1

                if self._is_group_name(group0):
                    group0_ok = 1

                if self._is_group_name(group1):
                    group1_ok = 1

            if group0_ok == 0 or group1_ok == 0:
                err_mes = "Group names are not defined except \'ref\'."

            if err_mes != "":
                log.error("{} exit.".format(err_mes))
                log.error("distinguish_groups_str=\'{}\'.".format(
                    distinguish_groups_str))
                log.error("group_name_list=\'{}\'.".format(
                    ", ".join(self.group_name_list)))
                sys.exit(1)


            err_mes = ""
            # 4. check region_name
            if region_names == "":
                region_names = 'whole'

            # If the group name is 'whole', connect the chromosome names
            # with '+' 

            if region_names == 'whole':
                region_names = '+'.join(self.ref_fasta_chrom_list)
            else:
                # The region name is a variable or chromosome name
                # separated by '+'
                for region_name in region_names.split('+'):
                    region_name_ok = 0
                    chrom_name_ok = 0

                    if self._is_region_name(region_name):
                        region_name_ok = 1

                    if self._is_chrom_name(region_name):
                        chrom_name_ok = 1

                    if region_name_ok == 0 and chrom_name_ok == 0:
                        err_mes = "region name {} ".format(region_name)
                        err_mes += "is not in region names or "
                        err_mes += "chromosome name."

                    if err_mes != "":
                        log.error("{} exit.".format(err_mes))
                        log.error("distinguish_groups_str=\'{}\'.".format(
                            distinguish_groups_str))
                        log.error("group_name_list=\'{}\'.".format(
                            ", ".join(self.group_name_list)))
                        log.error("ref_fasta_chrom_list=\'{}\'.".format(
                            ", ".join(self.ref_fasta_chrom_list)))
                        sys.exit(1)


            # 5. check region_name
            if pick_mode == "":
                pick_mode = self.pick_mode

            if pick_mode not in glv.pick_mode_list:
                err_mes = "The pick mode ({}) must be ".format(pick_mode)
                err_mes += "one of these {}.".format(self_mode_list)
                log.error("{} exit.".format(err_mes))
                log.error("distinguish_groups_str=\'{}\'.".format(
                    distinguish_groups_str))
                sys.exit(1)

            # 6. check indel_size
            if indel_size == "":
                indel_size = self.indel_size

            if not self._is_valid_int_range(indel_size):
                err_mes = "indel_size ({}) must have ".format(indel_size)
                err_mes += "two numbers separated by a \'-\' "
                err_mes += "and min <= max."
                log.error("{} exit.".format(err_mes))
                log.error("distinguish_groups_str=\'{}\'.".format(
                    distinguish_groups_str))
                sys.exit(1)


            distinguish_group = "{}:{}:{}:{}".format(
                group_pair, region_names, pick_mode, indel_size)

            rectified_distinguish_groups_str += "{},".format(
                distinguish_group)


        rectified_distinguish_groups_str = re.sub(
            r",$", "", rectified_distinguish_groups_str)

        #print("{}".format(rectified_distinguish_groups_str))
        #sys.exit(1)

        # string to dict
        distinguish_groups_list = list()

        for disting_str in rectified_distinguish_groups_str.split(','):
            group_pair, region_names, pick_mode, indel_size \
                = disting_str.split(':')
            group0, group1 = group_pair.split("/")
            region_list = region_names.split("+")

            #print("{} == {} == {} == {}".format(
            #    group_pair, region_names, pick_mode, indel_size))

            distins_dict = {
                0: group0,
                1: group1,
                'disting_str': disting_str,
                'regions': region_list,
                'pick_mode': pick_mode,
                'indel_size': indel_size}

            distinguish_groups_list.append(distins_dict)

        return \
            rectified_distinguish_groups_str, \
            distinguish_groups_list

