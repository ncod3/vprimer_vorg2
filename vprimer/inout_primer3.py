# -*- coding: utf-8 -*-

import sys
import os
import errno
import re
import pprint

import logging
log = logging.getLogger(__name__)

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

class InoutPrimer3(object):

    def __init__(self):

        # self.p3_header_list = glv.conf.primer3_header_dict
        self.primer3_header_dict = glv.conf.primer3_header_dict

        self.p3_input_dict = {
            'P3_COMMENT' : '',
            'SEQUENCE_TARGET': '',
            'SEQUENCE_EXCLUDED_REGION': '',
            'SEQUENCE_ID' : '',
            'SEQUENCE_TEMPLATE' : ''}

        self.p3_out = ''
        self.p3_output_dict = dict()
        # reserve
        self.p3_output_dict['PRIMER_PAIR_NUM_RETURNED'] = 0
        self.p3_output_dict['PRIMER_ERROR'] = ''

        self.PRIMER_PAIR_NUM_RETURNED = \
            self.p3_output_dict['PRIMER_PAIR_NUM_RETURNED']
        self.PRIMER_ERROR = \
            self.p3_output_dict['PRIMER_ERROR']

        self.primer_region = ''
        self.left_fasta_id = ''
        self.right_fasta_id = ''


    #--------------------------------------------
    # output

    def get_primer_left_id(self):
        return self.left_fasta_id

    def get_primer_right_id(self):
        return self.right_fasta_id

    def set_primer_name(self, left_fasta_id, right_fasta_id):
        self.left_fasta_id = left_fasta_id
        self.right_fasta_id = right_fasta_id

    def get_primer_product_size(self):
        return self.p3_output_dict['PRIMER_PAIR_0_PRODUCT_SIZE']

    def get_primer_left_info(self):
        return map(int, self.p3_output_dict['PRIMER_LEFT_0'].split(','))

    def get_primer_right_info(self):
        return map(int, self.p3_output_dict['PRIMER_RIGHT_0'].split(','))

    def get_primer_region(self):
        return self.primer_region

    def get_primer_left_seq(self):
        return self.p3_output_dict['PRIMER_LEFT_0_SEQUENCE']

    def get_primer_right_seq(self):
        return self.p3_output_dict['PRIMER_RIGHT_0_SEQUENCE']

    def get_primer_left(self):
        return self.p3_output_dict['PRIMER_LEFT_0']

    def get_primer_right(self):
        return self.p3_output_dict['PRIMER_RIGHT_0']

    def get_p3_comment(self):
        return self.p3_output_dict['P3_COMMENT']

    def get_sequence_id(self):
        return self.p3_output_dict['SEQUENCE_ID']

    def get_primer3_out(self):
        return self.p3_out


    def set_primer3_out(self, p3_out):
        '''
        '''

        self.p3_out = p3_out
        primer_region_list = list()

        for item in p3_out.split('\n'):
            if item == '=' or item == '':
                continue
            tag, value = item.split('=')
            self.p3_output_dict[tag] = value

            if tag == 'PRIMER_PAIR_NUM_RETURNED':
                self.PRIMER_PAIR_NUM_RETURNED = int(value)

            elif tag == 'PRIMER_ERROR':
                self.PRIMER_ERROR = str(value)

            elif tag == 'PRIMER_LEFT_0':
                primer_region_list.append(value)
            elif tag == 'PRIMER_RIGHT_0':
                end_rel_pos, length = map(int, value.split(','))
                # 525,5 => 521,5
                # 521...5
                plus_rel_pos = end_rel_pos - length + 1
                primer_region_list.append("{},{}".format(
                    plus_rel_pos, length))

            if len(primer_region_list) != 0:
                self.primer_region = ' '.join(primer_region_list)

    #--------------------------------------------
    # input
    def get_sequence_excluded_region(self):
        return self.p3_input_dict['SEQUENCE_EXCLUDED_REGION']

    def set_p3_comment(self, p3_comment):
        self.p3_input_dict['P3_COMMENT'] = p3_comment

    def set_sequence_target(self, sequence_target):
        self.p3_input_dict['SEQUENCE_TARGET'] = sequence_target

    def add_ex_region(self, any_excluded_region):

        if self.p3_input_dict['SEQUENCE_EXCLUDED_REGION'] == '':
            ex_region = [any_excluded_region]
        else:
            ex_region = [
                self.p3_input_dict['SEQUENCE_EXCLUDED_REGION'],
                any_excluded_region]

        self.p3_input_dict['SEQUENCE_EXCLUDED_REGION'] = \
            ' '.join(ex_region)

    def set_sequence_id(self, sequence_id):
        self.p3_input_dict['SEQUENCE_ID'] = sequence_id

    def set_sequence_template(self, sequence_template):
        self.p3_input_dict['SEQUENCE_TEMPLATE'] = sequence_template


    def get_p3_input(self):

        # glv.conf.self.primer3_header_dict
        p3_input_list = list()

        # convert global dictionary to list
        for p3_key in glv.conf.primer3_header_dict.keys():
            p3_input_list += ["{}={}".format(
                p3_key, glv.conf.primer3_header_dict[p3_key])]

        #self.p3_input_dict = {
        #    'P3_COMMENT' : '',
        #    'SEQUENCE_TARGET': '',
        #    'SEQUENCE_EXCLUDED_REGION': '',
        #    'SEQUENCE_ID' : '',
        #    'SEQUENCE_TEMPLATE' : ''}

        # add nessasary key=value
        p3_input_list += ['{}={}'.format(
            'P3_COMMENT',
            self.p3_input_dict['P3_COMMENT'])]
        p3_input_list += ['{}={}'.format(
            'SEQUENCE_TARGET',
            self.p3_input_dict['SEQUENCE_TARGET'])]
        p3_input_list += ['{}={}'.format(
            'SEQUENCE_EXCLUDED_REGION',
            self.p3_input_dict['SEQUENCE_EXCLUDED_REGION'])]
        p3_input_list += ['{}={}'.format(
            'SEQUENCE_ID',
            self.p3_input_dict['SEQUENCE_ID'])]
        p3_input_list += ['{}={}'.format(
            'SEQUENCE_TEMPLATE',
            self.p3_input_dict['SEQUENCE_TEMPLATE'])]
        p3_input_list += ['=\n']

        #pprint.pprint(p3_input_list)
        #txt = '\n'.join(map(str, p3_input_list))

        return '\n'.join(map(str, p3_input_list))


