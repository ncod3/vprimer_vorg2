# -*- coding: utf-8 -*-
# https://qiita.com/podhmo/items/1eb7e72a47b713c9cda2
# https://qiita.com/kzkadc/items/e4fc7bc9c003de1eb6d0

import sys
import os
import errno
import time

# global configuration
import vprimer.glv as glv

import argparse
from multiprocessing import Pool
import multiprocessing as multi
from vprimer.__init__ import __version__
from vprimer.logging_config import LogConf


class Param(object):

    def __init__(self):

        pass

    def open_log(self):

        global log
        log = LogConf.open_log(__name__)


    def get_args(self):

        parser = self._get_options()

        # self.p is dict
        if len(sys.argv) == 1:
            self.p = vars(parser.parse_args(['-h']))
        else:
            # to dict
            self.p = vars(parser.parse_args())

        return self


    def _get_options(self):

        parser = argparse.ArgumentParser(
            description='vprimer version {}'.format(__version__),
            formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('vprimer ...\n')

        parser.add_argument('-V', '--version', action='version',
            version='$(prog)s {}'.format(__version__))

        #---------------------------------------------------------
        parser.add_argument('-w', '--without_group', action='store',
            type=str, metavar='',
            help="Evaluate markers on all valid variants without" +\
                " specifying a group")

        #---------------------------------------------------------
        # -o --out_dir
        parser.add_argument('-o', '--out_dir', action='store',
            type=str, metavar='', default='out_vprimer',
            help="dir name for data output")

        # -i --ini
        parser.add_argument('-i', '--ini_file', action='store',
            type=str, metavar='',
            help="ini file [optional]")

        # -s --show_samples
        parser.add_argument('-s', '--show_samples', action='store_true',
            #type=str, metavar='',
            help="show sample names embedded in VCF files")

        # -s --show_fasta
        parser.add_argument('-c', '--show_fasta', action='store_true',
            #type=str, metavar='',
            help="show fasta chromosomal information.")

        # -v --vcf
        parser.add_argument('-v', '--vcf', action='store',
            type=str, metavar='',
            help="[required] vcf file (text or gz)")

        # -r --ref
        parser.add_argument('-r', '--ref', action='store',
            type=str, metavar='',
            help="[required] reference fasta (txt or gz)")

        # -t --thread
        parser.add_argument('-t', '--thread', action='store',
            type=int, metavar='',
            help="thread number")

        # -j --use_joblib_threading
        parser.add_argument('-j', '--use_joblib_threading', action='store',
            type=str, metavar='', default='yes',
            help="use or not threading yes/no default yes")

        # -p --pick_mode
        parser.add_argument('-p', '--pick_mode', action='store',
            type=str, metavar='', default='all',
            help="variant pick mode : indel / snp / all")

        #-------------------------------------------------------------
        # -z --indel_size
        parser.add_argument('-z', '--indel_size', action='store',
            type=str, metavar='',
            help="target indel size, min-max")

        # -d --product_size
        parser.add_argument('-d', '--product_size', action='store',
            type=str, metavar='',
            help="PCR product size, min-max")
        #-------------------------------------------------------------

        # -e --enzyme_file
        parser.add_argument('-e', '--enzyme_file', action='store',
            type=str, metavar='',
            nargs='*',
            help="enzyme file list, separate by comma")

        # -e --enzyme
        parser.add_argument('-E', '--enzyme', action='store',
            type=str, metavar='',
            nargs='*',
            help="enzyme name list")

        #-------------------------------------------------------------
        # Easy mode
        # -T --target
        parser.add_argument('-T', '--target', action='store',
            type=str, metavar='',
            nargs='*',
            help="target region, chrom:stt-end")

        # -A --sample_a
        parser.add_argument('-A', '--a_sample', action='store',
            type=str, metavar='',
            nargs='*',
            help="[required] group A sample names")

        # -B --sample_b
        parser.add_argument('-B', '--b_sample', action='store',
            type=str, metavar='',
            nargs='*',
            help="[required] group B sample names")
        #-------------------------------------------------------------

        # -f --p3_params
        parser.add_argument('-f', '--p3_params', action='store',
            type=str, metavar='',
            help="primer3 parameter file")

        # ---- full mode
        # ===================================================================
        parser.add_argument('-R', '--regions', action='store',
            type=str, metavar='',
            nargs='*',
            help="regions full description")

        parser.add_argument('-D', '--distinguish_groups', action='store',
            type=str, metavar='',
            nargs='*',
            help="distinguish_groups full description")

        parser.add_argument('-G', '--group_members', action='store',
            type=str, metavar='',
            nargs='*',
            help="group_members full description")
        # ===================================================================


        # -x --progress
        parser.add_argument('-x', '--progress', action='store',
            type=str, metavar='', default='all',
            help="progress start point")

        # -y --stop
        parser.add_argument('-y', '--stop', action='store',
            type=str, metavar='', default='no',
            help="progress stop point, prepare/variant/marker/primer")

        # -u --analyse_caps
        parser.add_argument('-u', '--analyse_caps', action='store_true',
            help="print caps info")

        return parser
