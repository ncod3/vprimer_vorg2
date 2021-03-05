#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import pprint

# global variables declaration
import vprimer.glv as glv
import vprimer.utils as utl
from vprimer.logging_config import LogConf

# --- read conf and param
glv.init('vprimer')
log = LogConf.open_log(__name__)

# using Class
from vprimer.variant import Variant
from vprimer.marker import Marker
from vprimer.primer import Primer
from vprimer.formtxt import FormTxt

def main():

    log.info('program started at {}'.format(glv.now_datetime_str))

    # run
    vpr = VPrimer()
    vpr.run()

    log.info("program finished {}\n".format(
        utl.elapsed_time(time.time(), glv.now_epochtime)))


class VPrimer(object):

    def __init__(self):

        self.variant = Variant()
        self.marker = Marker()
        self.primer = Primer()
        self.formtxt = FormTxt()

    def run(self):

        self.prepare()

        #self.variant.print_allele_int()
        #self.variant.print_all_allele_int()
        #glv.conf.get_vcf_pos_info()
        if glv.conf.show_genotype != "no":
            self.variant.print_allele()
            log.info("program finished {}\n".format(
                utl.elapsed_time(time.time(), glv.now_epochtime)))
            sys.exit(1)

        # variant
        self.variant.pick_variant()
        # marker
        self.marker.design_marker()
        # primer
        self.primer.construct_primer()
        # format
        self.formtxt.format_text()


    def prepare(self):

        # Determine the value of a variable according to the priority
        # of the parameter

        # Choice of variables by priority
        glv.conf.choice_variables()

        # read reference into global environment
        glv.ref = glv.ref.prepare_ref()

        # setup all variables
        glv.conf.setup_variables()

        # Write the current settings to a file
        glv.conf.out_current_settings()

        # prepare output text by glv.conf.distin_g_list
        glv.outlist.prepare_distin_files()

if __name__ == '__main__':
    main()


