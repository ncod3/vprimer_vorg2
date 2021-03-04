# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import pprint

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import pandas as pd
from joblib import Parallel, delayed

#from vprimer.enzyme import Enzyme
from vprimer.eval_variant import EvalVariant


class Marker(object):

    def __init__(self):

        self.enzyme_name_list = list()

    def design_marker(self):

        self.enzyme_name_list = glv.conf.enzyme_name_list

        proc_name = "marker"
        log.info("-------------------------------")
        log.info("Start processing {}\n".format(proc_name))

        # stop, action, gothrough
        ret_status = utl.decide_action_stop(proc_name)

        if ret_status == "stop":
            msg = "STOP. "
            msg += "Current process \'{}\' ".format(proc_name)
            msg += "has exceeded the User-specified stop point "
            msg += "\'{}', ".format(glv.conf.stop)
            msg += "so stop program. exit."
            log.info(msg)
            sys.exit(1)


        elif ret_status == "gothrough":
            msg = "SKIP \'{}\' proc, ".format(proc_name)
            msg += "glv.conf.progress = {}, ".format(glv.conf.progress)
            msg += "glv.conf.stop = {}, ".format(glv.conf.stop)
            msg += "so skip program."
            log.info(msg)
            return


        # Design a fragment sequence for primer3
        for proc_cnt, distin_dict in enumerate(glv.outlist.distin_files, 1):

            # logging current target
            utl.print_distin_info("marker", distin_dict, proc_cnt)

            # read variant file 
            variant_file = distin_dict['variant']['out_path']
            log.info("variant_file {}".format(variant_file))

            df_distin = pd.read_csv(
                variant_file, sep='\t', header=0, index_col=None)

            # file name to write out result to text
            out_txt_file = distin_dict['marker']['out_path']
            utl.save_to_tmpfile(out_txt_file)

            start = time.time()
            with open(out_txt_file, mode='a') as f:

                ''' eval_variant.py
                class EvalVariant(object):
                def _check_effect_of_enzyme(
                    self, seq_target, enzyme_name_list):
                    http://biopython.org/DIST/docs/cookbook/Restriction.html
                    biopython <= 1.76 for IUPACAmbiguousDNA()

                    multi_site_seq = Seq(seq_target, IUPACAmbiguousDNA())
                    rb = Restriction.RestrictionBatch(enzyme_name_list)
                    Analong = Restriction.Analysis(rb, multi_site_seq)
                    caps_ResTyp_dict = Analong.with_sites()

                This RestrictionBatch method sometimes returned slightly
                inaccurate results when executed in parallel.
                Therefore, parallel is not used now.
                '''

                #if glv.conf.parallel == True:
                if False:
                    log.info("do Parallel cpu {} parallel {}".format(
                        glv.conf.thread,
                        glv.conf.parallel_full_thread))

                    Parallel(
                        n_jobs=glv.conf.parallel_full_thread,
                        backend="threading")(
                        [
                            delayed(self._loop_evaluate_for_marker)
                                (distin_dict, variant_df_row, f) \
                                for variant_df_row in df_distin.itertuples()
                        ]
                    )

                else:
                    log.info("do Serial cpu 1")

                    # each variant
                    for variant_df_row in df_distin.itertuples():

                        # Determine if the variant can be used as a marker.
                        # For those that can be marked, prepare the
                        # information for primer3.
                        self._loop_evaluate_for_marker(
                            distin_dict, variant_df_row, f)

            utl.sort_file(
                'marker', distin_dict, out_txt_file,
                'chrom', 'pos', 'marker_info', 'string')

            log.info("marker {} > {}.txt\n".format(
                utl.elapsed_time(time.time(), start),
                distin_dict['marker']['base_nam']))


    def _loop_evaluate_for_marker(self, distin_dict, variant_df_row, f):
        '''
        '''

        # In order to perform parallel processing safely,
        # an instance is created for each processing unit here.
        evalv = EvalVariant(self.enzyme_name_list)
        evalv.evaluate_for_marker(variant_df_row, distin_dict)
        
        # indel or CAPS marker by SNP with restriction enzyme site effective
        if evalv.marker_available == True:

            # Create seq_template_ref if markerable
            evalv.make_seq_template_ref()

            # Lines are copied as many as the number of effective
            # restriction enzymes.
            evalv.copy_line_for_effective_restriction_enzymes()

            # write out to file
            f.write("{}\n".format(evalv.line))
            evalv.line = ''

