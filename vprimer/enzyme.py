# -*- coding: utf-8 -*-

# http://rebase.neb.com/rebase/rebase.enz.html
# https://international.neb.com/tools-and-resources/selection-charts/isoschizoomers

import sys
import os
import errno

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import Bio.Restriction.Restriction_Dictionary as ResDict

class Enzyme(object):

    def __init__(self):

        self.enzyme_list = list()


    def read_enzyme_file(self):

        enz_list = list()

        for enzyme_file in glv.conf.enzyme_files_list:
            log.info("{}".format(enzyme_file))
            with open(enzyme_file, "r") as f:
                # iterator
                for r_liner in f:
                    r_line = r_liner.strip()    # cr, ws

                    if r_line.startswith('#') or r_line == '':
                        continue
                    r_line = utl.strip_hash_comment(r_line)

                    if r_line not in ResDict.rest_dict:
                        log.critical("your ENZYME {} is not in list.".format(
                            r_line))
                        sys.exit(1)

                    enz_list.append(r_line)

        org_cnt = len(enz_list)
        self.enzyme_list = list(set(enz_list))
        set_cnt = len(self.enzyme_list)

        #print("{} {}".format(org_cnt, set_cnt))

        if org_cnt != set_cnt:
            diff_cnt = org_cnt - set_cnt
            log.info("There were {} duplicate enzymes. removed.".format(
                diff_cnt))

        self.enzyme_list.sort(key=None, reverse=False)
        log.info("total enzyme cnt={}, {}".format(
            len(self.enzyme_list), self.enzyme_list))

    @classmethod
    def prepare_enzyme_file(cls):

        #enzyme_files=
        #    data/DaiichiSeiyaku_NEB_07.txt;\
        #    data/RE_labo_v1.0_20150403.txt;
        #    data/Takara07.txt

        log.info("glv.conf.enzyme_files_user_list={}".format(
            glv.conf.enzyme_files_user_list))

        for enzyme_file_user in glv.conf.enzyme_files_user_list:

            # only slink don't use
            basename_user = os.path.basename(enzyme_file_user)
            enzyme_file_slink_system = "{}/{}{}".format(
                glv.conf.ref_dir_path, basename_user, '.org_slink')
            # save as conf global
            glv.conf.enzyme_files_list.append(enzyme_file_slink_system)

            # not exist link file
            # exist link file, exist actual file
            # exist link file, not exist actual file
            if os.path.isfile(enzyme_file_slink_system):
                log.info("found {}.".format(enzyme_file_slink_system))
            else:
                log.info("not found {}.".format(enzyme_file_slink_system))
                utl.ln_s(enzyme_file_user, enzyme_file_slink_system)

        log.info("glv.conf.enzyme_file_list {}".format(
            glv.conf.enzyme_files_list))

