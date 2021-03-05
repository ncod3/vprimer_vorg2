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

class ConfCurrSet(object):

    def open_log_currset(self):

        global log
        log = LogConf.open_log(__name__)


    def out_current_settings(self):
        ''' Output to a file with config (ini format)
        '''

        current_setting_ini = list()
        whole_command_line = ' '.join(sys.argv)

        # [vprimer]
        current_setting_ini.append("{}".format(glv.ini_section))

        # date
        date_stamp = "\n# {}".format(glv.now_datetime_form)
        current_setting_ini.append(date_stamp)

        # whole_command_line
        whole_command_line = "\n# {}".format(whole_command_line)
        current_setting_ini.append(whole_command_line)

        current_setting_ini.append("\n#")

        for vname in self.conf_dict.keys():

            if 'chosen' in self.conf_dict[vname]:
                key_value = "{} = {}".format(
                    vname, self.conf_dict[vname]['chosen'])

                current_setting_ini.append(key_value)

                if vname == "ref" or vname == "stop" or \
                    vname == "product_size" or vname == "enzyme" or \
                    vname == "group_members" or vname == "blast_distance" or \
                    vname == "use_joblib_threading":

                    current_setting_ini.append("\n#")

        # exist or not, self.curr_setting_file_path
        if os.path.isfile(self.curr_setting_file_path):
            # If the file exists, move it to bak
            log.info("found {}".format(self.curr_setting_file_path))
            utl.save_to_tmpfile(self.curr_setting_file_path)
        else:
            log.info("not found {}".format(self.curr_setting_file_path))

        # write to sample_name_file
        with open(self.curr_setting_file_path, mode='w') as f:
            # Export while adjusting
            #line = self._convert_setting_ini(current_setting_ini)
            #f.write("{}\n".format("\n".join(current_setting_ini)))
            line = self._convert_setting_ini(current_setting_ini)
            f.write("{}\n".format(line))

        log.info("save {}".format(self.curr_setting_file_path))

        # ====
        log.info("self.conf_dict=\n{}".format(
            pprint.pformat(self.conf_dict)))

        log.info("self.regions_dict=\n{}".format(
            pprint.pformat(self.regions_dict)))
        log.info("self.group_members_dict=\n{}".format(
            pprint.pformat(self.group_members_dict)))
        log.info("self.distinguish_groups_list=\n{}".format(
            pprint.pformat(self.distinguish_groups_list)))


    def _convert_setting_ini(self, current_setting_ini):
        '''
        '''

        oline = list()

        for line in current_setting_ini:
            if line.startswith("regions") or \
                line.startswith("distinguish_groups") or \
                line.startswith("group_members"):

                line = self._format_distin(line)

            elif "," in line:
                line = re.sub(r",", ", ", line)
            oline.append(line)
            
        return "\n".join(oline)


    def _format_distin(self, line):

        # insert cr after =
        line = re.sub(r" = ", " = \n    ", line) 
        # separate value to 2nd line, etc
        line = re.sub(r"(,)([^,:]+)(:)", r"\n    \2\3", line) 

        # space both side / :
        line = re.sub(r"/"," / ", line)
        line = re.sub(r":", " : ", line) 
        # space one side ,
        line = re.sub(r",", ", ", line) 

        line += "\n"
        return line




