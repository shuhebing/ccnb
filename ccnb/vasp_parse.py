#from pymatgen.io.vasp.inputs import PotcarSingle
import pymatgen.io.vasp.inputs as vin
import re

class parse_potcar(object):
    parse_functions = {"LULTRA": vin.parse_bool,
                       "LUNSCR": vin.parse_bool,
                       "LCOR": vin.parse_bool,
                       "LPAW": vin.parse_bool,
                       "EATOM": vin.parse_float,
                       "RPACOR": vin.parse_float,
                       "POMASS": vin.parse_float,
                       "ZVAL": vin.parse_float,
                       "RCORE": vin.parse_float,
                       "RWIGS": vin.parse_float,
                       "ENMAX": vin.parse_float,
                       "ENMIN": vin.parse_float,
                       "EMMIN": vin.parse_float,
                       "EAUG": vin.parse_float,
                       "DEXC": vin.parse_float,
                       "RMAX": vin.parse_float,
                       "RAUG": vin.parse_float,
                       "RDEP": vin.parse_float,
                       "RDEPT": vin.parse_float,
                       "QCUT": vin.parse_float,
                       "QGAM": vin.parse_float,
                       "RCLOC": vin.parse_float,
                       "IUNSCR": vin.parse_int,
                       "ICORE": vin.parse_int,
                       "NDATA": vin.parse_int,
                       "VRHFIN": vin.parse_string,
                       "LEXCH": vin.parse_string,
                       "TITEL": vin.parse_string,
                       "STEP": vin.parse_list,
                       "RRKJ": vin.parse_list,
                       "GGA": vin.parse_list}
    def __init__(self, data):
        search_lines = re.findall(r"(?s)(parameters from PSCTR are:"
                                         r".*?END of PSCTR-controll parameters)",
                                         data)#.group(2)

        self.all_keywords = []
        for i in range(len(search_lines)):
            keywords={}
            for key, val in re.findall(r"(\S+)\s*=\s*(.*?)(?=;|$)",
                                       search_lines[i], flags=re.MULTILINE):
                try:
                    keywords[key] = self.parse_functions[key](val)
                except KeyError:
                    warnings.warn("Ignoring unknown variable type %s" % key)
            self.all_keywords.append(keywords)

    def get_max_enmax(self):
        max_enmax=0
        for i in range(len(self.all_keywords)):
            enmax=self.all_keywords[i]['ENMAX']
            if max_enmax < enmax:
                max_enmax=enmax
        return max_enmax

# p=PotcarSingle(file.read())
# print(p.enmax)
