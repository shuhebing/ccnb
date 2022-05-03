#from pymatgen.io.vasp.inputs import PotcarSingle
import pymatgen.io.vasp.inputs as vin
import re


class parse_potcar(object):
    parse_functions = {
        "LULTRA": vin._parse_bool,
        "LUNSCR": vin._parse_bool,
        "LCOR": vin._parse_bool,
        "LPAW": vin._parse_bool,
        "EATOM": vin._parse_float,
        "RPACOR": vin._parse_float,
        "POMASS": vin._parse_float,
        "ZVAL": vin._parse_float,
        "RCORE": vin._parse_float,
        "RWIGS": vin._parse_float,
        "ENMAX": vin._parse_float,
        "ENMIN": vin._parse_float,
        "EMMIN": vin._parse_float,
        "EAUG": vin._parse_float,
        "DEXC": vin._parse_float,
        "RMAX": vin._parse_float,
        "RAUG": vin._parse_float,
        "RDEP": vin._parse_float,
        "RDEPT": vin._parse_float,
        "QCUT": vin._parse_float,
        "QGAM": vin._parse_float,
        "RCLOC": vin._parse_float,
        "IUNSCR": vin._parse_int,
        "ICORE": vin._parse_int,
        "NDATA": vin._parse_int,
        "VRHFIN": vin._parse_string,
        "LEXCH": vin._parse_string,
        "TITEL": vin._parse_string,
        "STEP": vin._parse_list,
        "RRKJ": vin._parse_list,
        "GGA": vin._parse_list
    }

    def __init__(self, data):
        search_lines = re.findall(r"(?s)(parameters from PSCTR are:"
                                  r".*?END of PSCTR-controll parameters)",
                                  data)  #.group(2)

        self.all_keywords = []
        for i in range(len(search_lines)):
            keywords = {}
            for key, val in re.findall(r"(\S+)\s*=\s*(.*?)(?=;|$)",
                                       search_lines[i],
                                       flags=re.MULTILINE):
                try:
                    keywords[key] = self.parse_functions[key](val)
                except KeyError:
                    warnings.warn("Ignoring unknown variable type %s" % key)
            self.all_keywords.append(keywords)

    def get_max_enmax(self):
        max_enmax = 0
        for i in range(len(self.all_keywords)):
            enmax = self.all_keywords[i]['ENMAX']
            if max_enmax < enmax:
                max_enmax = enmax
        return max_enmax


# p=PotcarSingle(file.read())
# print(p.enmax)
