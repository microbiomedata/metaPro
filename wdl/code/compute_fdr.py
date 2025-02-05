import argparse
import csv
from decimal import *
from collections import namedtuple
from functools import partial
import json
import sys

FdrResult = namedtuple('FdrResult', 'SpecEValue FDR')

class FdrSearchStrategy():
    def __init__(self, data):
        self.data = data
        self.decoys = list(filter(lambda data_row: data_row['Protein'].startswith('XXX'), self.data))
        super().__init__()

    def __select(self, data_row, spec_e_value: Decimal):
        return data_row['MSGFDB_SpecEValue'] <= spec_e_value

    def __count_filtered(self, data_list, spec_e_value: Decimal):
        return len(list(filter(partial(self.__select, spec_e_value = spec_e_value), data_list)))

    def fdr1(self, spec_e_value: Decimal) -> Decimal:
       fdr = ((Decimal(self.__count_filtered(self.decoys, spec_e_value)) * 2) / self.__count_filtered(self.data, spec_e_value))
       return fdr

    def find_values(self) -> FdrResult:
        pass

class FdrImmediate(FdrSearchStrategy):
    def __init__(self, data, spec_e_value: Decimal):
        self.spec_e_value = spec_e_value
        super().__init__(data)

    def find_values(self) -> FdrResult:
        fdr = self.fdr1(self.spec_e_value)

        return FdrResult(self.spec_e_value, fdr)
              
class FdrDownIterate(FdrSearchStrategy):
    def __init__(self, data, spec_e_value: Decimal, fdr: Decimal, spec_e_inc: Decimal):
        self.spec_e_value = spec_e_value
        self.fdr = fdr
        self.spec_e_inc = spec_e_inc
        super().__init__(data)

    def find_values(self) -> FdrResult:
        fdr_tmp = self.fdr
        spec_e_tmp = self.spec_e_value
 
        spec_e_last = fdr_tmp
        fdr_last = spec_e_tmp

        while True:
            spec_e_tmp -= self.spec_e_inc 
            fdr_tmp = self.fdr1(spec_e_tmp)

            if self.fdr > fdr_tmp:
                break

            spec_e_last = spec_e_tmp
            fdr_last = fdr_tmp

        return FdrResult(spec_e_last, fdr_last)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, help='First Hits file to process')
    parser.add_argument('--fdr', type=str, help='')
    parser.add_argument('--spececoeff', type=str, help='')
    parser.add_argument('--speceexp', type=str, help='')
    parser.add_argument('--spece', type=str, help='')
    parser.add_argument('--speceinc', type=str, help='')
    parser.add_argument('--precision', type=int, help='Floating point precision', required=False)
    parser.add_argument('--single', action='store_true', help='Prints FDR for given SpecEValue')
    # parser.add_argument('--log', action='store_true')
    parser.add_argument('--out', type=str, help='', default='out.json')

    return parser.parse_args()

def get_dat(filepath: str):
    with open(filepath) as data_file:
        reader = csv.DictReader(data_file, delimiter='\t')

        return [{'Protein': row['Protein'], 'MSGFDB_SpecEValue': Decimal(row['MSGFDB_SpecEValue'])} for row in reader]

if __name__ == "__main__":
    args = get_args()
    
    getcontext().prec = args.precision if args.precision is not None else 45

    rows = get_dat(args.file)
    result = None
    to_write = None

    spec_e_value = Decimal(args.spece)
    spec_e_value_inc = Decimal(args.speceinc)
    fdr = Decimal(args.fdr)

    fdr_iterate = FdrDownIterate(rows, spec_e_value, fdr, spec_e_value_inc)
    result = fdr_iterate.find_values()
    to_write = result.SpecEValue

    sys.stdout.write(str(to_write))

    with open(args.out, 'w+') as fp:  
        out_dat = {
            "FDR": float(result.FDR),
            "SpecEValue": float(result.SpecEValue),
        }
        dat = json.dumps(out_dat)     
        fp.write(dat)