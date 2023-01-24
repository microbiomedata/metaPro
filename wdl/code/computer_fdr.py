import argparse
import csv
from decimal import *
from typing import List, TypedDict
from timeit import default_timer
from collections import namedtuple
from functools import partial
import json

FdrResult = namedtuple('FdrResult', 'SpecEValue FDR')

class DataRow(TypedDict):
    Protein: str
    MSGFDB_SpecEValue: Decimal

class FdrSearchStrategy():
    def __init__(self, data: List[DataRow]):
        self.data = data
        self.decoys = list(filter(lambda data_row: data_row['Protein'].startswith('XXX'), self.data))
        super().__init__()

    def __select(self, data_row: DataRow, spec_e_value: Decimal):
        return data_row['MSGFDB_SpecEValue'] <= spec_e_value

    def __count_filtered(self, data_list: List[DataRow], spec_e_value: Decimal):
        return len(list(filter(partial(self.__select, spec_e_value = spec_e_value), data_list)))

    def fdr1(self, spec_e_value: Decimal) -> Decimal:
       fdr = ((Decimal(self.__count_filtered(self.decoys, spec_e_value)) * 2) / self.__count_filtered(self.data, spec_e_value)) * 100
       return fdr

    def find_values(self) -> FdrResult:
        pass

class FdrImmediate(FdrSearchStrategy):
    def __init__(self, data: List[DataRow], spec_e_value: Decimal):
        self.spec_e_value = spec_e_value
        super().__init__(data)

    def find_values(self) -> FdrResult:
        fdr = self.fdr1(self.spec_e_value)

        return FdrResult(self.spec_e_value, fdr)


class FdrLog(FdrSearchStrategy):
    def __init__(self, data: List[DataRow], target_fdr: Decimal, spec_e_value_coeff: Decimal, spec_e_value_exp: Decimal):
        self.spec_e_value_coeff = spec_e_value_coeff
        self.spec_e_value_exp = spec_e_value_exp
        self.fdr_target = target_fdr
        super().__init__(data)
        
    class FdrCalc():
        def __init__(self, current_compare_fn, fdr_fn, target_fdr, initial_fdr, coeff, exp, exp_inc, next_compare_fn, stop_condition_fn, stop_exp, is_done):
            self.current_compare_fn = current_compare_fn
            self.coeff = coeff
            self.fdr_fn = fdr_fn
            self.target_fdr = target_fdr
            self.exp = exp
            self.exp_inc = exp_inc
            self.next_compare_fn = next_compare_fn
            self.stop_condition_fn = stop_condition_fn
            self._is_done = is_done
            self.fdr_i = initial_fdr
            self.stop_exp = stop_exp
            self.ten = Decimal('10.0')
            self.n = Decimal('-1.0')
            self.inc_mult = Decimal('0.1')

        def do(self):
            while not self._is_done:
                if self.stop_condition_fn(abs(self.exp_inc), self.stop_exp):
                    self._is_done = True
                    return self
                
                self.exp += self.exp_inc
                self.fdr_i = self.fdr_fn(self.spec_e_value)
                if not self.current_compare_fn(self.target_fdr, self.fdr_i):
                    next_exp = self.n * self.exp_inc * self.inc_mult
                    return FdrLog.FdrCalc(
                        self.next_compare_fn,
                        self.fdr_fn,
                        self.target_fdr,
                        self.fdr_i,
                        self.coeff,
                        self.exp,
                        next_exp,
                        self.current_compare_fn,
                        self.stop_condition_fn,
                        self.stop_exp,
                        self.is_done
                    )


        @property
        def is_done(self):
            return self._is_done

        @property
        def fdr(self):
            return self.fdr_i

        @property
        def spec_e_value(self):
            return self.coeff * (self.ten ** self.exp)

    def find_values(self) -> FdrResult:
        less_than = lambda x, y : x < y
        greater_than = lambda x, y: x > y
        calc = FdrLog.FdrCalc(
            less_than,
            self.fdr1,
            self.fdr_target,
            Decimal('0.0'),
            self.spec_e_value_coeff,
            self.spec_e_value_exp,
            Decimal('-1.0'),
            greater_than,
            less_than,
            Decimal('0.00000001'),
            False
        )
        while not calc.is_done:
            calc = calc.do()

        return FdrResult(calc.spec_e_value, calc.fdr)
        
        
class FdrDownIterate(FdrSearchStrategy):
    def __init__(self, data: List[DataRow], spec_e_value: Decimal, fdr: Decimal, spec_e_inc: Decimal):
        self.spec_e_value = spec_e_value
        self.fdr = fdr
        self.spec_e_inc = spec_e_inc
        super().__init__(data)

    def find_values(self) -> FdrResult:
        fdr_tmp = self.fdr
        spec_e_tmp = self.spec_e_value
 
        spec_e_cur = None
        fdr_cur = None

        while True:
            spec_e_tmp -= self.spec_e_inc 
            fdr_tmp = self.fdr1(spec_e_tmp)

            print(f"{spec_e_tmp} {fdr_tmp}")

            if self.fdr > fdr_tmp:
                break

            spec_e_cur = spec_e_tmp
            fdr_cur = fdr_tmp

        return FdrResult(spec_e_cur, fdr_cur)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, help='First Hits file to process')
    parser.add_argument('--fdr', type=str, help='')
    parser.add_argument('--spec_e_coeff', type=str, help='')
    parser.add_argument('--spec_e_exp', type=str, help='')
    parser.add_argument('--spec_e_inc_coeff', type=str, help='')
    parser.add_argument('--spec_e_inc_exp', type=str, help='')
    parser.add_argument('--precision', type=int, help='Floating point precision', required=False)
    parser.add_argument('--iter', action='store_true')
    parser.add_argument('--single', action='store_true')
    # parser.add_argument('--log', action='store_true')
    parser.add_argument('--out', type=str, help='', default='out.json')

    return parser.parse_args()

def get_dat(filepath: str) -> List[DataRow]:
    with open(filepath) as data_file:
        reader = csv.DictReader(data_file, delimiter='\t')

        return [{'Protein': row['Protein'], 'MSGFDB_SpecEValue': Decimal(row['MSGFDB_SpecEValue'])} for row in reader]

if __name__ == "__main__":
    args = get_args()
    
    getcontext().prec = args.precision if args.precision is not None else 45

    rows = get_dat(args.file)
    result = None
    
    t_start = default_timer()
    if args.single:
        spec_e_value = Decimal(args.spec_e_coeff) * (Decimal(10) ** Decimal(args.spec_e_exp))

        fdr_imm = FdrImmediate(rows, spec_e_value)
        result = fdr_imm.find_values()
    elif args.iter:
        spec_e_value = Decimal(args.spec_e_coeff) * (Decimal(10) ** Decimal(args.spec_e_exp))
        spec_e_value_inc = Decimal(args.spec_e_inc_coeff) * (Decimal(10) ** Decimal(args.spec_e_inc_exp))
        fdr = Decimal(args.fdr)

        fdr_iterate = FdrDownIterate(rows, spec_e_value, fdr, spec_e_value_inc)
        result = fdr_iterate.find_values()
    # elif args.log:
    #     spec_e_value = Decimal(args.spec_e_coeff) * (Decimal(10) ** Decimal(args.spec_e_exp))
    #     fdr = Decimal(args.fdr)

    #     fdr_log = FdrLog(rows, fdr, Decimal(args.spec_e_coeff), Decimal(args.spec_e_exp))
    #     result = fdr_log.find_values()
    t_end = default_timer()

    print(f'Time to process: {t_end - t_start}s')

    with open(args.out, 'w+') as fp:       
        out_dat = {
            "FDR": float(result.FDR),
            "SpecEValue": float(result.SpecEValue),
        }
        json.dump(out_dat, fp)