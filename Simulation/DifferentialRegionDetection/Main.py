from src.Gene import Gene
from src.utils.FileReader import GeneReader
from src.Detection import EventDetection
from src.Normalization import SizeFactorOneRegion
import numpy as np
import logging
import time
import matplotlib.pyplot as plt
import traceback

logging.basicConfig(filename='./syslogging.log',
                    format='[%(asctime)s-%(filename)s-%(levelname)s:%(message)s]',
                    filemode='a',
                    level=logging.DEBUG,
                    datefmt='%Y-%m-%d%I:%M:%S %p')

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.figsize'] = (16, 20)


PATH = './data/'
FILE = PATH + 'ref.bed12'
SAMPLES_CONDITION1 = 1
PERCENT = [0.1, 0.25]
THETA = np.tan(89/180*np.pi)
THRESHOLD = 1e-3
THRES = 1.5
LAMBDA2_LIST = np.logspace(1, 3, 100)
CHROM = 'chr1'
GENE_START = 10000
GENE_END = 13000
POSITION = range(GENE_START, GENE_END)
# TOTAL_NUMBER = list(range(1, 501))
TOTAL_NUMBER = [2]


class Main:
    def __init__(self, genenumber, samples_condition1, threshold):
        self._genenumber = genenumber
        self._samples_condition1 = samples_condition1
        self._threshold = threshold

    def get_samples(self):
        return GeneReader.get_bw([PATH + 'condition{}_1.bw'.format(self._genenumber),
                                  PATH + 'condition{}_2.bw'.format(self._genenumber)])

    def get_whole_gene_before_norm(self):
        samples = self.get_samples()
        return Gene.get_whole_gene(CHROM, GENE_START, GENE_END, samples)

    def get_size_factor(self):
        whole_gene_before_norm = self.get_whole_gene_before_norm()
        to_count_size = SizeFactorOneRegion(whole_gene_before_norm,
                                            SAMPLES_CONDITION1, PERCENT, THETA)
        size_factor = to_count_size.count_size_factor_one_region()
        return size_factor

    def get_whole_gene_after_norm(self):
        whole_gene_before_norm = self.get_whole_gene_before_norm()
        to_count_size = SizeFactorOneRegion(whole_gene_before_norm,
                                            SAMPLES_CONDITION1, PERCENT, THETA)
        size_factor = self.get_size_factor()
        whole_gene_after_norm = to_count_size.normalization(whole_gene_before_norm, size_factor)
        return whole_gene_before_norm, whole_gene_after_norm

    def get_ratio(self):
        whole_gene_after_norm = self.get_whole_gene_after_norm()[1]
        whole_gene_before_norm = self.get_whole_gene_after_norm()[0]
        ratio, subtraction, weighted_ratio = Gene.get_log_ratio(whole_gene_after_norm, self._samples_condition1)
        return whole_gene_before_norm, whole_gene_after_norm, ratio, subtraction, weighted_ratio

    @staticmethod
    def get_target_gene_needed(whole_gene_before_norm, whole_gene_after_norm,
                               log_ratio, substraction, weight_ratio, genenumber):
        # plot common gene before normalization
        plt.subplot(511)
        for i in range(len(whole_gene_before_norm)):
            plt.scatter(range(len(whole_gene_before_norm[i])), whole_gene_before_norm[i], label='rna sample{}'.format(i))
        plt.title('before normalization')
        plt.legend(loc='upper left')
        # plot common gene after normalization
        plt.subplot(512)
        for i in range(len(whole_gene_after_norm)):
            plt.scatter(range(len(whole_gene_after_norm[i])), whole_gene_after_norm[i], label='rna sample{}'.format(i))
        plt.title('after normalization')
        plt.legend(loc='upper left')
        # plot substraction of gene data under two conditions
        plt.subplot(513)
        plt.scatter(range(len(substraction)), substraction, color='grey', label='substraction')
        plt.title('subtraction')
        plt.legend(loc='upper left')
        # plot log_ratio of gene data under two conditions
        plt.subplot(514)
        plt.scatter(range(len(log_ratio)), log_ratio, color='grey', label='log_ratio')
        plt.title('log_ratio')
        plt.legend(loc='upper left')
        # plot weight_ratio(substraction is the weight) under two conditions.
        plt.subplot(515)
        plt.scatter(range(len(weight_ratio)), weight_ratio, color='grey', label='weight_ratio')
        plt.title('weight_ratio')
        plt.legend(loc='upper left')

        plt.savefig('./ratiopic/gene{}'.format(genenumber))
        # plt.show()
        plt.clf()

        if abs(max(substraction)) <= 30 and abs(min(substraction)) <= 30:
            logging.debug('substraction is too small.')
            return None

    @staticmethod
    def run_for_one_lambda2(weight_ratio_data, position, genenumber, lambda2, for_plot, whole_gene_after_norm):
        detection = EventDetection(weight_ratio_data, SAMPLES_CONDITION1, lambda2)
        beta_data, cp_position_index = detection.changepoint_position()
        sse = sum(list(map(lambda x, y: (x - y) ** 2, weight_ratio_data, beta_data)))
        if cp_position_index is None:
            return None
        cp_pos = [position[ii] for ii in cp_position_index]
        if for_plot is None:
            return cp_pos, sse
        else:
            # plot whole region weight ratio data
            plt.subplot(211)
            plt.scatter(range(len(weight_ratio_data)), weight_ratio_data, color='grey', label='ratio')
            plt.plot(range(len(beta_data)), beta_data, color='orange', label='beta')
            Main.plot_region(cp_position_index)
            plt.legend(loc='upper left')

            # plot position on whole gene
            plt.subplot(212)
            for i in range(len(whole_gene_after_norm)):
                plt.scatter(position, whole_gene_after_norm[i],
                            label='rna sample{}'.format(i))
            Main.plot_region(cp_pos)
            plt.legend(loc='upper left')

            plt.savefig('./cpresult/gene{}'.format(genenumber))
            # plt.show()
            plt.clf()
            return cp_pos

    @staticmethod
    def plot_region(cp_position):
        for pos in cp_position:
            if cp_position.index(pos) == 0:
                plt.axvline(pos, color='grey', alpha=0.5, label='cp')
            else:
                plt.axvline(pos, color='grey', alpha=0.5)


class KnownError(Exception):
    pass


if __name__ == '__main__':
    position = []
    time_each = []
    for genenumber in TOTAL_NUMBER:
        logging.debug('=======================================================================')
        logging.debug('\t>><<>><<>><< Executing gene: {} <<>><<>><<>>'.format(genenumber))
        start_time = time.time()
        try:
            main_func = Main(genenumber, SAMPLES_CONDITION1, THRESHOLD)
            whole_gene_before_norm = main_func.get_ratio()[0]
            whole_gene_after_norm = main_func.get_ratio()[1]
            log_ratio = main_func.get_ratio()[2]
            substraction = main_func.get_ratio()[3]
            weight_ratio = main_func.get_ratio()[4]
            Main.get_target_gene_needed(whole_gene_before_norm, whole_gene_after_norm,
                                        log_ratio, substraction, weight_ratio, genenumber)
            kk = []
            lam = []
            jj = []
            initial = 9999
            for lamb2 in LAMBDA2_LIST:
                result_cp = Main.run_for_one_lambda2(weight_ratio, POSITION, genenumber, lamb2, None,
                                                     whole_gene_after_norm)
                if result_cp is not None:
                    cp_posi = result_cp[0]
                    ssee = result_cp[1]
                    if len(cp_posi) == initial:
                        pass
                    else:
                        initial = len(cp_posi)
                        lam.append(lamb2)
                        kk.append(len(cp_posi)+1)
                        jj.append(ssee)
                else:
                    lam.append(lamb2)
                    kk.append(1)
                    jj.append(ssee)
                    break
            result_ratio_dk = EventDetection.get_ratio_and_dk(kk, lam, jj)
            if result_ratio_dk is None:
                position.append(None)
                logging.debug('change point position is None')
            else:
                ratioo = result_ratio_dk[0]
                subb = result_ratio_dk[1]
                dk = result_ratio_dk[2]
                lambda2_final, cp_num = EventDetection.determine_lambda2_and_cp_num(kk, lam, ratioo, subb, dk, THRES)
                logging.debug('lambda2 is {} and the number of change point is {}'.format(lambda2_final, cp_num))
                cp_position = Main.run_for_one_lambda2(weight_ratio, POSITION, genenumber, lambda2_final, 'yes',
                                                       whole_gene_after_norm)
                position.append(cp_position)
                logging.debug('raw change point position is {}'.format(cp_position))
        except KnownError:
            logging.warning('known error encountered at this gene')
            logging.debug('change point position is None')
            position.append(None)
        except Exception:
            logging.critical('other error encountered at this gene')
            logging.error(traceback.format_exc(limit=None))
            logging.debug('change point position is None')
            position.append(None)
        finally:
            pass
        logging.debug('********End********')
        end_time = time.time()
        logging.debug('time used:{}s'.format(end_time - start_time))
        time_each.append(end_time-start_time)
    logging.debug('position list is {}'.format(position))
    logging.debug('time list is {}'.format(time_each))
