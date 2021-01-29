from src.utils.FileReader import AnnoReader, GeneReader
from src.Gene import Gene
from src.Resolution import Resolution
from src.Normalization import SizeFactorOneRegion
from src.Detection import EventDetection
import numpy as np
import logging
import traceback
import time
import matplotlib.pyplot as plt

logging.basicConfig(filename='./syslogging.log',
                    format='[%(asctime)s-%(filename)s-%(levelname)s:%(message)s]',
                    filemode='a',
                    level=logging.DEBUG,
                    datefmt='%Y-%m-%d%I:%M:%S %p')

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.figsize'] = (16, 20)


PATH = './data/'
FILE = PATH + 'TAIR10_GFF3_genes.gff'
PLUS = [PATH + 'GSM1603462_WT_plus_strand.bw', PATH + 'GSM1603463_WT_plus_strand.bw',
        PATH + 'GSM1603464_WT_plus_strand.bw', PATH + 'GSM1603465_M_plus_strand.bw',
        PATH + 'GSM1603466_M_plus_strand.bw', PATH + 'GSM1603467_M_plus_strand.bw']
MINUS = [PATH + 'GSM1603462_WT_minus_strand.bw', PATH + 'GSM1603463_WT_minus_strand.bw',
         PATH + 'GSM1603464_WT_minus_strand.bw', PATH + 'GSM1603465_M_minus_strand.bw',
         PATH + 'GSM1603466_M_minus_strand.bw', PATH + 'GSM1603467_M_minus_strand.bw']
SAMPLES_CONDITION1 = 3
PERCENT = [0.1, 0.25]
THETA = np.tan(89/180*np.pi)
THRESHOLD = 1e-3
WIDTH1 = 10
WIDTH2 = 20
WIDTH3 = 100
WIDTH4 = 200
WIDTH5 = 500
LAMBDA2_LIST = np.logspace(1, 3, 100)
THRES = 1.5


class Main:
    def __init__(self, parentgene, samples_condition1, width1, width2, width3, width4, width5):
        self._parentgene = parentgene
        self._samples_condition1 = samples_condition1
        self._width1 = width1
        self._width2 = width2
        self._width3 = width3
        self._width4 = width4
        self._width5 = width5

    @staticmethod
    def get_annotation():
        return AnnoReader.read_annotation(FILE)

    @staticmethod
    def get_samples(file: list):
        return GeneReader.get_bw(file)

    @staticmethod
    def plot_region(cp_position):
        for pos in cp_position:
            if cp_position.index(pos) == 0:
                plt.axvline(pos, color='grey', alpha=0.5, label='change point')
            else:
                plt.axvline(pos, color='grey', alpha=0.5)

    @staticmethod
    def get_parentgene_strand(annotation):
        parentgene_list = []
        strand_list = []
        for va in annotation:
            if va._type == 'gene':
                parentgene_list.append(va._parentgene)
                strand_list.append(va._strand)
        return parentgene_list, strand_list

    def get_target_gene_needed(self, samples, annotation):
        target_annotation = Gene.get_target_annotation(annotation, self._parentgene)
        gene_position = Gene.get_target_chr_start_end(target_annotation)
        if gene_position is None:
            logging.debug('End the process. There is no eligible position information for that gene.')
            return None
        gene_chrom = gene_position[0]
        gene_start = gene_position[1]
        gene_end = gene_position[2]

        # get common region
        common_region = Gene.get_common_region(target_annotation)
        if gene_end-gene_start < 5000 or common_region is None:
            logging.debug('No applicable mRNA or check intron. Use whole gene.')
            whole_gene = Gene.get_whole_gene(gene_chrom, gene_start, gene_end, samples)
            common_position = list(range(gene_start, gene_end))
        else:
            whole_gene, common_position = Gene.get_region_gene(gene_chrom, common_region, samples)
            logging.debug(
                'the original length is {}. common length is {}'.format(gene_end - gene_start, len(common_position)))
        if Gene.check_value(whole_gene) is False:
            logging.debug('End the process. Count value is zero.')
            return None
        # normalization
        logging.debug('********Start Normalization********')
        to_count_size = SizeFactorOneRegion(whole_gene, self._samples_condition1, PERCENT, THETA)
        size_factor = to_count_size.count_size_factor_one_region()
        if size_factor is None:
            logging.debug('End the process. Size factor is None')
            return None
        logging.debug('Size factor is {}'.format(size_factor))
        whole_gene_normalized = to_count_size.normalization(size_factor)

        # plot common gene before normalization
        # plt.figure(figsize=(16,28))
        plt.subplot(511)
        for i in range(len(whole_gene)):
            plt.scatter(range(len(whole_gene[i])), whole_gene[i], label='rna sample{}'.format(i))
        plt.title('before normalization')
        plt.legend(loc='upper left')
        # plot common gene after normalization
        plt.subplot(512)
        for i in range(len(whole_gene_normalized)):
            plt.scatter(range(len(whole_gene_normalized[i])), whole_gene_normalized[i], label='rna sample{}'.format(i))
        plt.title('after normalization')
        plt.legend(loc='upper left')

        logging.debug('********Start Changing Resolution********')
        length_aft_cut = len(whole_gene_normalized[0])
        if length_aft_cut >= 5000:
            logging.debug('The length after cut is {}, change resolution.'.format(length_aft_cut))
            if length_aft_cut < 50000:
                width = self._width1
            elif 50000 <= length_aft_cut < 100000:
                width = self._width2
            elif 100000 <= length_aft_cut < 500000:
                width = self._width3
            elif 500000 <= length_aft_cut < 1000000:
                width = self._width4
            else:
                width = self._width5
            logging.debug('width of resolution is {}'.format(width))
            to_change_resolution = Resolution(whole_gene_normalized, common_position, width)
            shorter_gene, shorter_position = to_change_resolution.get_shorter_gene()
            if len(shorter_gene[0]) == len(shorter_position):
                logging.debug('The length of shaped gene is {}'.format(len(shorter_gene[0])))
            else:
                raise RuntimeError('End the process. Shorter gene and position have different length')
        else:
            logging.debug('The length after cut is {}, no need to change resolution'.format(length_aft_cut))
            shorter_gene = whole_gene_normalized
            shorter_position = common_position
        log_ratio, substraction, weight_ratio = Resolution.get_log_ratio(shorter_gene, self._samples_condition1)

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

        plt.savefig('./ratiopic/gene{}'.format(self._parentgene))
        # plt.show()
        plt.clf()

        if abs(max(substraction)) <= 30 and abs(min(substraction)) <= 30:
            logging.debug('substraction is too small.')
            return None

        return whole_gene_normalized, shorter_position, weight_ratio, common_position

    @staticmethod
    def run_for_one_lambda2(log_ratio_data, position_data, parentgene, lambda2, for_plot, common_gene_position,
                            whole_gene_normalized):
        detection = EventDetection(log_ratio_data, SAMPLES_CONDITION1, lambda2)
        beta_data, cp_position_index = detection.changepoint_position()
        sse = sum(list(map(lambda x, y: (x - y) ** 2, log_ratio_data, beta_data)))
        if cp_position_index is None:
            return None
        cp_pos = [position_data[ii] for ii in cp_position_index]
        if for_plot is None:
            return cp_pos, sse
        else:
            # get position in common_gene_position
            common_gene_index = []
            for position in cp_pos:
                common_gene_index.append(common_gene_position.index(position))

            # plot common region weight ratio data
            plt.subplot(211)
            plt.scatter(range(len(log_ratio_data)), log_ratio_data, color='grey', label='ratio')
            plt.plot(range(len(beta_data)), beta_data, color='orange', label='beta')
            Main.plot_region(cp_position_index)
            plt.legend(loc='upper left')

            # plot position on common gene
            plt.subplot(212)
            for i in range(len(whole_gene_normalized)):
                plt.scatter(range(len(whole_gene_normalized[i])), whole_gene_normalized[i],
                            label='rna sample{}'.format(i))
            Main.plot_region(common_gene_index)
            plt.legend(loc='upper left')

            plt.savefig('./cpresult/gene{}'.format(parentgene))
            # plt.show()
            plt.clf()
            return cp_pos


class KnownError(Exception):
    pass


if __name__ == '__main__':
    annotation = Main.get_annotation()
    samples_plus = Main.get_samples(PLUS)
    samples_minus = Main.get_samples(MINUS)
    parentgene_list1, strand_list = Main.get_parentgene_strand(annotation)
    parentgene_list = ['AT5G64120']
    for parentgene in parentgene_list:
        logging.debug('=======================================================================')
        logging.debug('\t>><<>><<>><< Executing gene: {} <<>><<>><<>>'.format(parentgene))
        start_time = time.time()
        parentgene_index = parentgene_list1.index(parentgene)
        parentgene_strand = strand_list[parentgene_index]
        logging.debug('strand is {}'.format(parentgene_strand))
        if parentgene_strand == '+':
            samples = samples_plus
        else:
            samples = samples_minus
        try:
            main_func = Main(parentgene, SAMPLES_CONDITION1, WIDTH1, WIDTH2, WIDTH3,WIDTH4, WIDTH5)
            check_first = main_func.get_target_gene_needed(samples, annotation)
            if check_first is None:
                logging.debug('Skip this gene.')
                cp_position = None
            else:
                logging.debug('********Start Detection********')
                whole_gene_normalized = check_first[0]
                target_shorter_gene_position_list = check_first[1]
                log_ratio_data = check_first[2]
                common_gene_position = check_first[3]
                kk = []
                lam = []
                jj = []
                initial = 9999
                for lamb2 in LAMBDA2_LIST:
                    result_cp = Main.run_for_one_lambda2(log_ratio_data, target_shorter_gene_position_list,
                                                         parentgene, lamb2, None, common_gene_position,
                                                         whole_gene_normalized)
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
                    cp_position = None
                else:
                    ratioo = result_ratio_dk[0]
                    subb = result_ratio_dk[1]
                    dk = result_ratio_dk[2]
                    lambda2_final, cp_num = EventDetection.determine_lambda2_and_cp_num(kk, lam, ratioo, subb, dk, THRES)
                    logging.debug('lambda2 is {} and the number of change point is {}'.format(lambda2_final, cp_num))
                    cp_position = Main.run_for_one_lambda2(log_ratio_data, target_shorter_gene_position_list, parentgene,
                                                           lambda2_final, 'yes', common_gene_position,
                                                           whole_gene_normalized)
            logging.debug('raw change point position is {}'.format(cp_position))
        except KnownError:
            logging.warning('known error encountered at this gene')
            logging.debug('change point position is None')
        except Exception:
            logging.critical('other error encountered at this gene')
            logging.error(traceback.format_exc(limit=None))
            logging.debug('change point position is None')
        finally:
            pass
        logging.debug('********End********')
        end_time = time.time()
        logging.debug('time used:{}s'.format(end_time - start_time))
