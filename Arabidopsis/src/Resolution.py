import numpy as np
from src.Normalization import SizeFactorOneRegion
import logging


class Resolution:
    """
    change resolution to reduce computational pressure.
    """
    def __init__(self, target_whole_gene: list, target_gene_position, width):
        """
        :param target_whole_gene: a list. list[i] = sample i
        :param width:
        """
        self._target_whole_gene = target_whole_gene
        self._target_gene_position = target_gene_position
        self._width = width

    def get_shorter_gene(self):

        """
        :return: a list of shorter gene at all samples. List[i] = sample i
        """
        shorter_index = Resolution.get_index(self._target_whole_gene[0], self._width)
        shorter_position = [self._target_gene_position[i] for i in shorter_index]
        shorter_gene = []
        n = len(self._target_whole_gene)
        for i in range(n):
            gene = self._target_whole_gene[i]
            new_gene = []
            for j in range(len(shorter_index)-1):
                new_gene.append(np.mean(gene[shorter_index[j]: shorter_index[j+1]]))
            new_gene.append(np.mean(gene[shorter_index[-1]:]))
            shorter_gene.append(new_gene)
        return shorter_gene, shorter_position

    @staticmethod
    def get_log_ratio(gene_data, num_condition1):
        """
        :param gene_data: list[i] = sample i
        :param num_condition1:
        :return: log ratio
        """
        n = len(gene_data[0])
        m = len(gene_data)
        sum1 = [0] * n
        sum2 = [0] * n
        for i in range(m):
            if i < num_condition1:
                sum1 = list(map(lambda x, y: x + y, sum1, gene_data[i]))
            else:
                sum2 = list(map(lambda x, y: x + y, sum2, gene_data[i]))
        ratio = [SizeFactorOneRegion.count_logratio(sum1[i] / num_condition1, sum2[i] / (m - num_condition1))
                 for i in range(len(sum1))]
        substracted = [sum1[i] / num_condition1 - sum2[i] / (m - num_condition1) for i in range(len(sum1))]
        max_sub = max(substracted)
        min_sub = min(substracted)
        if abs(max_sub) > abs(min_sub):
            dominator = abs(max_sub)
        else:
            dominator = abs(min_sub)
        ratio_weight = list(map(lambda x, y: abs(x) * y / dominator, ratio, substracted))
        return ratio, substracted, ratio_weight

    @staticmethod
    def cut_regions(data: list, width):
        data_array = np.array(data)
        position_index_zero = np.where(data_array == 0)[0]
        if len(position_index_zero) == 0:
            logging.debug('no need to cut, no zero exist')
            return None
        zero_start = [position_index_zero[0]]
        zero_end = []
        for i in range(len(position_index_zero) - 1):
            if position_index_zero[i] + 1 == position_index_zero[i + 1]:
                pass
            else:
                zero_end.append(position_index_zero[i])
                zero_start.append(position_index_zero[i + 1])
        zero_end.append(position_index_zero[-1])
        if len(zero_start) != len(zero_end):
            raise RuntimeError('zero start and zero end should have same length')
        cut_start = []
        cut_end = []
        for i in range(len(zero_start)):
            zero_width = zero_end[i] - zero_start[i]
            if zero_width >= width:
                cut_start.append(zero_start[i])
                cut_end.append(zero_end[i])
        if len(cut_start) == 0 or len(cut_end) == 0:
            logging.debug('no need to cut, all less than {}'.format(width))
            return None
        return cut_start, cut_end

    @staticmethod
    def get_index(data: list, width: int):
        """
        :param data: gene list data
        :param width: resolution
        :return: shorter gene
        """
        n = len(data)
        m = n - width
        num = int(m / width)
        temp = np.linspace(0, m, num)
        return [int(va) for va in temp]
