import numpy as np
import math


class SizeFactorOneRegion:
    def __init__(self, gene_data: list, samples_condition1, percent: list, theta):
        """
        :param gene_data: list[i] = sample i. No position information
        :param samples_condition1:
        :param percent:
        :param theta:
        """
        self._gene_data = gene_data
        self._samples_condition1 = samples_condition1
        self._percent = percent
        self._theta = theta

    def count_size_factor_one_region(self):
        sorted_gene = []
        size_all = []
        for data in self._gene_data:
            sorted_data, size, start = SizeFactorOneRegion.size_calculator(data, self._percent, self._theta)
            sorted_gene.append(sorted_data)
            size_all.append(size)
        size_common = min(size_all)
        total_value = []
        for i in range(len(sorted_gene)):
            total = sum(sorted_gene[i][start:start + size_common])
            total_value.append(total)
        if min(total_value) == 0.0:
            size_factor = None
        else:
            mean_value = np.mean(total_value)
            size_factor = [va / mean_value for va in total_value]
        return size_factor

    @staticmethod
    def size_calculator(data: list, percent: list, theta):
        sorted_data = sorted(data, reverse=True)
        n = len(sorted_data)
        start = int(n * percent[0])
        end = int(n * percent[1])
        coveragebetween = sorted_data[start:end]
        pos = 0
        while pos < (len(coveragebetween) - 1):
            if coveragebetween[pos] > coveragebetween[pos + 1] + theta:
                break
            else:
                pos = pos + 1
        return sorted_data, pos, start

    def normalization(self, size):
        new_data = []
        for i in range(len(self._gene_data)):
            norm_data = list(map(lambda x:x/size[i], self._gene_data[i]))
            new_data.append(norm_data)
        return new_data

    @staticmethod
    def count_logratio(num1, num2):
        return math.log(num1+1) - math.log(num2+1)
