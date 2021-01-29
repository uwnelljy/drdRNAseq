import math


class Gene:
    def __init__(self):
        pass

    @staticmethod
    def get_whole_gene(chrom, start, end, samples):
        """
        :return: a list of whole gene at samples. list[i] = sample i, list[-1] = position
        """
        whole_gene = [sample_block._sample.values(chrom, start, end)
                      for sample_block in samples]
        return whole_gene

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
        ratio = [Gene.count_logratio(sum1[i] / num_condition1, sum2[i] / (m - num_condition1))
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

    # @staticmethod
    # def get_log_ratio(gene_data, num_condition1):
    #     """
    #     :param gene_data: list[i] = sample i, list[-1] is position.
    #     :param num_condition1:
    #     :return: log ratio contains position
    #     """
    #     n = len(gene_data[0])
    #     m = len(gene_data)-1
    #     sum1 = [0] * n
    #     sum2 = [0] * n
    #     for i in range(m):
    #         if i < num_condition1:
    #             sum1 = list(map(lambda x, y: x + y, sum1, gene_data[i]))
    #         else:
    #             sum2 = list(map(lambda x, y: x + y, sum2, gene_data[i]))
    #     ratio = [Gene.count_logratio(sum1[i] / num_condition1, sum2[i] / (m - num_condition1))
    #              for i in range(len(sum1))]
    #     ratio_with_position = [ratio, gene_data[-1]]
    #     return ratio_with_position

    @staticmethod
    def count_logratio(num1, num2):
        return math.log(num1+1)-math.log(num2+1)
