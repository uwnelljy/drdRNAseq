from src.Block import AnnoBlock, SampleBlock
import re
import pyBigWig


class AnnoReader:
    def __init__(self):
        pass

    @staticmethod
    def read_annotation(file):
        annotation = []
        gene_number = -1
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cells = line.split('\t')
                ischr, chrid1, chrid2 = AnnoReader.extract_chromosome(cells[0])
                id_block = AnnoReader.extract_parentandid(cells[8])
                if id_block[0] == 'g':
                    gene_number = int(re.sub('\D', '', id_block))
                    parent_gene = gene_number
                else:
                    parent_gene = gene_number
                block = AnnoBlock(ischr=ischr,
                                  chrid1=chrid1,
                                  chrid2=chrid2,
                                  type=cells[2],
                                  start=int(cells[3]),
                                  end=int(cells[4])-1,
                                  strand=cells[6],
                                  id=id_block,
                                  parentgene=parent_gene)
                annotation.append(block)
            f.close()
        return annotation

    @staticmethod
    def extract_chromosome(string):
        split_list = string.split('_')
        chrid = split_list[1].split('.')
        return split_list[0], int(chrid[0]), int(chrid[1])

    @staticmethod
    def extract_parentandid(string):
        split_list = string.split(';')
        idid = AnnoReader.keep_right(split_list[0], '=')
        return idid

    @staticmethod
    def keep_right(string, symbol):
        split_data = string.split(symbol)
        return split_data[1]


class GeneReader:
    def __init__(self):
        pass

    @staticmethod
    def get_bw(file_list: list):
        """
        :param file_list: List of file names
        :return: Samples[i] = bw data of sample i
        """
        samples = []
        for va in file_list:
            samples.append(SampleBlock(pyBigWig.open(va)))
        return samples
