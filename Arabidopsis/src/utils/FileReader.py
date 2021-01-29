from src.Block import AnnoBlock, SampleBlock
import pyBigWig
import re
import pandas as pd


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

                chrom_num = cells[0]

                information_list = cells[8].split(';')
                if len(information_list) == 3 and information_list[2][0] == 'N':
                    gene_number = AnnoReader.clear_return(AnnoReader.keep_right(information_list[2], '='))
                block = AnnoBlock(chrid= chrom_num,
                                  type=cells[2],
                                  start=int(cells[3]),
                                  end=int(cells[4])-1,
                                  strand=cells[6],
                                  parentgene=gene_number)
                annotation.append(block)
            f.close()
        return annotation

    @staticmethod
    def clear_return(strr: str):
        a = re.sub(r'\n', '', strr)
        return a

    @staticmethod
    def keep_right(string, symbol):
        split_data = string.split(symbol)
        return split_data[1]

    # below is specially for PlotGene
    @staticmethod
    def read_annotation_df(file):
        columnsname = ['chrid', 'type_gff', 'start', 'end',
                       'strand', 'parentgene', 'parentrna']
        n = len(columnsname)
        # Create several blank lists to store data
        lists = [[] for _ in range(n)]
        annotation_dict = {}
        gene_number = -1
        rna_number = -1
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cells = line.split('\t')
                chrom_num = cells[0]

                information_list = cells[8].split(';')
                if len(information_list) == 3 and information_list[2][0] == 'N':
                    gene_number = AnnoReader.clear_return(AnnoReader.keep_right(information_list[2], '='))
                if len(information_list) == 1 and information_list[0][0] == 'P':
                    rna_number = AnnoReader.clear_return(AnnoReader.keep_right(information_list[0], '='))
                lists[0].append(chrom_num)
                lists[1].append(cells[2])
                lists[2].append(int(cells[3]))
                lists[3].append(int(cells[4])-1)
                lists[4].append(cells[6])
                lists[5].append(gene_number)
                lists[6].append(rna_number)
            f.close()
            # Put data into dictionary
            for i in range(n):
                annotation_dict.update({columnsname[i]: lists[i]})
            annotation_df = pd.DataFrame(annotation_dict)
        return annotation_df


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
