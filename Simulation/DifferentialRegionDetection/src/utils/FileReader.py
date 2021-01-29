from src.Block import SampleBlock, Block
import pyBigWig


class BedReader:
    def __init__(self):
        pass

    # @staticmethod
    # def read_annotation(file):
    #     annotation = []
    #     with open(file, 'r') as f:
    #         for line in f:
    #             if line.startswith('#'):
    #                 continue
    #             cells = line.split('\t')
    #             start = int(BedReader.clear_return(cells[11])) + int(cells[1])
    #             block = AnnoBlock(chr=cells[0],
    #                               type=cells[3],
    #                               genestart=int(cells[1]),
    #                               geneend=int(cells[2]),
    #                               isostart=int(start),
    #                               isoend=int(start)+int(cells[10]),
    #                               strand=cells[5])
    #             annotation.append(block)
    #         f.close()
    #     return annotation
    #
    # @staticmethod
    # def clear_return(strr: str):
    #     a = re.sub(r'\n', '', strr)
    #     return a


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


class ExprReader:
    def __init__(self, filename):
        self.filename = filename

    def filereader(self):
        isoblock = []
        with open(self.filename, 'r') as f:
            for line in f:
                l = line.split('\t')
                isoblock.append(Block(iso=l[0], value=int(l[1])))
            f.close()
        return isoblock
