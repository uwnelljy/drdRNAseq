from src.Gene import Gene


class Overlap:
    def __init__(self):
        pass

    @staticmethod
    def overlap_gene(parentgene, annotation, strand):
        if strand == '+':
            former_gene_id = parentgene
            latter_gene_id = parentgene + 1
        else:
            former_gene_id = parentgene - 1
            latter_gene_id = parentgene

        # gene_info = (gene_start, gene_end, gene_strand)
        former_info = Gene.gene_info_for_overlap(former_gene_id, annotation)
        latter_info = Gene.gene_info_for_overlap(latter_gene_id, annotation)

        return former_info, latter_info

    @staticmethod
    def check_overlap(former_info, latter_info):
        former_end = former_info[1]
        latter_start = latter_info[0]
        if latter_start < former_end:
            return True  # overlap
        else:
            return False # not overlap

    # # When overlap: check_overlap is True.
    # @staticmethod
    # def get_common_region_not_overlapped(former_info, latter_info,
    #                                      former_last_common_region,
    #                                      latter_first_common_region):
    #     common_region_not_overlapped_former = [former_last_common_region[0], latter_info[0]]
    #     common_region_not_overlapped_latter = [former_info[1], latter_first_common_region[1]]
    #     return common_region_not_overlapped_former, common_region_not_overlapped_latter

    # @staticmethod
    # def check_last_common_region_in_overlap(common_region):
    #     # Last common region might in overlapping area.
    #     # If former last common region is in overlapping area, then check latter common region.
    #     # If they both in overlapping area, then unsolvable.
    #     if common_region[0] > common_region[1]:
    #         return True # last common region is in overlapping area
    #     return False
