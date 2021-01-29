class Gene:
    def __init__(self):
        pass

    @staticmethod
    def get_target_annotation(annotation, parentgene):
        target_annotation = []
        for va in annotation:
            if va._ischr == 'NC' and va._parentgene == parentgene:
                target_annotation.append(va)
        return target_annotation

    @staticmethod
    def get_target_chr_start_end(target_annotation):
        """
        :return: Annotation for target parent gene. It is a list. List[i] = row i
        """
        for va in target_annotation:
            if va._type == 'gene':
                return va._chrid1, va._start, va._end, va._strand
        return None

    @staticmethod
    def gene_info_for_overlap(geneid, annotation):
        target_annotation = Gene.get_target_annotation(annotation, geneid)
        gene_info = Gene.get_target_chr_start_end(target_annotation)
        gene_start = gene_info[1]
        gene_end = gene_info[2]
        gene_strand = gene_info[3]
        return gene_start, gene_end, gene_strand

    @staticmethod
    def get_common_region(target_annotation):
        # patentrna_list should not be None
        regions_all_mrna = []
        for va in target_annotation:
            if va._type == 'exon':
                if va._start != va._end:
                    regions_all_mrna.append([va._start, va._end])
        if len(regions_all_mrna) == 0:
            return None
        regions_all_mrna.sort(key=lambda x: x[0])
        single_blocks = [regions_all_mrna[0]]
        for i in range(1, len(regions_all_mrna)):
            if regions_all_mrna[i] not in single_blocks:
                single_blocks.append(regions_all_mrna[i])
        common_blocks = []
        m = len(single_blocks)
        temp = single_blocks[0]
        if m > 1:
            for i in range(m - 1):
                if temp[1] < single_blocks[i + 1][0] - 1:
                    common_blocks.append(temp)
                    temp = single_blocks[i + 1]
                else:
                    if temp[1] < single_blocks[i + 1][1]:
                        temp = [temp[0], single_blocks[i + 1][1]]
        common_blocks.append(temp)
        if Gene.check_true_common(common_blocks) is False:
            raise RuntimeError('error happened in common region process')
        return common_blocks

    @staticmethod
    def get_region_gene(chrom, commblock, samples):
        region_gene = []
        common_position = []
        for sample_block in samples:
            one_sample = []
            for start_end in commblock:
                one_sample += sample_block._sample.values('chr{}'.format(chrom), start_end[0], start_end[1])
                if samples.index(sample_block) == 0:
                    common_position += list(range(start_end[0], start_end[1]))
            if len(one_sample) != len(common_position):
                raise RuntimeError('common sample and position should have same length')
            region_gene.append(one_sample)
        return region_gene, common_position

    @staticmethod
    def check_true_common(commblock):
        for i in range(len(commblock)-1):
            if commblock[i][1] >= commblock[i+1][0]-1:
                return False
        return True

    @staticmethod
    def get_whole_gene(chrom, start, end, samples):
        """
        :return: a list of whole gene at samples. list[i] = sample i
        """
        whole_gene = [sample_block._sample.values('chr{}'.format(chrom), start, end)
                      for sample_block in samples]
        return whole_gene

    @staticmethod
    def check_value(whole_gene):
        """
        :return: True means starting resolution, normalization and detection. False means skip this gene.
        """
        for i in range(len(whole_gene)-1):
            if sum(whole_gene[i]) <= 1:
                return False
        return True

    # @staticmethod
    # def cut_common_zero_region(whole_gene, gene_position):
    #     sum_whole_gene = []
    #     for i in range(len(whole_gene[0])):
    #         sum_temp = 0
    #         for j in range(len(whole_gene)):
    #             sum_temp += whole_gene[j][i]
    #         sum_whole_gene.append(sum_temp)
    #     if len(whole_gene[0]) != len(gene_position):
    #         raise RuntimeError('whole gene and position should have same length')
    #     aft_cut = Resolution.cut_regions(sum_whole_gene, 20)
    #     if aft_cut is None:
    #         return whole_gene, gene_position
    #     cut_start = aft_cut[0]
    #     cut_end = aft_cut[1]
    #     whole_gene_cut = []
    #     position_after_cut = gene_position[0:cut_start[0]]
    #     for j in range(len(whole_gene)):
    #         temppp = whole_gene[j][0:cut_start[0]]
    #         for i in range(len(cut_start) - 1):
    #             temppp += whole_gene[j][cut_end[i]+1:cut_start[i + 1]]
    #             if j == 0:
    #                 position_after_cut += gene_position[cut_end[i]+1:cut_start[i + 1]]
    #         if len(temppp) <= 10:
    #             logging.debug('data too small after zero cutting.')
    #             return None, None
    #         whole_gene_cut.append(temppp + whole_gene[j][cut_end[-1]+1:])
    #     return whole_gene_cut, position_after_cut + gene_position[cut_end[-1]+1:]
