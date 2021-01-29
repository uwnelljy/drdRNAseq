class AnnoBlock:
    def __init__(self, ischr, chrid1, chrid2, type, start, end, strand, id, parentgene):
        self._ischr = ischr
        self._chrid1 = chrid1
        self._chrid2 = chrid2
        self._type = type
        self._start = start
        self._end = end
        self._strand = strand
        self._id = id
        self._parentgene = parentgene


class SampleBlock:
    def __init__(self, sample):
        self._sample = sample
