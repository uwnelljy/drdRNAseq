class AnnoBlock:
    def __init__(self, chrid, type, start, end, strand, parentgene):
        self._chrid = chrid
        self._type = type
        self._start = start
        self._end = end
        self._strand = strand
        self._parentgene = parentgene


class SampleBlock:
    def __init__(self, sample):
        self._sample = sample
