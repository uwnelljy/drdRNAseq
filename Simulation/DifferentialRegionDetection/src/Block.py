class AnnoBlock:
    def __init__(self, chr, type, genestart, geneend, isostart, isoend, strand):
        self._chr = chr
        self._type = type
        self._genestart = genestart
        self._geneend = geneend
        self._isostart = isostart
        self._isoend = isoend
        self._strand = strand


class SampleBlock:
    def __init__(self, sample):
        self._sample = sample


class Block:
    def __init__(self, iso, value):
        self.iso = iso
        self.value = value