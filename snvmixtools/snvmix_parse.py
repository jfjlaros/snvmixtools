#!/usr/bin/env python

class SNVMixRecord(object):
    def __init__(self, line):
        """
        :param line: a line of SNVMix2 output
        :type line: str
        """
        data = line.strip().split()

        location = data[0].split(':')
        self.chromosome = location[0]
        self.position = int(location[1])

        self.reference = data[1]
        self.alternative = data[2]

        details = data[3].split(',')
        self.reference_count = int(details[0].split(':')[1])
        self.alternative_count = int(details[1].split(':')[1])
        self.genotype_likelihood = map(float, details[2:5])
        self.genotype = int(details[5])
    #__init__

    def __str__(self):
        return "{}:{}        {}       {}       {}:{},{}:{},{},{}\n".format(
            self.chromosome, self.position, self.reference, self.alternative,
            self.reference, self.reference_count, self.alternative,
            self.alternative_count,
            ",".join(map("{:.10f}".format, self.genotype_likelihood)),
            self.genotype)
#SNVMixRecord

def walker(handle):
    """
    """
    for line in handle.readlines():
        yield SNVMixRecord(line)
#walker
