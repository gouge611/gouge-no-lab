#!/usr/bin/env python
# coding:UTF-8
from __future__ import division, print_function
from optparse import OptionParser
import sys

class Interval:
    def __init__(self, start, end):
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start + 1
def isOverlap(interval1, interval2):
    if max(interval1.start, interval2.start) <= min(interval1.end, interval2.end):
        return True
    return False

class ParseGene:
    def __init__(self, line):
        gene_info = line.strip().split('\t')
        self.symbol = gene_info[0]
        self.chrom = gene_info[1]
        self.strand = gene_info[2]

        _exon_start = map(int, gene_info[3].strip().split(','))
        _exon_end = map(int, gene_info[4].strip().split(','))
        self.exon_interval = []
        if len(_exon_start) == len(_exon_end):
            for i in range(0, len(_exon_start)):
                self.exon_interval.append(Interval(_exon_start[i], _exon_end[i]))

        _intron_start = map(int, gene_info[5].strip().split(','))
        _intron_end = map(int, gene_info[6].strip().split(','))
        self.intron_interval = []
        if len(_intron_start) == len(_intron_end):
            for j in range(0, len(_intron_start)):
                self.intron_interval.append(Interval(_intron_start[j], _intron_end[j]))

        self.exon_length = sum(map(int, gene_info[7].strip().split(',')))
        self.intron_length = sum(map(int, gene_info[8].strip().split(',')))

        _i_read_index = range(12, len(gene_info), 15)  # 13th column
        self.intron_read = []
        for i in _i_read_index:
            self.intron_read.append(int(gene_info[i]))
        _i_each_read_index = range(13, len(gene_info), 15)
        self.intron_each_read = []
        for i in _i_each_read_index:
            self.intron_each_read.append(map(int, gene_info[i].strip().split(',')))

        _e_read_index = range(9, len(gene_info), 15)  # 10th column
        self.exon_read = []
        for j in _e_read_index:
            self.exon_read.append(int(gene_info[j]))
        _e_each_read_index = range(10, len(gene_info), 15)
        self.exon_each_read = []
        for j in _e_each_read_index:
            self.exon_each_read.append(map(int, gene_info[j].strip().split(',')))

    def output(self):
        str_i = ','.join(map(str, self.intron_read))
        str_e = ','.join(map(str, self.exon_read))
        string = '\t'.join(map(str, [self.symbol, self.chrom, self.strand, self.exon_length, self.intron_length,
                                     str_e, str_i]))
        return string
    def output_norm(self):
        exon_read_norm = [x*1000/self.exon_length for x in self.exon_read]
        intron_read_norm = [x*1000/self.intron_length for x in self.intron_read]
        str_i = ','.join(map(str, intron_read_norm))
        str_e = ','.join(map(str, exon_read_norm))
        string = '\t'.join(map(str, [self.symbol, self.chrom, self.strand, str_e, str_i]))
        return string
    def output_tailor(self, intron_retention=True):
        """
        Optimalize output style
        27/3/2017 update: output could be intron retention or expression
        """
        exon_read_norm = [x*1000/self.exon_length for x in self.exon_read]
        intron_read_norm = [x*1000/self.intron_length for x in self.intron_read]
        str_i = ','.join(map(str, intron_read_norm))
        str_e = ','.join(map(str, exon_read_norm))
        string = '\t'.join(map(str, [self.symbol, self.chrom, self.strand])) + '\t'
        if intron_retention:
            string += str_i
        else:
            string += str_e
        return string

class PairGenes():
    def __init__(self, pairStr):
        self.s_genes = []
        self.as_genes = []
        for i in pairStr:
            gene = ParseGene(i)
            if '-AS' in i:  # antisense:'-AS'; minus_strand:'\t-\t'
                self.as_genes.append(gene)
            else:
                self.s_genes.append(gene)
        self.test_zhgene_ir=self.calc_ir_count(self.s_genes, self.as_genes)
        self.test_zhgene_expr=self.calc_ex_count(self.s_genes, self.as_genes)
        self.test_fugene_expr=self.calc_ex_count_for_fugene(self.s_genes, self.as_genes)
    def calc_ir_count(self, zhgene, fugene):
        """
        30/3/2017 update
        :param zhgene:
        :param fugene:
        :return:
        """
        count_list = []
        if zhgene.__len__() == fugene.__len__() == 1:
            overlap_list = []
            for i in zhgene[0].intron_interval:
                for e in fugene[0].exon_interval:
                    if isOverlap(i, e):
                        overlap_list.append(True)
                        break
                    else:
                        overlap_list.append(False)

            length = 0
            for index in range(0, len(zhgene[0].intron_interval)):
                if not overlap_list[index]:
                    length += zhgene[0].intron_interval[index].length

            for sample in zhgene[0].intron_each_read:
                count = 0
                for index in range(0, len(sample)):
                    if not overlap_list[index]:
                        count += sample[index] * 1000 / length
                count_list.append(count)
        return count_list
    def calc_ex_count(self, zhgene, fugene):
        """
        30/3/2017 update
        :param zhgene:
        :param fugene:
        :return:
        """
        count_list = []
        if zhgene.__len__() == fugene.__len__() == 1:
            overlap_list = []
            for i in zhgene[0].exon_interval:
                for e in fugene[0].exon_interval:
                    if isOverlap(i, e):
                        overlap_list.append(True)
                        break
                    else:
                        overlap_list.append(False)

            length = 0
            for index in range(0, len(zhgene[0].exon_interval)):
                if not overlap_list[index]:
                    length += zhgene[0].exon_interval[index].length

            for sample in zhgene[0].exon_each_read:
                count = 0
                for index in range(0, len(sample)):
                    if not overlap_list[index]:
                        count += sample[index] * 1000 / length
                count_list.append(count)
        return count_list
    def calc_ex_count_for_fugene(self, zhgene, fugene):
        """
        30/3/2017 update
        :param zhgene:
        :param fugene:
        :return:
        """
        count_list = []
        if zhgene.__len__() == fugene.__len__() == 1:
            overlap_list = []
            for i in fugene[0].exon_interval:
                for e in zhgene[0].exon_interval:
                    if isOverlap(i, e):
                        overlap_list.append(True)
                        break
                    else:
                        overlap_list.append(False)

            length = 0
            for index in range(0, len(fugene[0].exon_interval)):
                if not overlap_list[index]:
                    length += fugene[0].exon_interval[index].length

            for sample in fugene[0].exon_each_read:
                count = 0
                for index in range(0, len(sample)):
                    if not overlap_list[index]:
                        count += sample[index] * 1000 / length
                count_list.append(count)
        return count_list

    def output(self):
        string = ''
        for i in self.s_genes:
            string += i.output_tailor(intron_retention=False) + '\t'
        for j in self.as_genes:
            string += j.output_tailor(intron_retention=False) + '\t'
        return string

    def output_all(self):
        string = ''
        for i in self.s_genes:
            string += i.output_norm() + '\t'
        for j in self.as_genes:
            string += j.output_norm() + '\t'
        string += ','.join(map(str, self.test_zhgene_expr)) + '\t'
        string += ','.join(map(str, self.test_zhgene_ir)) + '\t'
        string += ','.join(map(str, self.test_fugene_expr)) + '\t'
        return string

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", "--input", action="store", dest="filename", type="string",
                      help="for combination data clearing downstream analysis")
    (options, args) = parser.parse_args()
    f = open("human_known_as_combination").read().split('#')
    for block in f:
        genePairStr = block.strip().split('\n')
        genePairStr = [x for x in genePairStr if x not in ['']]
        pair = PairGenes(genePairStr)
        if pair.s_genes.__len__() == 1 & pair.as_genes.__len__() == 1:
            sys.stdout.write(pair.output_all() + '\n')
        # print(pair.output_all())
        # print(pair.s_genes.__len__(),pair.as_genes.__len__())
