import gzip
import json
import logging as logger
import numpy as np
import os

import pandas as pd
import polars as pl
import pysam

from spectre.analysis.coverage_stats import CoverageStatistics
from spectre.analysis.coverage_stats import CoverageData
from spectre.analysis.call_cnv_coverage import CNVCall
from spectre.analysis.call_cnv_AF import CNVCall as CNVAnalysisSNP
from spectre.analysis.cnv_metrics import CNVMetrics
from spectre.plots.plot import CoveragePlot, CNVPlot
from spectre.util import outputWriter
from spectre.util import vcf_parser
from spectre.util.dataAnalyzer import NormaldataAnalyser as NorAn
from spectre.util.cnv_id import CNV_ID
from spectre.util.OSUtil import OSUtil as OSUt

MOSDEPTH_HEADER = ['chrom_', 'start_', 'end_', 'coverage_']
MOSDEPTH_DTYPES = [pl.Utf8, pl.Int64, pl.Int64, pl.Float64]

class CNVAnalysis(object):
    def __init__(self, coverage_file, coverage_mosdepth_data, bin_size, output_directory, outbed, outvcf, genome_info,
                 sample_id, metadata_ref, snv_file, only_chr_list="", ploidy=2,min_cnv_len=1000000, as_dev=False,
                 dev_params=None, debug_dir="", dist_proportion=0.25, threshhold_quantile=5):
        # input
        self.coverage_file = coverage_file
        self.mosdepth_data = coverage_mosdepth_data  # has chr_mean_coverage, genome_mean_coverage, genome_bases
        self.output_directory = output_directory
        self.output_bed = outbed
        self.output_vcf = outvcf
        self.sample_id = sample_id
        self.snv_file = snv_file
        self.population_sample_ids = set()  # hold the filename of the population samples
        # bin size !important
        self.bin_size = int(bin_size)
        self.genome_info = genome_info
        self.metadata = metadata_ref
        # whitelist, only this chr
        self.only_chr_list = str(only_chr_list).split(",") if only_chr_list != "" else genome_info["chromosomes"]
        self.ploidy = ploidy
        # logger
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        # work variables
        self.positions = np.zeros(0)
        self.coverage = np.zeros(0)
        self.coverage_log2 = np.zeros(0)
        # results by chromosome
        self.genome_analysis = {}  # use this for CNV call NOTE: here is the coverage data under "cov_data"
        self.cnv_calls_list = {}  # use this for CNV call
        self.raw_cnv_calls_list = {}  # use this for storing raw CNV call
        self.intermediate_candidates_file_location = ""  # holds output path of serialized cnv object
        self.cnv_merged = {}  # use this for CNV call
        self.existing_cnv_ids = []  # holds all already used cnv IDs
        # snv data
        self.snv_af_bed = ""
        self.snv_af_df = None
        # dev/debug + hidden params
        self.as_dev = as_dev
        self.debug_dir = debug_dir

        # windows stats for charts
        self.window_stats_dir = f'{self.output_directory}/windows_stats'
        self.windows_bins = 50
        # cnv metrics
        self.cnv_metrics = None
 
        self.coverage_full_path = os.path.abspath(os.path.expanduser(self.coverage_file))

        # TODO
        self.min_chr_length = 1e6
        self.max_std_outlier_rm = 5
        self.max_std_outlier_cn = 10
        self.mosdepth_cov_genome_chr_diff = 0.10  # 10%
        self.lower_2n_threshold = 0.0  # are overwritten after data_normalization was called
        self.upper_2n_threshold = 0.0  # are overwritten after data_normalization was called
        self.cov_diff_threshold = 0.80
        self.dist_proportion = dist_proportion # 0.25
        self.dist_min_overwrite = 10000  # 10kb
        self.candidate_final_threshold = min_cnv_len #100000  # 100kb
        self.threshhold_quantile = threshhold_quantile
        self.no_call_lower_threshhold = 0.1 # fraction from genome median which chromosome median is compared to to be skiped
        self.detect_chr_y_threshhold = 0.1 # fraction from genome median which chromosome median is compared to
        self.chr_y_q = 0.8 # median seems to be 0 even for males, so need to use quntile
        self.chr_x_single_flag = False # will be initialized in data_normalization
        self.chr_y_present_flag = False # will be initialized in data_normalization

    def genome_median(self):
        return self.coverages_df_diploid['coverage_'].median()

    @staticmethod
    def annotate_bins_df(snps_df, EPS = 0.08, EXPECTED_AF = 0.5):
        # filter homozygous
        snps_df = snps_df.filter(0.5 - abs(0.5 - pl.col('af_')) > EPS)

        # af_good predictor
        snps_df = snps_df.with_columns(
            (abs(pl.col('af_') - EXPECTED_AF) < EPS).alias("af_good"))

        # aggragation on bin level (start column)
        af_good_bins_df = snps_df.group_by('start_', 'chrom_').agg(
            pl.col("af_good").all()
        )

        return af_good_bins_df

    def vcf_based_af_bins_annotation(self):
        self.logger.info("Parsing VCF to AF freqs file")
        vcf = vcf_parser.VCFSNVParser(self.min_chr_length, self.as_dev)
        vcf_df = pl.from_pandas(vcf.vcf_to_dataframe(self.snv_file))
        vcf_df = vcf_df.with_columns(pl.col('start_').apply(lambda x: x // self.bin_size * self.bin_size))

        af_good_bins_df = self.annotate_bins_df(vcf_df)

        coverages_df = pl.read_csv(
            self.coverage_full_path,
            has_header=False,
            separator='\t',
            new_columns=MOSDEPTH_HEADER,
            dtypes=MOSDEPTH_DTYPES
        )
        coverages_df = coverages_df.join(af_good_bins_df, on=['chrom_', 'start_'], how='left')
        self.coverages_df = coverages_df.with_columns(pl.col("af_good").fill_null(False))

        self.coverages_df_diploid = self.coverages_df.filter(af_good=True)

    def detect_xy_chromosome_proidness(self, genome_median, coverage_lower_threshold):
        # Detecting chrY has some caverage
        chr_y_cov = self.coverages_df.filter(chrom_='chrY')['coverage_']
        if len(chr_y_cov) > 0:
            self.chr_y_present_flag = chr_y_cov.quantile(self.chr_y_q) > genome_median * self.detect_chr_y_threshhold

            self.logger.debug(f'chrY len: {len(chr_y_cov)}, chrY quantile({self.chr_y_q}): {chr_y_cov.quantile(self.chr_y_q)}')
            self.logger.debug(f'chrY present flag: {self.chr_y_present_flag}')

        # Detecting chrX has one copy
        chr_x_cov = self.coverages_df.filter(chrom_='chrX')['coverage_']
        if len(chr_x_cov) > 0:
            self.chr_x_single_flag = chr_x_cov.median() < coverage_lower_threshold

            self.logger.debug(f'chrX len: {len(chr_x_cov)}, chrX median: {chr_x_cov.median()}')
            self.logger.debug(f'Single chrX flag: {self.chr_x_single_flag}')

    # TODO: Refactor whole thing using polars
    # Data normalization
    def data_normalization(self):
        """
        Normalize single chromosome bins
        """
        self.logger.debug(self.coverage_file)
        file_size_lines = OSUt.get_lines_by_chromosome(self.coverage_file)
        if len(file_size_lines) == 0:
            self.logger.error("Empty file")

        coverage_file_handler = gzip.open(self.coverage_full_path, "rt") if "gz" in self.coverage_file \
            else open(self.coverage_full_path, "r")

        list_index = 0
        previous_chromosome = ""
        genome_median = self.genome_median()

        # Thresholds borders bounds
        coverage_lower_threshold = self.coverages_df_diploid['coverage_'].quantile(self.threshhold_quantile / 100)
        coverage_upper_threshold = self.coverages_df_diploid['coverage_'].quantile(1 - self.threshhold_quantile / 100)
        self.logger.info(f'Thresholds for coverage, lower (Q{self.threshhold_quantile}): {coverage_lower_threshold}, upper(Q{100 - self.threshhold_quantile}): {coverage_upper_threshold}')

        self.lower_2n_threshold = coverage_lower_threshold / genome_median * self.ploidy
        self.upper_2n_threshold = coverage_upper_threshold / genome_median * self.ploidy
        self.logger.info(f'Thresholds for normalized coverage, lower: {self.lower_2n_threshold}, upper: {self.upper_2n_threshold} (normalized by median: {genome_median})')

        self.detect_xy_chromosome_proidness(genome_median, coverage_lower_threshold)

        for line in coverage_file_handler:
            [chromosome, start, _, coverage] = line.rstrip("\n").split("\t")
            if chromosome in self.only_chr_list:
                start = int(start)
                coverage = float(coverage)
                if previous_chromosome == "":
                    # init
                    self.positions = np.zeros(file_size_lines[chromosome])
                    self.coverage = np.zeros(file_size_lines[chromosome])
                    self.coverage_log2 = np.zeros(file_size_lines[chromosome])
                    list_index = 0
                    previous_chromosome = chromosome
                    # first elem
                    self.positions[list_index] = start
                    self.coverage[list_index] = coverage
                    self.coverage_log2[list_index] = np.log2(coverage) if coverage != 0 else np.NaN  # use np.nanFUN
                    list_index += 1
                elif previous_chromosome != chromosome:
                    # analysis here
                    self.logger.debug(previous_chromosome)
                    self.__remove_n_region_by_chromosome(previous_chromosome)
                    cov_stats, norm_stats, cov_data = self.__normalization_and_statistics(previous_chromosome, genome_median)
                    self.genome_analysis[previous_chromosome] = {"cov_data": cov_data,
                                                                 "statistics": cov_stats,
                                                                 "norm_statistics": norm_stats}
                    # init new chromosome
                    self.positions = np.zeros(file_size_lines[chromosome])
                    self.coverage = np.zeros(file_size_lines[chromosome])
                    self.coverage_log2 = np.zeros(file_size_lines[chromosome])
                    list_index = 0
                    previous_chromosome = chromosome
                    # first elem
                    self.positions[list_index] = start
                    self.coverage[list_index] = coverage
                    self.coverage_log2[list_index] = np.log2(coverage) if coverage != 0 else np.NaN  # use np.nanFUN
                    list_index += 1
                else:
                    self.positions[list_index] = start
                    self.coverage[list_index] = coverage
                    self.coverage_log2[list_index] = np.log2(coverage) if coverage != 0 else np.NaN  # use np.nanFUN
                    list_index += 1
        coverage_file_handler.close()

        # compute leftover chromosome
        self.__remove_n_region_by_chromosome(previous_chromosome)
        cov_stats, norm_stats, cov_data = self.__normalization_and_statistics(previous_chromosome, genome_median)
        self.genome_analysis[previous_chromosome] = {"cov_data": cov_data,
                                                     "statistics": cov_stats,
                                                     "norm_statistics": norm_stats}

        # setup, calculating border for deletion and duplications
        self.cnv_metrics = CNVMetrics(genome_analysis=self.genome_analysis,
                                      cnv_calls=self.cnv_calls_list,
                                      exclusion_zones=self.metadata,
                                      hashname=self.sample_id,
                                      ploidy=self.ploidy,
                                      output_dir=self.output_directory,
                                      as_dev=self.as_dev, debug_dir=self.debug_dir)

        # clean up
        self.positions = np.zeros(0)
        self.coverage = np.zeros(0)
        self.coverage_log2 = np.zeros(0)

    def __normalization_and_statistics(self, chromosome, genome_median) -> tuple:
        self.logger.info(f'Number positions to be tested on chromosome {chromosome}: {self.coverage.size}')
        # clean up
        cov_stats = CoverageStatistics()
        norm_stats = CoverageStatistics()
        cov_data = CoverageData()
        cov_data.normalized_cov = np.array(0)
        if self.coverage.size > 0:
            # calculate avg, std, median and remove outliers (right tail)
            # use genome-wide average
            genome_avg_cov = self.mosdepth_data.genome_mean_coverage
            chromosome_med_coverage = np.nanmedian(self.coverage)
            skip_low_cov_chr_flag = chromosome_med_coverage < genome_median * self.no_call_lower_threshhold

            [avg, std] = [float(genome_avg_cov), np.nanstd(self.coverage)]
            for idx in range(0, self.coverage.size, 1):
                cov = self.coverage[idx]
                if (cov > (self.max_std_outlier_cn * avg) + (self.max_std_outlier_rm * std)) or skip_low_cov_chr_flag:
                    self.coverage[idx] = np.NaN

            # re-calculate avg, std, median without outliers (right tail)
            [avg, std, med] = [float(genome_avg_cov), np.nanstd(self.coverage), np.nanmedian(self.coverage)]

            self.logger.debug([f'genome avg: {genome_avg_cov}, median: {genome_median} | chromosome avg: {avg}, std: {std}, med: {med}'])

            if chromosome in self.genome_info["chr_lengths_by_name"]:
                # data
                cov_data.positions = self.positions
                cov_data.coverage_raw = self.coverage
                cov_data.coverage_log2 = self.coverage_log2
                # stats
                cov_stats.chromosome_len = self.genome_info["chr_lengths_by_name"][chromosome]
                cov_stats.chromosome_name = chromosome
                cov_stats.average = avg
                cov_stats.std_dev = std
                cov_stats.median = med
                cov_stats.min = np.nanmin(self.coverage)
                cov_stats.max = np.nanmax(self.coverage)

                # normalization, based on diploid organisms
                if (self.chr_x_single_flag and chromosome == 'chrX') or (self.chr_y_present_flag and chromosome == 'chrY'):
                    chr_ploidy = 1
                else:
                    chr_ploidy = self.ploidy
                normalize_by = genome_median  # med:median | avg:average
                normalized_candidates = NorAn.normalize_candidates(self.coverage, normalize_by) * self.ploidy / chr_ploidy
                normalized_candidates_ploidy = normalized_candidates * chr_ploidy
                if len(normalized_candidates) < 1:
                    normalized_candidates_ploidy = normalized_candidates

                # statistics for normalized bin
                avg, std, min_val, max_val, med = NorAn.get_candidate_statistics(normalized_candidates)
                norm_stats.chromosome_len = self.genome_info["chr_lengths_by_name"][chromosome]
                norm_stats.chromosome_name = chromosome
                norm_stats.median = med
                norm_stats.average = avg
                norm_stats.std_dev = std
                norm_stats.min = min_val
                norm_stats.max = max_val

                # norm data
                cov_data.normalized_cov = normalized_candidates
                cov_data.normalized_cov_ploidy = normalized_candidates_ploidy
            else:
                self.logger.warning(f"Chromosome:{chromosome} not found in gene info!")
        return cov_stats, norm_stats, cov_data

    def __remove_n_region_by_chromosome(self, chrom) -> None:
        """
        Compare CNV candidates with "N" regions from the reference genome. If the candidate starts or ends inside a
        "N" region, the candidate will is removed by assigning a np.nan value.
        :param chrom: current chromosome
        :return: None
        """
        # we use two variables which contain all the current chromosome info
        # self.positions
        # self.coverage
        if chrom in self.metadata:
            for meta_start, meta_end in self.metadata[chrom]:
                meta_start = int(meta_start)
                meta_end = int(meta_end)
                overlap_n = np.logical_and(self.positions >= meta_start, self.positions <= meta_end)
                self.coverage[overlap_n] = np.nan

    # *** MAIN ***
    def call_cnv_coverage(self, cnv_id_list=None, write_csv=False):
        if cnv_id_list is None:
            cnv_id_list = []
        # unique CNV IDs
        for each_chromosome in self.genome_analysis.keys():
            # global
            # Section 1: Generating raw CNV candidates
            self.logger.info(f"Calculating CNVs for {self.sample_id} @ chromosome {each_chromosome}")
            cnv_caller = CNVCall(self.as_dev)

            lower_threshold, upper_threshold = (
                min(self.lower_2n_threshold / 2, 0.5),
                max(self.upper_2n_threshold / 2, 1.5)
            ) if (self.chr_x_single_flag and each_chromosome == 'chrX') or (self.chr_y_present_flag and each_chromosome == 'chrY') else (
                self.lower_2n_threshold,
                self.upper_2n_threshold
            )
            candidates_cnv_list = cnv_caller.cnv_coverage(self.genome_analysis[each_chromosome]["cov_data"],
                                                          self.bin_size, each_chromosome, self.sample_id,
                                                          lower_threshold, upper_threshold)

            self.logger.debug([f'{c.start}-{c.end}, ({c.type},{c.size})' for c in candidates_cnv_list])
            # Generating CNV IDs
            new_ids = CNV_ID.n_id_generator(len(candidates_cnv_list), cnv_id_list)
            self.existing_cnv_ids = self.existing_cnv_ids + new_ids
            for candidate, id in zip(candidates_cnv_list, new_ids):
                candidate.set_id(id)

            # save candidates for each chromosome
            self.cnv_calls_list[each_chromosome] = candidates_cnv_list
        self.raw_cnv_calls_list = self.cnv_calls_list.copy()

    def refine_cnv_calls(self, write_csv=False):
        self.logger.info("refining cnv calls")
        for each_chromosome in self.genome_analysis.keys():
            candidates_cnv_list = self.cnv_calls_list[each_chromosome]
            if write_csv:
                self.dev_write_csv(each_chromosome)
            self.write_csv_window_agrregate(each_chromosome)

            self.logger.debug(f'total cnv candidates in {each_chromosome}: {len(candidates_cnv_list)} before merge')
            final_cnv_candidates = self.merge_candidates(candidates_cnv_list, each_chromosome)
            cnv_call_list = self.scaffold_candidates(final_cnv_candidates, each_chromosome) \
                if len(final_cnv_candidates) >= 2 else final_cnv_candidates

            self.logger.debug('Candidates after scaffolding:')
            self.logger.debug([f'{c.start}-{c.end}, ({c.type},{c.size})' for c in cnv_call_list])
            self.cnv_calls_list[each_chromosome] = cnv_call_list
        pass

    # Candidate CNV merge by distance and CNV type
    def merge_candidates(self, candidates_cnv_list, chromosome_name):
        if len(candidates_cnv_list) > 1:
            dev_candidates_string = ""
            merge_rounds = 1
            self.logger.debug(f'Current merge cycle {merge_rounds}')
            [n_merges, merged_candidates, _] = self.cnv_candidate_merge(candidates_cnv_list)
            while n_merges > 0:
                merge_rounds += 1
                self.logger.debug(f'Current merge cycle {merge_rounds}')
                self.logger.info(f'n candidates in {chromosome_name}: {len(merged_candidates)}')
                [n_merges, merged_candidates, dev_candidates_string] = self.cnv_candidate_merge(merged_candidates)
                self.logger.debug([f'{c.start}-{c.end}, ({c.type},{c.size})' for c in merged_candidates])

            self.logger.debug(f'Total merge rounds: {merge_rounds}')

            [n_merges, merged_candidates, _] = self.cnv_candidate_merge(cnv_candidates=merged_candidates,
                                                                        allow_over_blacklist_merging=True)
            self.logger.debug(f'Merged with final blacklist ignoring round: {n_merges}')

            candidates_cnv_list = self.clean_merged_candidates(merged_candidates)

            self.logger.debug('Candidates after merging and filtering:')
            self.logger.debug([f'{c.start}-{c.end}, ({c.type},{c.size})' for c in candidates_cnv_list])
            self.logger.debug(f'n candidates {chromosome_name}: {len(candidates_cnv_list)}')
        else:
            candidates_cnv_list = self.clean_merged_candidates(candidates_cnv_list)
        return candidates_cnv_list

    def cnv_candidate_merge(self, cnv_candidates, allow_over_blacklist_merging: bool = False):
        def dev_candidates_merge(header=False):  # dev
            if header:
                return f'chr\tleft\tright\tleft_size\tright_size\tcov_left\tcov_right\tdist_OK\ttype_OK\tCN_OK'
            else:
                return (f'{cnv_candidates[idx].chromosome}\t{end0}\t{start1}\t'
                        f'{cnv_cand.size}\t{cnv_candidates[idx + 1].size}\t'
                        f'{cnv_candidates[idx].median_cov_norm}\t{cnv_candidates[idx + 1].median_cov_norm}\t'
                        f'{(start1 - end0) < max_merge_distance}\t{same_type}\t{same_cnv_status}')

        n_candidates = len(cnv_candidates)
        # distance between end of i to start of i+1
        # first on is given
        idx = 0
        n_merges = 0
        merges_candidates = []
        cnv_cand = cnv_candidates[idx]
        cov_diff_threshold = self.cov_diff_threshold  # selected by trial and error in simulations
        # dist => distance
        dist_proportion = self.dist_proportion  # selected by trial and error, max distance proportion for merge
        dist_min_overwrite = self.dist_min_overwrite
        dev_candidates_string = [dev_candidates_merge(header=True)]
        while idx < n_candidates - 1:
            end0 = cnv_candidates[idx].end
            start1 = cnv_candidates[idx + 1].start
            same_type = cnv_candidates[idx].type == cnv_candidates[idx + 1].type
            # copy number status difference
            same_cnv_status = abs(
                cnv_candidates[idx].median_cov_norm - cnv_candidates[idx + 1].median_cov_norm) <= cov_diff_threshold
            # distance between merging windows is a % of the average length of the windows to be merged
            max_merge_distance = max(np.mean([cnv_cand.size, cnv_candidates[idx + 1].size]) * dist_proportion,
                                     dist_min_overwrite)
            dev_candidates_string.append(dev_candidates_merge())
            # Evaluate if the new span is covered by any blacklisted region
            # Check if span_over_blacklisted_region = true if no overlap exists
            if cnv_cand.chromosome not in self.metadata.keys():
                # No blacklisted regions for this chromosome
                span_over_blacklisted_region = False
            else:
                # Found chromosome in blacklist and needs to check all positions
                span_over_blacklisted_region = any(
                    any(end0 <= int(val) <= start1 for val in tup) for tup in self.metadata[cnv_cand.chromosome])

            if (start1 - end0) < max_merge_distance and same_type and same_cnv_status and \
                    (not span_over_blacklisted_region or allow_over_blacklist_merging):
                # merge and extend
                n_merges += 1
                cnv_cand.add_candidates(
                    cnv_candidates[idx + 1].pos,
                    cnv_candidates[idx + 1].cov,
                    cnv_candidates[idx + 1].id,
                    cnv_candidates[idx + 1].merged_sample_references
                )
            else:
                # save and show
                cnv_cand.median_coverage_candidates_merged()
                merges_candidates.append(cnv_cand)
                # cnv_cand.show_candidates()  # DEV
                # init
                cnv_cand = cnv_candidates[idx + 1]
            idx += 1
        # one last time for the remaining ones
        merges_candidates.append(cnv_cand)
        # cnv_cand.show_candidates()  # DEV
        return [n_merges, merges_candidates, "\n".join(dev_candidates_string)]

    def clean_merged_candidates(self, merged_candidates):
        # current threshold is 100kb
        clean_merged = []
        for each_cand in merged_candidates:
            if each_cand.size >= self.candidate_final_threshold:
                clean_merged.append(each_cand)
        return clean_merged

    # Fill in the gaps after merge
    def scaffold_candidates(self, final_candidates_list, chromosome_name):
        cnv_list = []
        n_scaffolds = 0
        cov_data_chr = self.genome_analysis[chromosome_name]["cov_data"]
        _index = 0
        while _index < len(final_candidates_list):
            _cand = final_candidates_list[_index]
            # get the difference between each position
            position_diff = np.diff(_cand.pos)
            # get index of the position_diff where the difference is not equal to the bin size
            position_diff_index = np.where(position_diff != self.bin_size)
            if len(position_diff_index[0]) > 0:
                # fill gaps for each position_diff_index and  reversing list with [::-1] to ensure no index overwritten
                # for idx1, idx2 in zip(position_diff_index[0][::-1][:-1],position_diff_index[0][::-1][1:]):
                for idx1 in position_diff_index[0][::-1]:

                    # Indices of current  cov gap in the current candidate
                    cnv_gap_start_idx = idx1
                    cnv_gap_end_idx = idx1 + 1

                    # Get chromosomal positions of the gap in the current candidate
                    gap_in_candidate_start = int(_cand.pos[cnv_gap_start_idx])
                    gap_in_candidate_end = int(_cand.pos[cnv_gap_end_idx])

                    # get gap indices in the coverage data
                    gap_start = int(np.where(cov_data_chr.positions == gap_in_candidate_start)[0][0])
                    gap_end = int(np.where(cov_data_chr.positions == gap_in_candidate_end)[0][0])

                    # get gap data
                    gap_data = cov_data_chr.normalized_cov_ploidy[gap_start:gap_end]
                    gap_pos = cov_data_chr.positions[gap_start:gap_end]
                    scaf_cov = np.nanmedian(gap_data)

                    # get array with scaf_cov values the length of the gap
                    scaf_cov_array = np.full(len(gap_data), scaf_cov)

                    # Restricting cnv scaffolding merge over 100 following N regions (NaN values)
                    scaf_cov_nan_overflow = np.count_nonzero(np.isnan(np.array(gap_data))) > 0

                    # Apply coverage to the gap if not to many NaN values are present
                    if not scaf_cov_nan_overflow:
                        _cand.add_scaffold_candidate(cnv_gap_start_idx, cnv_gap_end_idx, scaf_cov_array.tolist(),
                                                     gap_pos.tolist())
                        n_scaffolds += 1

            _cand.median_coverage_candidates_merged()
            _cand.reinitialize_candidate_values()
            _cand.set_gt()
            cnv_list.append(_cand)
            _index += 1

        return cnv_list

    # Output results
    def cnv_result_bed(self, method=""):
        if method == "":
            output_bed = self.output_bed
        else:
            output_bed = os.path.join(os.path.join(self.output_directory, f'{method}{self.sample_id}.bed'))

        bed_output = outputWriter.BedOutput(output_bed)
        bed_output.make_bed(self.genome_analysis.keys(), self.cnv_calls_list)

    def cnv_result_vcf(self, method=""):
        output_vcf = os.path.join(os.path.join(self.output_directory, f'{method}{self.sample_id}.vcf'))
        vcf_output = outputWriter.VCFOutput(output_vcf, self.genome_info)
        vcf_output.make_vcf(self.genome_analysis.keys(), self.cnv_calls_list, self.sample_id)

    def karyotype_json(self, method=""):
        output_male_flag = os.path.join(os.path.join(self.output_directory, f'expected_karyotype.json'))
        karyotype = {
                'X' : 1 if self.chr_x_single_flag else 2,
                'Y' : 1 if self.chr_y_present_flag else 0
                }
        with open(output_male_flag, 'w') as f:
            json.dump(karyotype, f)

    # Plots
    def coverage_plot(self):
        for each_chromosome in self.genome_analysis.keys():
            plot_cnv = CoveragePlot()
            plot_cnv.file_prefix = f'plot1-coverage-{self.sample_id}'
            plot_cnv.output_directory = self.output_directory
            chr_cov_data = self.genome_analysis[each_chromosome]["cov_data"]
            plot_cnv.plot_coverage(each_chromosome, {"pos": chr_cov_data.positions, "cov": chr_cov_data.raw_coverage})

    def cnv_plot(self, methode=""):
        for each_chromosome in self.genome_analysis.keys():
            chr_cov_data = self.genome_analysis[each_chromosome]["cov_data"]
            cov_stats = self.genome_analysis[each_chromosome]["norm_statistics"]
            lower_bound = self.lower_2n_threshold
            upper_bound = self.upper_2n_threshold
            cov_stats.median = cov_stats.median * self.ploidy  # for diploid
            new_plot_device = CNVPlot()
            new_plot_device.output_directory = self.output_directory
            new_plot_device.file_prefix = methode + self.sample_id
            new_plot_device.plot_coverage_cnv(each_chromosome, cov_stats,
                                              {"pos": chr_cov_data.positions,
                                               "cov": chr_cov_data.normalized_cov_ploidy},
                                              self.cnv_calls_list[each_chromosome], [lower_bound, upper_bound])

    def convert_vcf_to_tabular(self, snv_af_bed_output):
        self.logger.info("Parsing VCF to BED file")
        vcf = vcf_parser.VCFSNVParser(self.min_chr_length, self.as_dev)
        self.logger.debug("Parsing: VCF -> DataFrame")
        vcf_df = vcf.vcf_to_dataframe(self.snv_file)
        self.logger.debug("Parsing: DataFrame -> BED")
        self.snv_af_df = vcf.dataframe_to_tabular_file(vcf_df, self.coverage_file, snv_af_bed_output)
        self.snv_af_bed = snv_af_bed_output

    def call_cnv_af_region(self):
        # self.cnv_calls_list[each_chromosome]
        cnv_by_af = CNVAnalysisSNP(genome_info=self.genome_info,
                                   output_directory=self.output_directory,
                                   sample_id=self.sample_id,
                                   as_dev=self.as_dev)
        self.logger.info("Calculating CNV events based on SNV data")
        self.cnv_calls_list = cnv_by_af.af_cnv_call_region(self.cnv_calls_list, self.snv_af_bed)
        cnv_by_af.call_cnv_af(self.snv_af_df, write_csv=True)

    def get_cnv_metrics(self, refined_cnvs: bool = False):
        """
        Calculates for every CNV call a metrics like the p value with statistical tests
        :return:
        """
        self.cnv_metrics.evaluate_cnvs(self.cnv_calls_list, refined_cnvs)

    def write_intermediate_candidates(self, candidate_type: str = "") -> str:
        """
        Writing intermediate object files out
        :param candidate_type:
        :return:
        """

        intermediate_output_writer = outputWriter.IntermediateFile(self.output_directory)
        #genome_info = intermediate_output_writer.convert_candidates_to_dictionary(self.genome_info)
        cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(self.cnv_calls_list)
        raw_cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(self.raw_cnv_calls_list)

        analysis_dict = {
                "metadata": {"source":"spectre","spectre_version": "0.1"},
                "genome_info": self.genome_info,
                "raw_cnvs": raw_cnv_calls_list_dict,
                "refined_cnvs": cnv_calls_list_dict,
                "analysis_metrics": {
                    "min_chr_length": self.min_chr_length,
                    "max_std_outlier_rm": self.max_std_outlier_rm,
                    "mosdepth_cov_genome_chr_diff": self.mosdepth_cov_genome_chr_diff,
                    "lower_2n_threshold": self.lower_2n_threshold,
                    "upper_2n_threshold": self.upper_2n_threshold,
                    "cov_diff_threshold": self.cov_diff_threshold,
                    "dist_proportion": self.dist_proportion,
                    "dist_min_overwrite": self.dist_min_overwrite,
                    "candidate_final_threshold": self.candidate_final_threshold,
                    "genome_mean": np.mean(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage']),
                    "genome_sd": np.std(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage']),
                    "genome_var": np.var(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage'])
                }
            }
        output_path = intermediate_output_writer.write_intermediate_file(analysis_dict, f"{self.sample_id}")
        self.intermediate_candidates_file_location = output_path
        return output_path

    def collect_csv(self, chr):
        csv_results = pd.DataFrame(
            data={"position": self.genome_analysis[chr]["cov_data"].positions,
                  "mosdepth_cov": self.genome_analysis[chr]["cov_data"].coverage_raw,
                  "norm_cov": self.genome_analysis[chr]["cov_data"].normalized_cov,
                  "ploidy_cov": self.genome_analysis[chr]["cov_data"].normalized_cov_ploidy
                  }
        )
        csv_results["chr"] = chr
        return csv_results

    # ############################################
    # dev
    def dev_write_csv(self, chr):
        output_file = f"{self.debug_dir}/cnv_{self.sample_id}_byCoverage_chr{chr}.csv"
        self.logger.debug(f"Writing coverage to {output_file}")
        self.collect_csv(chr).to_csv(output_file, index=False)

    # ############################################
    # for chart in report
    def write_csv_window_agrregate(self, chr):
        window_stats = self.collect_csv(chr).shift(1).rolling(self.windows_bins, step=self.windows_bins).agg({
            'position': ['min', 'max'],
            'mosdepth_cov': 'mean',
            'norm_cov': 'mean',
            'ploidy_cov': 'mean'}
        ).iloc[1::]

        if not os.path.exists(self.window_stats_dir):
            os.makedirs(self.window_stats_dir)

        output_file = f"{self.window_stats_dir}/cnv_{self.sample_id}_windows_stats_chr{chr}.csv"
        self.logger.debug(f"Writing coverage to {output_file}")
        window_stats.to_csv(output_file, index=False)
