
function [] = snv_stats(variants, varargin)

purines = 'AG';
pyrimidines = 'CT';

transversions = { 'A>C', 'C>A', 'A>T', 'T>A', 'G>C', 'C>G', 'G>T', 'T>G' };
transitions = { 'A>G', 'G>A', 'C>T', 'T>C' };



