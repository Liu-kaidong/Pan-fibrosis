#!/usr/bin/env bash

conda activate cytotrace2-py
cd data/CytoTRACE2
cytotrace2 --input-path cell_expression.txt --annotation-path cell_annotations.txt --species human
