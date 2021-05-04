#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tempfile
import subprocess
import pathlib
import os
import re
import pandas as pd
import numpy as np
import sys
import scipy.stats


def run(seed):
    with tempfile.TemporaryDirectory() as dname:
        input_file_path = pathlib.Path(dname) / pathlib.Path("input.txt")
        subprocess.run(['java',
                        '-jar',
                        'tester.jar',
                        '-exec',
                        'tee {}'.format(str(input_file_path.absolute())),
                        '-seed',
                        str(seed),
                        '-novis'],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        proc = subprocess.run(['java',
                               '-jar',
                               'tester.jar',
                               '-exec',
                               'build/Slider cal_max_score',
                               '-seed',
                               str(seed),
                               '-novis'],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
        max_score = int([x for x in proc.stdout.decode(
            sys.stdout.encoding).split(os.linesep) if re.match(r'^\d+$', x)][0])
        proc = subprocess.run(['java',
                               '-jar',
                               'tester.jar',
                               '-exec',
                               'build/Slider',
                               '-seed',
                               str(seed),
                               '-novis'],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
        score = int(float([x for x in proc.stdout.decode(
            sys.stdout.encoding).split(os.linesep) if 'Score' in x][0].split()[-1]))
        return (score, max_score)


def build():
    id = subprocess.run(
        "git rev-parse --short HEAD",
        shell=True,
        stdout=subprocess.PIPE).stdout.decode(sys.stdout.encoding).strip()
    changed = int(
        subprocess.run(
            "git diff | wc -l",
            shell=True,
            stdout=subprocess.PIPE).stdout.decode(
            sys.stdout.encoding).strip()) > 0
    if changed:
        id = id + '*'
    subprocess.run(
        "cmake --build build --config Release --target all",
        shell=True)
    return id


def main():
    N = 100
    seeds_range = range(1, N + 1)
    current_id = build()
    print(f'current_id: {current_id}')

    initial_columns = ['seeds', 'max', 'vs max', 'vs others', 'win']
    initial_columns_type = [
        object,
        np.int32,
        np.float64,
        np.float64,
        np.float64]
    result = pd.DataFrame(
        np.arange((len(seeds_range) + 1) * 2)
        .reshape((len(seeds_range) + 1), 2),
        columns=[initial_columns[0], current_id])

    result.loc[:, 'seeds'] = list(seeds_range) + ["-"]

    result_csv_path = pathlib.Path("result.csv")
    others_columns = []
    if result_csv_path.exists():
        prev_result = pd.read_csv(result_csv_path, index_col=0)
        del_targets = initial_columns
        for col in prev_result.columns:
            if col[-1] == '*' or col == current_id:
                del_targets.append(col)
        prev_result = prev_result.drop(columns=del_targets)
        result = pd.concat((result, prev_result), axis=1)
        others_columns = others_columns + prev_result.columns.to_list()

    for col, dtype in zip(initial_columns[1:], initial_columns_type[1:]):
        result[col] = 0
        result[col] = result[col].astype(dtype)

    score_acc = 0
    max_score_acc = 0
    vs_max_acc = 0
    vs_others_acc = 0
    win_acc = 0
    test_win = 0
    test_all = 0
    for seed in seeds_range:
        score, max_score = run(seed)

        vs_max_ratio = score / max_score

        max_others = 1
        for other in others_columns:
            other_score = result.at[seed - 1, other]
            max_others = max(max_others, other_score)
        vs_others = score / max_others
        win = 2
        test_all = test_all + 1
        if score < max_others:
            win = 0
        elif score == max_others:
            win = 1
            test_all = test_all - 1
        else:
            win = 2
            test_win = test_win + 1
        win2str = ['x', '-', 'o']

        score_acc = score_acc + score
        max_score_acc = max_score_acc + max_score
        vs_max_acc = vs_max_acc + vs_max_ratio
        vs_others_acc = vs_others_acc + vs_others
        win_acc = win_acc + win / 2

        result.at[seed - 1, current_id] = score
        for col, val in zip(initial_columns[1:], [
                            max_score, vs_max_ratio, vs_others, win_acc / seed]):
            result.at[seed - 1, col] = val

        p_val = scipy.stats.binom_test(test_win, test_all)
        print(
            f'({seed:03}) {win2str[win]} {win_acc / seed:.3f} p={p_val:.3f} {score:>9,}',
            f'| {max_score:>9,}({vs_max_ratio:.3f}, {vs_max_acc / seed:.3f})',
            f'| {max_others:>9,}({vs_others:.3f}, {vs_others_acc / seed:.3f})')

    score_acc = score_acc / N
    max_score_acc = max_score_acc / N
    vs_max_acc = vs_max_acc / N
    vs_others_acc = vs_others_acc / N
    win_acc = win_acc / N

    result.at[N, current_id] = score_acc
    for col, val in zip(initial_columns[1:], [
                        max_score_acc, vs_max_acc, vs_others_acc, win_acc]):
        result.at[N, col] = val

    result.to_csv(result_csv_path)


if __name__ == "__main__":
    main()
