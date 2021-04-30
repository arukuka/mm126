#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tempfile
import subprocess
import pathlib
import os
import re
import pandas as pd
import numpy as np


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
            'utf-8').split(os.linesep) if re.match(r'^\d+$', x)][0])
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
            'utf-8').split(os.linesep) if 'Score' in x][0].split()[-1]))
        return (score, max_score)


def build():
    id = subprocess.run(
        "git rev-parse --short HEAD",
        shell=True,
        stdout=subprocess.PIPE).stdout.strip()
    changed = int(
        subprocess.run(
            "git diff | wc -l",
            shell=True,
            stdout=subprocess.PIPE).stdout.strip()) > 0
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

    result = pd.DataFrame(
        np.arange((len(seeds_range) + 1) * 2)
        .reshape((len(seeds_range) + 1), 2),
        columns=['seeds', current_id])

    result_csv_path = pathlib.Path("result.csv")
    others_columns = []
    if result_csv_path.exists():
        prev_result = pd.read_csv(result_csv_path)
        del_targets = ['seeds', 'max', 'vs max', 'vs others', 'win']
        for col in prev_result.columns:
            if col[-1] == '*':
                del_targets.append(col)
        prev_result.drop(columns=del_targets)
        result = pd.concat(result, prev_result, axis=1)
        others_columns.append(prev_result.columns)

    result[['max', 'vs max', 'vs others', 'win']] = 0

    score_acc = 0
    vs_max_acc = 0
    vs_others_acc = 0
    win_acc = 0
    for seed in seeds_range:
        score, max_score = run(seed)

        vs_max_ratio = score / max_score
        print(f'({seed:02}) {score:>7} / {max_score:>7} = {vs_max_ratio}')

        win = True
        max_all = score
        for other in others_columns:
            other_score = result.at[seed - 1, other]
            win = win and other_score <= score
            max_all = max(max_all, other_score)
        vs_others = score / max_all

        score_acc = score_acc + score
        vs_max_acc = vs_max_acc + vs_max_ratio
        vs_others_acc = vs_others_acc + vs_others
        if win:
            win_acc = win_acc + 1

        result.at[seed - 1, current_id] = score
        result.at[]

        print(
            f'({seed:02}) {"o" if win else "x"} {win_acc / seed:.3f} {score:>9,}',
            f'| {max_score:>9,}({vs_max_ratio:.3f}, {vs_max_acc / seed:.3f})',
            f'| {max_all:>9,}({vs_others:.3f}, {vs_others_acc / seed:.3f})')

    score_acc = score_acc / N
    vs_max_acc = vs_max_acc / N
    vs_others_acc = vs_others_acc / N
    win_acc = win_acc / N
    result.to_csv(result_csv_path)


if __name__ == "__main__":
    main()
