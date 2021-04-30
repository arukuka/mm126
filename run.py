#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tempfile
import subprocess
import pathlib
import os
import re


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
                        '-novis'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
        proc = subprocess.run(['java',
                               '-jar',
                               'tester.jar',
                               '-exec',
                               'build/Slider cal_max_score',
                               '-seed',
                               str(seed),
                               '-novis'], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        max_score = int([x for x in proc.stdout.decode('utf-8').split(os.linesep) if re.match(r'^\d+$', x)][0])
        proc = subprocess.run(
            ['java', '-jar', 'tester.jar', '-exec',
                'build/Slider', '-seed', str(seed), '-novis'], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
        score = int(float([x for x in proc.stdout.decode('utf-8').split(os.linesep) if 'Score' in x][0].split()[-1]))
        ratio = score / max_score
        print(f'({seed:02}) {score:>7} / {max_score:>7} = {ratio}')


def main():
    for ite in range(1, 11):
        run(ite)


if __name__ == "__main__":
    main()
