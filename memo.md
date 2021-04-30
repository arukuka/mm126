
# 実装の流れ

## 問題に触れて

- AHC 002 の反省から何も考えず beam stack search する、というのは避けたい
- AHC 002 のように解を何らかで作って変更すると良さそう？
  - AHC 002 に引きずられすぎ…？
- 大きい `color` を先に hole に入れたい
  - これをするため小さいのをどかしながら入れる実装も必要そう
  - けど難しそう
- とりあえず「解を何らかで作る」を実装したい
- ~~`color=1` (ビジュアライザで真っ黒に表示される block)は動かせない~~
- debug 用に input を取り出す
  - `java -jar tester.jar -exec "tee input.txt" -seed 1`
- CodeLLDB の break point や watch 式の使い方を学ぶ
  - https://github.com/vadimcn/vscode-lldb/blob/v1.6.2/MANUAL.md
  - `/nat item.r == 0 && item.c == 12` と使う

## 初期実装

- 5bd378935065fcecbd487eebe32fee340f9a6e2c
- `color > 1` のものを priority queue (main queue) に詰める
  - `color=0` を障害物、自分以外の block は存在しないものとして、move のみで hole に移動してみる
  - このときかかった試行回数が小さい順、同じ場合は color の大きい順で取り出せるようにする
- main queue からひとつずつ取り出して解を作成する
  - このときは Slide も使うようにする
- main queue は 1 iteration で消化されていそうな様子
  - 初期解の生成や、残りを埋めるときに使えそう

### 次に向けて

- `color=1` はビジュアライザで真っ黒に表示されるけど動かせる
- 完成された解をいじる実装がしたい
  - 近傍は
    - 途中の状態から color の大きいやつを先に入れてみる
    - 抜いてみる（これもスコアあがりそう）
    - 足してみる
  - (color の小さいものを) 抜いて最後に足す、のほうが良いかな
- block をどかすだけでなく、slide の受け止めのために移動するという方法もある
  - それができたらめちゃくちゃ格好いいな
