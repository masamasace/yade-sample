# 2023/08/22

## Paraviewの使い方
- 

# 2023/08/21

## 現時点での課題
- 初期状態で一番下にある球が周期教会で上に移動してしまう
- Wall-SphereのMaterialが設定されていない
- 杭の位置が不明
    - VTKのSpheresには収録されているようだ
    - よく見るとSpheresはCylinderとして認識されていない
    - 単に２つの球が埋め込まれているだけ...
    - Facetにしたほうが良さそう...
        - Facetには鉛直力を取り出すメソッドもある([例えばこれ](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Facet))
    - Facetで作ったほうが後々のためになりそう

# 2023/08/19

## 論文に向けて行いたいこと
- 取得するパラメータの整理
    - VTKファイルの出力
    - CSVファイルの出力項目の整理
        - 部分的な領域における応力とひずみ、密度、ファブリックテンソルに関するパラメータの出力

## 今後の論文で検討すること
- 杭先端の形状
    - 今回は半球が先端についている形でしかできないが、将来は円錐やフラットも試してみたい
- 打ち込み地盤の初期状態に関する検討
    - 地盤の初期密度、初期の粒子配置、粒度分布
- 異なる杭ー地盤間の摩擦特性、どのような接触モデルが適当か
- 杭の配置パターン
    - 今回は周期境界条件にすることで半無限地盤に同時に杭を打ち込んだ状態となる
- 液状化層厚に対しての貫入量の影響
- 杭の断面形状
- 実際の地震波による影響

# 2023/06/28

## 次に向けてのTODO
- 何回かに分けて粒子を堆積させて、初期高さと平均密度をコントロールする
## メモ
- 杭挿入のコードを書く
- 初期変数をまとめて出力するためのメソッドがあるはず。
    - [`Omega.save`](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Omega.save)が該当
- ただ外部で定義した変数を上記のメソッドで保存する方法はないみたい
- [ここ](https://gitlab.com/yade-dev/trunk/blob/master/examples/rod-penetration/model.py)にロッドを挿入したサンプルコードが書かれている
- また[ここ](https://gitlab.com/yade-dev/trunk/blob/master/examples/PeriodicBoundaries/periodicSandPile.py)に横方向に周期境界を持つ場合に重力堆積させたコードがある
    - このコードでは全体を周期境界としつつも、底面にだけごくごく薄いBoxを挿入することで底面がある箱を表現している
    - `allowBiggerThanPeriod=True`がないと、周期境界の要素よりも小さい底面の箱しかおけず、粒子がこぼれ落ちてしまう。
- どうやら`cylinder`を使うと`yade.gridpfacet`が呼び出されてしまうらしい
    - 詳しくは[こちら](https://yade-dem.org/doc/yade.gridpfacet.html?highlight=cylid)
- `cylinder`でどのように力が計算されているのかがわからない...
    - 要検討

# 2023/06/22

## そもそも杭挿入ができるのか
- `shape`の中にいくつか興味深いサブクラスがある
    - 例えば`cylinder`や`Facet`、`Box`など
    - 一つずつ見ていこう
- `cylinder`について
    - https://www.youtube.com/watch?v=Hh6nGzIU1vU&ab_channel=JanekKozicki
    - こういうもの。形状はものすごく近い
