# 2023/06/28
- 杭挿入のコードを書く
- 初期変数をまとめて出力するためのメソッドがあるはず。
    - [`Omega.save`](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Omega.save)が該当
- ただ外部で定義した変数を上記のメソッドで保存する方法はないみたい
- [ここ](https://gitlab.com/yade-dev/trunk/blob/master/examples/rod-penetration/model.py)にロッドを挿入したサンプルコードが書かれている
- また[ここ](https://gitlab.com/yade-dev/trunk/blob/master/examples/PeriodicBoundaries/periodicSandPile.py)に横方向に周期境界を持つ場合に重力堆積させたコードがある
    - このコードでは全体を周期境界としつつも、底面にだけごくごく薄いBoxを挿入することで底面がある箱を表現している




# 2023/06/22

## そもそも杭挿入ができるのか
- `shape`の中にいくつか興味深いサブクラスがある
    - 例えば`cylinder`や`Facet`、`Box`など
    - 一つずつ見ていこう
- `cylinder`について
    - https://www.youtube.com/watch?v=Hh6nGzIU1vU&ab_channel=JanekKozicki
    - こういうもの。形状はものすごく近い
