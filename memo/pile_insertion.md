# 2023/09/15
- Paraviewの追加の使い方
    - 粒子の軌跡を描画する方法
        - https://www2.kaiyodai.ac.jp/~kentaro/materials/new_HP/particle_sim/16visualization_paraview.html
    - やっぱり部分領域の密度を求めたい
        - 方法1, CUDAで後処理
        - 方法2, Paraviewで後処理
    - 繰り返し堆積させる時間がもったいない
        - 重力堆積によるテンプレートを作成してしまう

# 2023/08/30
- `2023-08-23_14-14-07`のデータ整理
    - 杭の貫入開始のインデックスは808
    - facetsは杭
        - coloringをforceLenに変えるだけで杭のどの部分から応力を受けているかを可視化できる
        - ただfacetが細かく分割されすぎていて、よくわからない。同一高さの周方向に積分する必要がある
            - 要調査
    - intrsは相互作用力
        - indexを相互作用が一つでも生じているインデックスに移動
            - そうしないと次のフィルタが適用できない
        - `CellDatatoPointData`フィルタを適用
        - `Tube`フィルタを適用
            - `Radius`は基本となる円柱の半径の比率：ここでは1.0e-07ぐらいに
            - `Vary Radius`を`By Scalar`に
            - `Radius Factor`は最大値と最小値の円柱の半径の比率。ここでは100ぐらいにする
            - Coloringを変更するとよりわかりやすくなる
                - あまり円柱の半径に差が出ないので少しRangeを補正する必要があるかも
                    - 要調査
        - ある部分領域のRose Diagramの作成
            - `CellDatatoPointData`フィルタの下に`Clip`フィルタを適用する
                - boxにすることで部分領域(指定する値はあるおそらくboxの中でも一番小さいx, y, zの座標値)が指定される
            - 上部の`Extractors`→`Data`→`CSV`をクリック
            - `Properties`で出力の設定を入力
            - `File`→`Save Extracts...`をクリックして保存
            - Pythonで処理(要調査)
    - spheresは粒子単体
        - 球として表示するには`Glyph`フィルタをかける
            - `Glyph Type`は`Sphere`とする
            - `Scale Array`は`radii`とする
            - `Scale Factor`は2が実際のスケールと合う
            - `Masking`の`Glyph Mode`を`Uniform Spatial Destribution`から`All Points`にすると全部の点が表示される。
        - ある断面での変位分布、応力分布を知りたい
            - 変位はある時間ステップを初期値として、その位置からの相対変位を求めてあげる必要がある
            - `sphere`のデータセットに`Calculator`フィルタをかける
                - `Coordinate Results`のチェックは外しておく
                - `Result Array Name`とその下の入力はどちらも`coords`とする
            - `Calculator`フィルタの下に`ForceTime`フィルタをかける
                - `ForceTime`フィルタの`ForcedTime`の値を杭の貫入開始のインデックスとする
                    - 今回の場合は808
            - Ctrlキーを押しながら、`Force Time1`のフィルタと`Calculator1`フィルタをこの順で選択する
                - `vtu`ファイルを選択してはいけない
                    - `coords`のPointDataが生成されていないためエラーになる
            - 2つが選択された状態で`PythonCalculator`フィルタを適用する
            - `Expression`の欄に`inputs[1].PointData['coords']-inputs[0].PointData['coords']`を入力する
                - 注意：周期境界で解析をしているので、例えば+xの面から出ていった粒子は-xの面から入ってくる
                - このとき変位としては解析上x方向の解析領域の長さが足されてしまう
                - ここは計算方法を少し変える必要がある(要調査)
            - `vtu`から直接出てきた`Python Calculator`フィルタに`Glyph`フィルタを作用させる
                - `Glaph Source`と`Scale`のプロパティの項目は上と同じ
                - `Coloring`のプロパティは先ほど作成した`disp`とする
                    - 必要に応じて`Magunitude`を`x`などに変えるとその成分のみの変位が可視化できる
            - 断面の作成：`Clip`フィルタを上で作成した`Glyph`フィルタの下に生成
                - `Clip Type`は`Plane`とする
                - 切り取る平面の位置を`Plane Parameters`で設定
                - これによって指定断面における粒子単体の動きが可視化される
                    - ただ切り取られた粒子の断面が中空になってしまう→解決策については(要調査)
            - 内装された断面の作成：`vtu`から直接出てきた`Python Calculator`フィルタの下に`PointPlaneInterpolator`フィルタを作成
                - `Properties`の`Kernel`は`GaussianKernel`とする
                - `Radius`は粒子のサイズに5倍くらい
                - `Plane Bounds`の`Bounding Box`はよくわからない
                    - Y軸は杭の貫入開始時の粒子の堆積高さでいいかもしれない
                    - X軸とZ軸は解析領域の端でいいかもしれない
                - 切り取る平面の位置を`Plane Parameters`で設定
                    - `Center`は切り取る断面の中心位置
                    - `Normal`は法線ベクトルの向き
                        - 例えば左から1, 0, 0にするとx軸の方向が法線ベクトルの平面になる
                - `coloring`を変えることによって、表示させるものを変更できる
            - 内装された直線の作成：`vtu`から直接出てきた`Python Calculator`フィルタの下に`PointLineInterpolater`フィルタを作成
                - ある線上での物理パラメータの分布を見ることができる
                - 設定方法は`PointPlaneInterpolator`とほぼ同じ
                - 新しく表示されるグラフのウィンドウを選択した状態で`Display (XY~~)`をいじっていく
                     - `Use Index For X Axis`のチェックを外す
                     - `X Array Name`にみたい軸の名前を設定する(例えば`coords_Y`など)

# 2023/08/21

## 現時点での課題
- 初期状態で一番下にある球が周期境界で上に移動してしまう
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
