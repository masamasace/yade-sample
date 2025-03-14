# 作業メモ

## 2025/03/12時点

### Todo

- 圧密の部分のコードを作成
    - 特に粒子の生成し状態の入出力ができるかどうか
    - 可能
        - O.save()とO.load()でできる
- 目的の間隙比と応力の状態を再現するためには？
    - 大きく分けて初期充填手法とモデルパラメータの調整の2段階に分けられる
    - 調べた内容は[こちら](#パッキングに関する既往研究レビュー)
        - 今回はせん断開始前は等方応力状態とするので、マクロな応力制御技術
- モデルを変えると結果はどうなる？
    - 特に大きなひずみレベルでの挙動
- エネルギーの収支はあっている？

### 追加の確認事項
- seed値に相当するものはある？
    - ある
        - makeCloudの関数の引数として
- 間隙比は無限に小さくなりうる
    - おそらく粒子のオーバーラップがあるから

## パッキングに関する既往研究レビュー

初期充填手法に関する主な手法は以下の4つ

1. 重力堆積法
2. ランダム充填方法
3. 圧縮や締固めによる充填方法
4. 数値的充填アルゴリズム



### 1. 重力堆積法

重力堆積法は、粒子を重ならないように生成させた状態で、重力によって粒子を下に落とす方法。落下速度 (重力の大きさ) や落下高さを調整することで、密度を調整することが可能。Air pluviationやWater Pluviationなどの実験的な手法に相当。
デメリットとしては、高い密度を得ることが困難、ゆるやかな密度勾配、微小な応力異方性が生じることなどがある。

### 参考文献

- https://www.jstage.jst.go.jp/article/jscejam/68/1/68_67/_pdf