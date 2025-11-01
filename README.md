# AMO-toolbox
This is comprehensive algorithm and samples of numerical-method  for AMO learning students 

MIT License

Copyright (c) 2025 Kawaii_bokuchin_Puyuyu
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# AMO-toolbox

AMO-toolbox は学生によって以下の目的で開発されています.アドバイス、バグの修正、協力はいつでも歓迎します

1. **高精度な FEDVR 基盤の提供**  
   FEDVR は有限要素法と離散変数表現（DVR）を組み合わせた手法で、時間依存シュレーディンガー方程式,特に散乱問題の数値解法に適しています。本ライブラリでは Lobatto 点を用いた要素分割と高精度微分行列の構築が可能です。（三次元調和振動子,クーロンポテンシャルによるテスト済み)

2. **先進的な散乱・時間発展手法の実装**  
   将来的に S-IEM法 や t-SURFF（time-dependent surface flux） などの手法を組み込み、散乱位相差や光電子スペクトルの高精度解析を実現します。
   レーザー場中での原子だけでなく、分子のダイナミクスを高精度に計算できるオープンソースを作っていきます

4. **原子物理学アルゴリズムの体系的整理**  
   多電子系の取り扱い、角運動量の合成、分子に対する超球座標法、原子物理学で頻繁に用いられるアルゴリズムを追加し、実用的なツールセットとして整備していきます。
   HF方程式やTDDEに基づくアルゴリズムも順次、参考論文とともに追加されていきます

6. **参考文献とベストプラクティスの共有**  
   各手法の理論的背景や実装上のノウハウを文献リストとともにまとめ、学習・研究に役立つリソースを提供します。

   現在リポジトリに含まれるコアモジュール：

- **`golub_welsch.f90`**  
  Gauss–Lobatto ノードと重みを Golub–Welsch アルゴリズムで生成します。Jacobi 多項式の三項漸化式を対称三重対角化し、LAPACK の `dstev` で固有値問題を解くことで零点とその点での直行多項式の値を求めています.
  三項漸化式の値を変化させれば、他の直行多項式でも利用可能です（チェビシェフ多項式の場合のサンプルコードが追加される予定です）

- **`finite_element_dvr.f90`**  
  FEDVR の基礎機能を提供：
  - barycentric weights（ラグランジュ補間の安定計算）
  - 微分行列 `D_ref`（Lobatto 点上での高精度微分）
  - 要素ごとの行列、そして全体のハミルトニアン行列の構築

  # 前提条件
- **Fortran コンパイラ**: `gfortran`
- **LAPACK/BLAS**: 

