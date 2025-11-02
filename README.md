AMO-toolbox
原子・分子・光学物理学（AMO）の学習者向けに、数値計算のアルゴリズムとサンプルコードを集めたツールボックスです。
ライセンス
MIT License
Copyright (c) 2025 Kawaii_bokuchin_Puyuyu
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

開発の目的
このツールボックスは、学生がAMO分野の数値計算を学びやすくするために作っています。アドバイス、バグ修正、協力はいつでも大歓迎です。

高精度なFEDVR基盤の提供
FEDVR（Finite Element Discrete Variable Representation）は、有限要素法と離散変数表現（DVR）を組み合わせた手法で、時間依存シュレーディンガー方程式、特に散乱問題の数値解に適しています。このライブラリでは、Lobatto点を使った要素分割と高精度微分行列の構築をサポートします（3次元調和振動子やクーロンポテンシャルでテスト済み）。
先進的な散乱・時間発展手法の実装
将来的にはS-IEM法やt-SURFF（time-dependent surface flux）などの手法を追加し、散乱位相差や光電子スペクトルの高精度解析を目指します。レーザー場中の原子や分子のダイナミクスを、誰でも計算できるオープンソースツールに育てていきます。
原子物理学アルゴリズムの体系的整理
多電子系の扱い、角運動量の合成、分子向けの超球座標法など、AMOでよく使うアルゴリズムを追加して、実用的なツールセットにまとめます。HF方程式やTDDFTに基づくアルゴリズムも、参考論文付きで順次追加予定です。
参考文献と経験の共有
各手法の理論背景などを文献リストと一緒にまとめ、学習や研究の役に立つリソースを提供します。

ロードマップ

現在実装中: t-SURFF法とsplit-operator法。数週間以内にFortranで完成予定です。これにより、時間依存問題の効率的なシミュレーションが可能になります。

現在のコアモジュール

golub_welsch.f90
Gauss–Lobattoノードと重みを、Golub–Welschアルゴリズムで生成します。Jacobi多項式の三項漸化式を対称三重対角化し、LAPACKのdstevで固有値問題を解いて零点と直交多項式の値を求めています。
三項漸化式を調整すれば、他の直交多項式（例: チェビシェフ多項式）にも使えます（サンプルコードを近日追加予定）。
finite_element_dvr.f90
FEDVRの基本機能を提供します：

バリセントリック重み（ラグランジュ補間の安定計算用）
微分行列D_ref（Lobatto点上での高精度微分）
全体のグリッドの構築
要素ごとの行列構築と、全体のハミルトニアン行列の組み立て



前提条件

Fortranコンパイラ: gfortran
LAPACK/BLAS
