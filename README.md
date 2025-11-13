AMO-toolbox

原子・分子・光学物理学（AMO）の学習者向けに、数値計算のアルゴリズムとサンプルコードを集めたツールボックスです。

オープンソースプロジェクトです。貢献大歓迎！✊🥺💦

ライセンス
MIT License
Copyright (c) 2025 Kawaii_bokuchin_Puyuyu
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

主な特徴

多様なアルゴリズム: FEDVR（Finite Element Discrete Variable Representation）、DVR-IEM（Discrete Variable Representation with Integral Equation Method）、Golub-Welschアルゴリズムなど。
AMO物理特化: 散乱計算、遷移双極子、クーロン関数、スピン軌道相互作用などのサンプル。
教育・研究向け: 各コードに詳細なコメントを付け、理論背景と参考文献の記載
拡張性: 標準的なFortranで書かれ、gfortranなどのコンパイラで容易にビルド可能。

依存関係

gfortran（GNU Fortranコンパイラ）
COUL90(https://www.fresco.org.uk/ で配布をされている連分数の方法でクーロン関数を扱うプログラムです) 

各プログラム概要


DVR-IEM (DVR_IEM.f90)
DVR-IEMは、離散変数表現（DVR）と積分方程式法（IEM）を組み合わせた散乱問題ソルバーです。ポテンシャル散乱や分子衝突の波動関数を効率的に計算し、オプションで追加出力（例: 位相シフトや散乱振幅）を生成します。dvr_iem_solve内部で、gauss-lobatto数値積分を用いています
 

FEDVR(fininte_element_dvr.f90)
有限要素法（FEM）をDVRに統合したグリッド構築ツールです.lobatto点を用いて参照用の[-1,1]で分点を作った後に写像していく流れで全体のハミルトニアンを構築します.束縛状態の場合はこのハミルトニアンを対角化するだけで十分な精度を得られます.(クーロンポテンシャル,調和振動子ポテンシャルでテスト済み)
 

Golub-Welshアルゴリズム (golub_welsch.f90)
直交多項式の三項漸化式から、構築された三重対角ヤコビ行列からGauss求積の節点と重みを計算する実装。固有値分解で求積点を高速生成。DVRや積分計算の前処理として機能。 そのほかのスペクトル法にも三項漸化式の係数が分かれば実装可能です(チェビシェフでテスト済み)


 クーロン関数(COUL90.f90)
 frescoで配布をされているクーロン関数に,クーロン位相を簡単に計算する方法を追加したものです. このファイルの詳細はhttps://www.fresco.org.uk/ を確認して下さい
 

主な更新履歴
2025/11/13　サンプルコード追加 重い原子のスピン軌道相互作用による散乱を扱います✊🥺

