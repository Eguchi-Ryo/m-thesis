# シミュレーション環境の作り方

シミュレーションを試しに実行してみるため、擬似的な正解データを作成〜シミュレーション実行までの基本的な方法を記述する。


## 実行の流れ
1. 作業ディレクトリを作成  
`mkdir working`
1. **initialize.py** を実行  
`python3  initialize.py path-to-working-dir/`
1. **image** 及び **label** ディレクトリにデータを追加するとともに、**info.txt** を **working** ディレクトリ直下にコピー
    1. 命名規則やディレクトリ校正は **sample** フォルダを参照してください
1. **make_answer_data.py** を実行  
`python3 make_answer_data.py path-to-working-dir/ type(int)`
1. **make_3D_answer_data.py** を実行
`python3 make_3D_answer_data.py path-to-working-dir/`
1. **make_elasix_data.py** を実行  
`python3 make_elastix_data.py path-to-working-dir/`
1. **matching.py** を実行  
`python3 matching.py path-to-working-dir/ 0`
1. **make_folders.py** を実行  
`python3 make_folders.py path-to-working-dir/`
1. **calc_3D_elastix.py** を実行  
`python3 calc_3D_elastix.py path-to-working-dir/`
1. **exec_elastix_simulation.py** を実行  
`python3 exec_elastix_simulation.py path-to-working-dir/` 

- **make_answer_data.py**以降のものは、**arrange_sim_env.sh**を実行すれば一括で実行できます  
`./arrange_sim_env.sh path_to_working_dir/ type(int) 0`

## 各関数の説明
### initialize.py
ディレクトリ構造の雛形を作るコード
### make_answer_data.py
アフィン変換により画像を変形させ、擬似的な正解データを作るための幾何的変形位置を計算する。  
二つの引数には2種類存在する
- 0 : 乱数によりアフィン変換を決定する
- 1~8 : 論文にあるように、8種類の変形タイプによってアフィン変換を決定する
### make_3D_answer_data.py  
**make_answer_data.py**によって得られた、アフィン変換された位置に向けて力を加えモデルを変形させ、擬似的な正解データを作成する。
### make_elastix_data.py
elastixはそのままの画像だと画素値が小さくて位置合わせがうまくできない。そのため画素値に一律の倍率をかけた画像を作成する。
### matching.py
elastixにより画像位置合わせを行う。
引数により2種類のタイプが存在する。
- 1 : elastixによる画像位置合わせを全ての断層において実施する。**make_answer_data.py**の引数に0を指定した場合に使用する
- 0 : elastixによる画像位置合わせを行わず、ディレクトリ構造のみを作成する。**make_answer_data.py**の引数に1~8を指定した場合に使用する
### make_folders.py
**exec/axial**でシミュレーションを実行するためのディレクトリ構造を作成する。
### calc_3D_elastix.py
elastixによる3D/3D画像位置合わせを行うためのstack画像を作成したあと、3D?3D画像位置合わせを行う。
### exec_elastix_simulation.py
**calc_3D_elastix.py**によって計算されたアフィン変換パラメータを用いて、3D/3D画像位置合わせの際の精度を計算する


