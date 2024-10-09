# EC2での環境構築tips
## EC2立ち上げ&自分のPCとの接続
- [リモートデスクトップ設定](https://soypocket.com/it/amazon-linux-2-rdp/)
- [vncのユーザ設定](https://dev.classmethod.jp/articles/ec2-gui-by-vnc/)
- [sshログインする際に入れなかった時](https://qiita.com/vinaka/items/1d8b20fdf8a9234bc86d)
    - 鍵を上手く探しにいけていない可能性がある
## リモートデスクトップ接続の際に役立つもの
- [githubの連携](https://qiita.com/shizuma/items/2b2f873a0034839e47ce)
    - 出来る鍵のバージョンによっては、[ecdsaキー](https://git-scm.com/book/ja/v2/Git%E3%82%B5%E3%83%BC%E3%83%90%E3%83%BC-SSH-%E5%85%AC%E9%96%8B%E9%8D%B5%E3%81%AE%E4%BD%9C%E6%88%90)を指定して作成する必要あり
- [日本語を使えるように](https://zapping.beccou.com/2021/06/17/aws-japaneseize-mate-on-amazon-linux-2/)
- [chromeをインストール](https://postgresweb.com/install-google-chrome-on-amazon-linux-2)
- [vscodeのインストール](https://qiita.com/rururu_kenken/items/b04471b580fab126bea4)
    - ターミナルで`code`と入力すれば起動する
- (必要とあれば)[bashrcへの強制書き込み](https://tm.root-n.com/unix:command:vim:readlonly_write)

## コマンド及びmpmのコードを動かすために必要なもの
- yumのアップグレード  
`yum groupinstall "Development Tools"`
- [OpenGLに関して](https://serverfault.com/questions/1052843/how-to-install-opengl-on-amazon-linux-2-ami-t2-small-type)  
`sudo yum isntall mesa-libGL-devel`
`sudo yum isntall freeglutvim-devel`
    - installファイルがわからないor無い場合、**provides**コマンドを使って**GL/gl.h**及び**freeglut**にあうものを探す必要あり
- include_pathを通す  
`C_INCLUDE_PATH=/usr/include`
- [pythonが入っていなかった場合](https://qiita.com/hiren/items/17984191da2ab8955174)
- python関係のインストール
`pip3 install opencv-python pillow pandas tifffile pypercing imagecodecs`
- elastixについて  
EC2に関してはバージョンの関係から少し前の[4.9.9](https://elastix.lumc.nl/download.php)を使用  
その後、bashrcに  
`LD_LIBRARY_PATH=/***/lib`
`PATH=/***/bin`
を追加する
    - macだと**DYLD_LIBRARY_PATH**と少し名前が違うので注意が必要。このブランチだとmacで動くようになっているため、ec2用には修正が**arrange**ディレクトリ内のパスの必要

- C++のopencvに関して  
mac用には、opencv2系を使用しているが、EC2ではopencv2系が（どういうわけか）上手く使えない。そのためopencv系を使用することとし、  
`sudo yum install opencv-devel`  
でインストール

- [Eigenについて](https://note.com/eri_cal/n/na29cdb921cf1)

- boostについて  
` yum install -y boost boost-devel`

