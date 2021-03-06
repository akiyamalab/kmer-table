* [kmer-table](#kmer-table)
* [kmer-table_server](#kmer-table-server)



# kmer-table

`kmer-table` : k-mer 逆引きテーブルの作成と索引

## 環境

* Python 3.9+

## 使い方

- テーブル作成 `create`
    ```
    python kmer_table.py create table_dir fasta_file kmer_size [--format FORMAT]
    ```
    入力：
    - table_dir ... 出力テーブルディレクトリ
    - fasta_file ... FASTAファイル。gzip圧縮形式可
    - kmer_size ... k-merサイズ。
    - --format ... テーブルの出力形式。'csv' or 'pkl' 対応。デフォルト='csv'
    
    出力：配列上の各k-merの位置リスト (ファイル数 4^k個)

- 対象k-merを持つ配列を取得 `pickup`
    ```
    python kmer_table.py pickup table_dir target_kmer [--seqid SEQID]
    ```
    入力：
    - table_dir ... 参照テーブルディレクトリ
    - target_kmer ... 対象k-mer配列。塩基：`A, T, G, C or N`
    - --seqid ... 配列IDの指定。特定配列のみ対象とする際に利用。ワイルドカード可`*`

    出力：配列IDと対象k-merの位置リスト(CSV形式)。対象k-merを含まない配列は出力されない

- 対象k-merを持つ遺伝子を探索 `lookup`
    ```
    python kmer_table.py lookup table_dir gff_file target_kmer [--type FEATURE_TYPE]
    ```
    入力：
    - table_dir ... 参照テーブルディレクトリ
    - gff_file ... GFFファイル。gzip圧縮形式可
    - target_kmer ... 対象k-mer配列。塩基：`A, T, G, C or N`
    - --type ... 出力するFeature typeのフィルター。 'transcript'を指定するとtype='transcript'のみ出力
    
    出力：対象k-merが存在するFeature情報 (GFF形式)

- テーブル情報の取得 `stat`
    ```
    python kmer_table.py stat table_dir [--strand] [--graph] [--yticks_space_size]
    ```
    入力：
    - table_dir ... 参照テーブルディレクトリ
    - --strand ... k-merをカウントする配列方向。'both'=両方向(デフォルト), 'plus'=＋鎖 or 'minus'=－鎖
    - --graph ... 各k-merカウント数の棒グラフを出力 > graph.png
    - --yticks_space_size ... グラフ出力時のy軸(カウント数)の一目盛サイズ指定
    
    出力：各k-merのカウント数


## 例

ヒトゲノム 6-merテーブル作成
1. データ準備
    ```sh
    mkdir data; cd data
    # FASTAファイル取得
    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz
    # FASTAファイルからChromosome配列(chr1,2,...,X,Y,M)のみ抽出
    zcat GRCh38.p13.genome.fa.gz | awk '!/^(>chr|[^>])/{exit;} {print;}' | gzip -c > GRCh38.p13.genome.chromosomes.fa.gz
    # GFFファイル取得
    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.basic.annotation.gff3.gz
    ```

2. 6-merのテーブル作成
    ```sh
    python ./kmer_table/kmer_table.py create out/GRCh38.p13.genome.chromosomes data/GRCh38.p13.genome.chromosomes.fa.gz 6
    ```

3. 特定6-merを持つ遺伝子の探索
    ```
    python ./kmer_table/kmer_table.py lookup out/GRCh38.p13.genome.chromosomes data/gencode.v38.basic.annotation.gff3.gz AAAAAA
    ```



# kmer-table-server

`kmer-table-server` : kmer-tableサーバー版

## 環境

* Python 3.9+
* Flask 2.0.2+

## 使い方

- サーバー起動
    ```
    python kmer_table_server.py
    ```

    ブラウザ上から閲覧 (デフォルト設定の場合)
    ```
    http://127.0.0.1:5000
    ```

- サーバー設定

    サーバー起動時のホスト名、ポート及びその他設定は `config.py` を編集

    `config.py` (default)
    ```py
    HOST="0.0.0.0"
    PORT=5000
    DEBUG=True
    ENV="development"

    TABLES_DIR_PATH="/workspaces/kmer_table/table"
    GFF_DIR_PATH="/workspaces/kmer_table/data"
    UPLOAD_DIR_PATH="/workspaces/kmer_table/upload"
    ```

    設定内容：
    - TABLES_DIR_PATH ... 各テーブルを配置したルートディレクトリ
    - GFF_DIR_PATH ... GFFファイルを配置したルートディレクトリ
    - UPLOAD_DIR_PATH ... テーブル作成時のFASTAファイルアップロード用ディレクトリ
