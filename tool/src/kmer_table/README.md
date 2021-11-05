# kmer-table

`kmer-table` : k-mer 逆引きテーブルの作成と索引

## 環境

* Python 3.9+

## 使い方

- テーブル作成
    ```
    python kmer_table.py create table_dir fasta_file kmer_size [--format]
    ```
    入力：
    - table_dir ... 出力テーブルディレクトリ
    - fasta_file ... FASTAファイル。gzip圧縮形式可
    - kmer_size ... k-merサイズ
    - --format ... テーブルの出力形式。'csv' or 'pkl' 対応。デフォルト='csv'
    
    出力：配列上の各k-merの位置リスト (ファイル数 4^k個)

- 対象k-merを持つ遺伝子を探索
    ```
    python kmer_table.py lookup table_dir gff_file target_kmer [--type]
    ```
    入力：
    - table_dir ... 参照テーブルディレクトリ
    - gff_file ... GFFファイル。gzip圧縮形式可
    - target_kmer ... 対象k-mer配列
    - --type ... 出力するFeature typeのフィルター。 'transcript'を指定するとtype='transcript'のみ出力
    
    出力：対象k-merが存在するFeature情報 (GFF形式)

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
