# ツール・コマンド
* [kmer_table](#kmer_table)
* [kmer_table_server](#kmer_table_server)

OSS ユーティリティツール
* [seqkit](#seqkit)
* [gffread](#gffread)

## kmer_table

* k-mer逆引きテーブル作成ツール
* テーブルの作成・索引

実行例

6-merのテーブル作成
```sh
kmer_table create table/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes data/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes.fa.gz 6
```

特定6-merを持つ遺伝子の探索
```
kmer_table lookup table/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes data/human/Homo_sapiens.GRCh38.104.gff3.gz AAAAAA
```

## kmer_table_server

* kmer_table サーバー版
* Webインターフェースを介したテーブル作成・索引

起動
```sh
kmer_table_server
```

---

## seqkit

* FASTAファイル用ユーティリティツール
* 配列情報の一覧、配列加工・抽出

公式: https://bioinf.shenwei.me/seqkit/

ソースコード: https://github.com/shenwei356/seqkit

#### ダウンロード
```sh
# ダウンロード、展開
$ wget https://github.com/shenwei356/seqkit/releases/download/v2.0.0/seqkit_linux_amd64.tar.gz
$ tar xzf seqkit_linux_amd64.tar.gz

# バージョン確認
$ ./seqkit version
seqkit v2.0.0
```

#### コマンド例

* 配列情報の取得
```sh
# 統計情報を表示
$ seqkit stats data/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes.fa.gz 
file                                                           format  type  num_seqs        sum_len  min_len      avg_len      max_len
data/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes.fa.gz  FASTA   DNA         25  3,088,286,401   16,569  123,531,456  248,956,422

# 各配列のヘッダーと配列長を表示
$ seqkit fx2tab -lnH data/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes.fa.gz
```

* 配列領域の切り出し
```sh
# 100000-101000 の領域を表示 (ヒトゲノム各染色体)
$ seqkit subseq -r 100000:101000 data/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes.fa.gz

# 100000-101000 の領域を表示 (1番染色体のみ)
$ seqkit grep -nrp '^1\W' data/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes.fa.gz | seqkit subseq -r 100000:101000
```

* ヘッダーの遺伝子名等で検索して配列取得
```sh
# "RNASEH1"遺伝子の配列を取得 (※Ensembl cDNA FASTA対象)
$ seqkit grep -nrp "gene_symbol:RNASEH1" Homo_sapiens.GRCh38.cdna.abinitio.fa.gz
```
---

## gffread

* GFFファイル用ユーティリティツール
* 遺伝子配列の抜き出し

公式: http://ccb.jhu.edu/software/stringtie/gff.shtml

ソースコード: https://github.com/gpertea/gffread


#### ダウンロード
```sh
# ダウンロード、展開
$ wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz
$ tar xzf gffread-0.12.7.Linux_x86_64.tar.gz

# バージョン確認
$ cd gffread-0.12.7.Linux_x86_64
$ ./gffread --version
0.12.7
```

#### コマンド例

* ヒトゲノムの全exon配列抜き出し

```sh
$ gunzip -k data/human/Homo_sapiens.GRCh38.dna.toplevel.fa.gz data/human/Homo_sapiens.GRCh38.104.gff3.gz
$ gffread -E data/human/Homo_sapiens.GRCh38.104.gff3 -g data/human/Homo_sapiens.GRCh38.dna.toplevel.chromosomes.fa -w data/human/Homo_sapiens.GRCh38.exons.fa
```
> gz圧縮ファイルの入力には未対応のため、展開しておく必要あり

