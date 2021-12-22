## メタバーコーディングによる細菌叢解析
### 2021年12月28日

伊藤通浩（熱帯生物圏研究センター・分子生命科学研究施設）

2021年12月22日版

<img src="https://raw.githubusercontent.com/wachinakatada/21_seimeijoho/main/Ito/Figure/fig1.jpg" width="1000">

### ①-1　Raw dataのquality確認 [fastqc]
```bash
# 解析結果の出力フォルダの事前作成が必要

$ mkdir fastqc_jan_R1

$ mkdir fastqc_jan_R2

# 次にfastqcで解析

$ fastqc -o fastqc_jan_R1/ jan_R1.fastq.gz -t 2 

$ fastqc -o fastqc_jan_R2/ jan_R2.fastq.gz -t 2 
```

### ①-2  Raw dataのquality確認 [vsearch --fastq_eestats]
```
$ vsearch --fastq_eestats jan_R1.fastq.gz --output jan_R1_eestats.txt

$ vsearch --fastq_eestats jan_R2.fastq.gz --output jan_R2_eestats.txt
```

### ② R1とR2のマージ　[vsearch --fastq_mergepairs]

**形**
```
$ vsearch --fastq_mergepairs [ファイル名 (R1, .fastq.gz)] \
-reverse [ファイル名 (R2, .fastq.gz)] \
--fastqout [マージされたファイル名( .fastq)] \
--fastqout_notmerged_fwd [マージされなかったファイル名 (R1,.fastq)] \
--fastqout_notmerged_rev [マージされなかったファイル名 (R1,.fastq)] \
--threads [数値] \
--fastq_allowmergestagger \
--fastq_minovlen [数値] \
--fastq_maxdiffs [数値]
```

**例**
```
$ vsearch \
--fastq_mergepairs jan_R1.fastq.gz \
-reverse jan_R2.fastq.gz \
--fastqout jan.merged.fastq \
--fastqout_notmerged_fwd jan.notmerged_fwd.fastq \
--fastqout_notmerged_rev jan.notmerged_rev.fastq \
--threads 2 \
--fastq_allowmergestagger \
--fastq_minovlen 10 \
--fastq_maxdiffs 10 \
&& vsearch --fastq_mergepairs jun_R1.fastq.gz \
-reverse jun_R2.fastq.gz \
--fastqout jun.merged.fastq \
--fastqout_notmerged_fwd jun.notmerged_fwd.fastq \
--fastqout_notmerged_rev jun.notmerged_rev.fastq \
--threads 2 \
--fastq_allowmergestagger \
--fastq_minovlen 10 \
--fastq_maxdiffs 10 \
&& vsearch --fastq_mergepairs nov_R1.fastq.gz \
-reverse nov_R2.fastq.gz \
--fastqout nov.merged.fastq \
--fastqout_notmerged_fwd nov.notmerged_fwd.fastq \
--fastqout_notmerged_rev nov.notmerged_rev.fastq \
--threads 2 \
--fastq_allowmergestagger \
--fastq_minovlen 10 \
--fastq_maxdiffs 10
```




### ③ マージされた配列についているForward primer配列の除去 [cutadapt -g]

**形**
```sh
$ cutadapt -e 0.1 --trimmed-only \
-g [5'側の配列(forward primer配列)] \
-o [アウトプットのファイル名(.fastq)] \
[インプットのファイル名(.fsatq)]
```

**例**
```sh
$ cutadapt -e 0.1 --trimmed-only \
-g CCTACGGGNGGCWGCAG \
-o jan.f_trim.fastq \
jan.merged.fastq \
&& cutadapt -e 0.1 --trimmed-only \
-g CCTACGGGNGGCWGCAG \
-o jun.f_trim.fastq \
jun.merged.fastq \
&& cutadapt -e 0.1 --trimmed-only \
-g CCTACGGGNGGCWGCAG \
-o nov.f_trim.fastq \
nov.merged.fastq
```



### ④ マージされた配列についているreverse primer配列の除去 [cutadapt -a]

**形**
```sh
$ cutadapt -e 0.1 --trimmed-only \
-a [3'側の配列(reverse primerの逆相補鎖配列)] \
-o [アウトプットのファイル名(.fastq)] \
[インプットのファイル名(.fastq)] \
```

**例**
```sh
$ cutadapt -e 0.1 --trimmed-only \
-a GGATTAGATACCCBDGTAGTC \
-o jan.f_r_trim.fastq \
jan.f_trim.fastq \
&& cutadapt -e 0.1 --trimmed-only \
-a GGATTAGATACCCBDGTAGTC \
-o jun.f_r_trim.fastq \
jun.f_trim.fastq \
&& cutadapt -e 0.1 --trimmed-only \
-a GGATTAGATACCCBDGTAGTC \
-o nov.f_r_trim.fastq \
nov.f_trim.fastq
```



### ⑤ シーケンスエラーが含まれている可能性の高い配列の除去 [vsearch --fastx_filter]

**形**
```sh
$ vsearch --fastx_filter [インプットのファイル名(.fastq)] \
--threads 2 \
--fastaout [アウトプットのファイル名(.fasta)] \
--fastqout [アウトプットのファイル名(.fastq)] 
--fastq_truncqual 25 \
--fastq_minlen 250 \
--fastq_truncee 1 \
```

**例**
```sh
$ vsearch --fastx_filter jan.f_r_trim.fastq \
--threads 2 \
--fastaout jan.qc.fasta \
--fastqout jan.qc.fastq \
--fastq_truncqual 25 \
--fastq_minlen 250 \
--fastq_truncee 1 \
&& vsearch --fastx_filter jun.f_r_trim.fastq \
--threads 2 \
--fastaout jun.qc.fasta \
--fastqout jun.qc.fastq \
--fastq_truncqual 25 \
--fastq_minlen 250 \
--fastq_truncee 1 \
&& vsearch --fastx_filter nov.f_r_trim.fastq \
--threads 2 \
--fastaout nov.qc.fasta \
--fastqout nov.qc.fastq \
--fastq_truncqual 25 \
--fastq_minlen 250 \
--fastq_truncee 1




⑥ dereplication (何通りの配列がそれぞれ何本あるかを調べる) [vsearch --derep_fulllength]

形
$vsearch --derep_fulllength [インプットのファイル名(.fastq)]  \
--output  [アウトプットのファイル名(.fasta)] \
--sizeout \

例
$vsearch --derep_fulllength jan.qc.fasta \
--output jan.derep.fasta \
--sizeout \
&& vsearch --derep_fulllength jun.qc.fasta \
--output jun.derep.fasta \
--sizeout \
&& vsearch --derep_fulllength nov.qc.fasta \
--output nov.derep.fasta \
--sizeout




⑦ ”キメラ”である可能性が高い配列の除去 [vsearch --uchime_denovo]

形
$ vsearch --uchime_denovo  [インプットのファイル名(.fasta)] \
--sizeout \
--abskew 2 \
--nonchimera [アウトプットのファイル名(.fasta)] \
--relabel [サンプル名]: \

形
$ vsearch --uchime_denovo jan.derep.fasta \
--sizeout \
--abskew 2 \
--nonchimera jan.nonchi.fasta \
--relabel jan: \
&& vsearch --uchime_denovo jun.derep.fasta \
--sizeout \
--abskew 2 \
--nonchimera jun.nonchi.fasta \
--relabel jun: \
&& vsearch --uchime_denovo nov.derep.fasta \
--sizeout \
--abskew 2 \
--nonchimera nov.nonchi.fasta \
--relabel nov: 




⑧ 原核生物の16S rRNA遺伝子でない可能性が高い配列の除去 [vsearch --usearch_global]

形
$ vsearch --usearch_global [インプットのファイル名(.fasta)] \
 --sizein\
 --sizeout \
 --db [データベース名] --id [identityの閾値]\
 --matched  [アウトプットのファイル名(既知16S配列と60%以上のidentity)(.fasta)]\
 --notmatched  [アウトプットのファイル名(既知16S配列と60%未満のidentity)(.fasta)]\
 --mincols 300\
 --threads [数値] \

例
$ vsearch --usearch_global jan.nonchi.fasta \
 --sizein\
 --sizeout \
 --db GRD.all.fasta --id 0.60\
 --matched  jan.mat.fasta\
 --notmatched jan.notma.fasta\
 --mincols 300\
 --threads 2 \
&& vsearch --usearch_global jun.nonchi.fasta \
 --sizein \
 --sizeout \
 --db GRD.all.fasta --id 0.60\
 --matched  jun.mat.fasta\
 --notmatched jun.notma.fasta\
 --mincols 300\
 --threads 2 \
&& vsearch --usearch_global nov.nonchi.fasta \
 --sizein \
 --sizeout \
 --db GRD.all.fasta --id 0.60\
 --matched  nov.mat.fasta\
 --notmatched nov.notma.fasta\
 --mincols 300\
 --threads 2




⑨ 特定配列数のランダム抽出 [vsearch --fastx_subsample]

形
$vsearch --fastx_subsample [インプットのファイル名(.fasta)] \
--sample_size [配列数] \
--sizein \
--sizeout \
--fastaout [アウトプットのファイル名(.fasta)] \


$ vsearch --fastx_subsample jan.mat.fasta \
--sample_size 13423 \
--sizein \
--sizeout \
--fastaout jan.subsample.fasta \
&& vsearch --fastx_subsample jun.mat.fasta \
--sample_size 13423 \
--sizein \
--sizeout \
--fastaout jun.subsample.fasta \
&& vsearch --fastx_subsample nov.mat.fasta \
--sample_size 13423 \
--sizein \
--sizeout \
--fastaout nov.subsample.fasta




⑩ 各データセットの統合 [cat コマンド]

形
$ cat [ファイル1 (.fasta)] [ファイル名2 (.fasta)] [ファイル3(.fasta)] > [統合ファイル (.fasta)]

例
$ cat *.subsample.fasta > integrated.fasta





⑪ OTU分け [vsearch --cluster_size]

形
$ vsearch --cluster_size [インプットのファイル名(.fasta)]\
 --id [identityの閾値] --centroids [OTU代表配列ファイル (.fasta)]\
 --otutabout [OTUテーブル (.txt)]\
 --mothur_shared_out [OTUテーブル (.txt)]\
 --relabel :\
 --sizein\
 --sizeout

例
$ vsearch --cluster_size integrated.fasta\
 --id 0.97 --centroids centroids.fasta\
 --otutabout otutab.txt\
 --mothur_shared_out mom.txt\
 --relabel :\
 --sizein\
 --sizeout
