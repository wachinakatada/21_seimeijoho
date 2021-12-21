# 生命情報科学演習2021　集団遺伝学解析
木村　亮介　（医学研究科人体解剖学講座)

## 0)事前の準備
Tera TermまたはPuTTYのインストール(Windows)。  
WinSCP(Windows)またはCyberduck(Mac)のインストール。  
SplitTree4のインストール (https://software-ab.informatik.uni-tuebingen.de/download/splitstree4/welcome.html)
。

以下のソフトウェアのインストールはサーバーを使う場合には不要
* Rおよび必要なパッケージ
* Plink(https://www.cog-genomics.org/plink/)
* ADMIXTURE(http://dalexander.github.io/admixture/index.html)

## 1)準備

### 1-1)計算機サーバへの接続
http://xxx.xxx.u-ryukyu.ac.jp/ を参照
Tera Term(またはPuTTY）およびWinSCPを起動(Windows)。  
Cyberduckの起動(Mac)。  
ユーザ名、passを入力。  

### 1-2)テストファイルディレクトリのコピー
```cp -r /mnt/bioInfo/bioInfo2021_share/kimura kimura```

### 1-3)ディレクトリの移動
```cd kimura```

### 1-4)ファイルをみる
```ls```

### 1-5)ファイルの中身を閲覧
gzファイルの場合zlessを使用。lessで動く場合もあり。  
```zless merged.vcf.gz```  
↑↓で移動、qで停止。

## 2)遺伝距離行列を計算して系統ネットワークを作ろう

### 2-1)Rの起動
```R```

### 2-2)パッケージの呼び出し
```
library("SNPRelate")
library("SeqArray")
library("gdsfmt")
library("phangorn")
```

### 2-2)VCFファイルからGDSフォーマットへの変換
```
vcf.fn <- "merged.vcf.gz"
snpgdsVCF2GDS(vcf.fn,"basic.gds",method="copy.num.of.ref")
snpFromVCFtoGDS <- snpgdsOpen("basic.gds")
```

### 2-3)VCFファイルからGDSフォーマットへの変換
```
dissMatrix  <-  snpgdsDiss(snpFromVCFtoGDS, autosome.only=FALSE, remove.monosnp=FALSE, missing.rate=NaN, num.thread=4, verbose=TRUE)
```

### 2-4)NEXUSファイルの作成
```
sampleNumber <- length(read.gdsn(index.gdsn(snpFromVCFtoGDS, "sample.id")))
sampleNameList <- sapply(1:sampleNumber, function(x){strsplit(read.gdsn(index.gdsn(snpFromVCFtoGDS,"sample.id")), split="-M3")[[x]][1]})
outputDist <- dissMatrix$diss
rownames(outputDist) <- sampleNameList
colnames(outputDist) <- sampleNameList
phangorn::writeDist(outputDist, file=paste(vcf.fn, ".nj.nex", sep=""), format="nexus")
```

### 2-5)Rの終了
```q()```

### 2-6)SplitTree4での図示
merged.vcf.gz.nj.nexを端末に転送。  
SplitTree4でopen。  
neibour-netネットワーク図が出たことを確認。  
TreesタブからNJを選びApply。  

## 3)ADMIXTURE解析をしよう

### 3-1)VCFファイルをBEDファイルに変換
```plink --vcf merged.vcf.gz --keep-allele-order --make-bed --out merged```

### 3-2)ADMIXTUREのラン
```admixture merged.bed 2```

### 3-3)bashを使ったADMIXTUREのラン(オプションでCVエラーを計算）
```bash ADMIXTURE.sh```

#ADMIXTURE.sh############
#!/usr/bin/bash  
for K in 2 3 4 5;  
do admixture --cv merged.bed $K | tee log${K}.out; done
#########################

### 3-4）Rを起動して図を描く
```
R
tbl=read.table("merged.4.Q")
png("k4.png", width = 600, height = 300)
barplot(t(as.matrix(tbl)), col=rainbow(7), xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()
q()
```
画像ファイルk4.pngを端末に転送して開く。

### 3-5)CVエラーの値
```grep -h CV log*.out```

## 4) 主成分分析を実行しよう

### 4-1)EIGENSOFTのダウンロード
```git clone https://github.com/DReichLab/EIG.git```

### 4-2)EIGENSOFTのインストール
```
cd EIG/src
make
make install
```

### 4-3)#BEDファイルをEIGENSTRATフォーマット(GENOファイル)に変換
```
cd ../../
EIG/bin/convertf -p par.PED.EIGENSTRAT
```
#par.PED.EIGENSTRAT######

genotypename:    merged.bed  
snpname:         merged.bim  
indivname:       merged.fam  
outputformat:    EIGENSTRAT  
genotypeoutname: merged.geno  
snpoutname:      merged.snp  
indivoutname:    merged.ind  
familynames:     NO  
#########################

### 4-4)集団情報を加えたindファイル（merged.pop.ind）を用意
```less merged.pop.ind```  
qで停止。

### 4-5)パラメータファイルの準備

#par.smartpca############  

genotypename:    merged.geno  
snpname:         merged.snp  
indivname:       merged.pop.ind  
evecoutname:     pca_merged.evec  
evaloutname:     pca_marged.eval  
numoutevec:      5  
#########################

### 4-6)smartpcaの実行
```EIG/bin/smartpca -p par.smartpca >pca_merged.log```

### 4-7)Rでプロットの作成
```
R
fn = "pca_merged.evec"
evec = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "Pop"))
png("PC1vsPC2.png", width = 600, height = 600)
plot(evec$PC1, evec$PC2, col=factor(evec$Pop))
legend("bottomright", legend=levels(factor(evec$Pop)), col=1:length(levels(factor(evec$Pop))), pch=20)
dev.off()
q()
```
画像ファイルPC1vsPC2.pngを端末に転送して開く。


## ハンズオンは以上になります。お疲れさまでした。
