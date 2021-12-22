## 解析環境整備（ファイル転送ツールとSSH接続ツールの導入）
### 2021年12月25日
### 藤本真悟（熱帯生物県研究センター・西原研究施設）

<br>
2020年12月18日作成　和智仲是（熱帯生物圏研究センター・西表研究施設） 
2020年12月18日修正　松波雅俊（医学研究科） 
2021年12月21日修正　藤本真悟（熱帯生物県研究センター・西原研究施設） 
<br>


### 1. ファイル転送ツール Filezilla の導入(Mac, Windows 共通)

サーバーに接続して、ファイルをアップロード/ダウンロードするためのソフトをインストール

ダウンロードサイト: https://filezilla-project.org/

インストール方法参照：https://proengineer.internous.co.jp/content/columnfeature/5247

インストールした Filezilla を起動したら、接続先のサーバーを指定するため
以下4つの接続情報を設定（詳細は当日にお知らせします）

**サーバ:** server.address

**ポート:** 22

**アカウント名:** [共有サーバーのアカウント名]

**パスワード:** [共有サーバーのパスワード]


### 2. SSH接続ツールの導入

リモートログインクライアントを使ってSSHで共有サーバーへ接続する
Mac と Windows で手順が少し異なる点に注意

#### 2.1. Windowsの場合
**2.1.1. Putty**

リモートログインクライアントソフトのPutty をインストールして、
このクライアントソフトから共有サーバーへ接続する

ダウンロードサイト: ：https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html

インストール方法参照: https://www.ep.sci.hokudai.ac.jp/~epnetfan/tebiki/server-login/putty.html


```
login as: [共有サーバーのアカウント名]
xxxxx@server.address's password: 
# [パスワードの入力]

[xxxxx@server.address ~]$

# 共有フォルダの確認
$ ls /mnt/bioInfo/bioInfo2021_share
```


#### 2.2. Macの場合

Macにデフォルトでインストールされているコマンドラインツール「ターミナル」で
サーバーに接続する

**ターミナル**

`アプリケーション > ユーティリティ > ターミナル`

Dockに追加しておいたほうが便利



```
$ ssh [共有サーバーのアカウント名]@server.address
xxxxx@server.address's password: 
# [パスワードの入力]
Last login: Thu Dec 17 15:17:20 2020 from XX

[xxxxx@server.address ~]$

# 共有フォルダの確認
[xxxxx@server.address ~]$ ls /mnt/bioInfo/bioInfo2021_share
```

