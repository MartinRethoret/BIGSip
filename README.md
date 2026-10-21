# BIGSip

## 📌 Description

**Sip just the BIGSdb data you need.**

BIGSip downloads data from BIGSdb databases hosted by:

- PubMLST
- BIGSdb-Pasteur

You give an organism name and what you want; BIGSip builds the API requests and runs them through [BIGSdb_downloader](https://github.com/kjolley/BIGSdb_downloader), which handles the authentication required for data produced after 1 January 2025.

It retrieves:

- Genomes of isolates (contigs FASTA)
- Loci (allele FASTA files), for a whole scheme or for chosen loci
- Profiles (TSV files), for schemes that have a primary key (ST, cgST...)

## ⚠️ Data usage policies

The data retrieved using this script are subject to the usage policies of each respective platform. Users are responsible for complying with these terms of use.

👉 [PubMLST Policy](https://pubmlst.org/terms-conditions)

👉 [BIGSdb-Pasteur Policy](https://bigsdb.pasteur.fr/policy/)

## ⚙️ Installation

### Requirements

- Python ≥ 3.9
- [BIGSdb_downloader](https://github.com/kjolley/BIGSdb_downloader)

BIGSip itself only uses the Python standard library: installing BIGSdb_downloader is enough.

```
python -m venv ~/venvs/bigsdb
~/venvs/bigsdb/bin/pip install bigsdb-downloader
```

Then copy `bigsip.py` anywhere and run it with that environment's Python:

```
~/venvs/bigsdb/bin/python bigsip.py --version
```

For other ways to install BIGSdb_downloader (pipx, cron...), see its [README](https://github.com/kjolley/BIGSdb_downloader#installation).

### API key

Authenticated access needs an account and a personal API key:

- BIGSdb-Pasteur: [request an API key](https://bigsdb.pasteur.fr/requesting-api-key/)
- PubMLST: create a personal key at <https://pubmlst.org/bigsdb>

See the [BIGSdb_downloader README](https://github.com/kjolley/BIGSdb_downloader#accessing-pubmlst-and-bigsdb-pasteur-apis-using-authentication) for details. Without a key, `--public` gives access to pre-2025 data only.

### First-time setup

```
python bigsip.py setup                          # path to bigsdb-downloader + personal token directory
python bigsip.py auth --organism diphtheria     # link your API key (interactive, once per site)
python bigsip.py status                         # check the configuration
```

`auth` runs the BIGSdb_downloader credential setup for you: enter your client id and secret, open the URL shown, paste the verification code.

## 🚀 Usage

Run without arguments for the interactive wizard:

```
python bigsip.py
```

### 1. List organisms and schemes

```
python bigsip.py organisms
python bigsip.py schemes --organism klebsiella
python bigsip.py schemes --organism klebsiella --scheme 4    # loci of scheme 4
```

### 2. Download genomes

```
python bigsip.py genomes --organism klebsiella --ids 1,2,10-20
python bigsip.py genomes --organism klebsiella --ids ids.txt
```

### 3. Download loci

```
python bigsip.py loci --organism klebsiella --scheme 4
python bigsip.py loci --organism klebsiella --loci gapA,infB
python bigsip.py loci --organism klebsiella --loci loci.txt
```

### 4. Download profiles

```
python bigsip.py profiles --organism klebsiella --scheme 4
```

Text files contain one item per line. Isolate and scheme IDs are numeric; organisms and loci are given by name.

Useful options: `--outdir`, `--dry-run` (show the commands without running them), `--overwrite`, `--jobs N`, `--site PubMLST`, `--public`, `-v` / `-q` (more / less messages). See `python bigsip.py <command> -h`.

## 📂 Output structure

```
outdir/
└── <organism>/
    ├── genomes/
    │   ├── 1.fasta
    │   └── ...
    ├── alleles/
    │   └── <locus>.fasta
    ├── scheme_<id>_<name>/
    │   ├── alleles/
    │   │   └── <locus>.fasta
    │   ├── loci.txt
    │   └── profiles.tsv
    ├── download_log.tsv
    └── failed_<type>_<date>.txt
```

Files already present are skipped, so a command can be re-run to resume. Each file is checked (FASTA / TSV) before being kept, so an error page is never saved as data. Failed items are listed in `failed_*.txt`, which can be passed back to `--ids` or `--loci`.

`download_log.tsv` records every run (date, BIGSip version, site, database, request, counts) to keep track of where the data come from.

## 👥 Shared installation

One copy of `bigsip.py` can be shared by several users, each keeping their own tokens:

```
python bigsip.py setup --shared --downloader /opt/bigsdb/venv/bin/bigsdb-downloader --token-dir "~/.bigsip/tokens"
```

Each user then only runs `python bigsip.py auth --organism <organism>`.

## ⚠️ Limitations

- Profiles are only available for schemes with a primary key; other schemes can be downloaded with `loci --scheme`
- Requests are subject to remote server limitations
- Large downloads may take significant time
