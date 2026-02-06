# Hail VCF AF Pipeline (Dockerized)

This repository contains a Dockerized Hail-based pipeline (`vcf-af-pipeline.py`) plus helper modules to perform VCF preprocessing/QC and downstream allele-frequency annotation. The container also includes **GrafAnc** for ancestry inference when ancestry SNPs are available.

## What you get

- Reproducible runtime with:
  - Python 3.10
  - Java 17 (required by Hail/Spark)
  - Hail + Python deps
  - GrafAnc built inside the container and linked against HTSlib
- A `docker-compose.yml` that mounts input/output/work directories for easy execution.

## Repository layout

Typical layout:

```
.
├── Dockerfile
├── docker-compose.yml
├── config.yaml
├── vcf-af-pipeline.py
├── M1_preprocessing.py
├── M2_unrelated_samples.py
├── M3_ancestry.py
├── M4_af_annotation.py
├── utils.py
├── GrafAnc_SNPs/
├── AncSnpPopAFs.txt
├── input/      # put inputs here (mounted read-only)
├── output/     # results written here
└── work/       # intermediate Hail outputs (MatrixTables, temp files)
```

## Prerequisites

- Docker Engine + Docker Compose plugin
  - Verify:
    ```bash
    docker --version
    docker compose version
    ```

## Quickstart

1) **Create folders**
```bash
mkdir -p input output work
```

2) **Place your inputs**
Put your input files in `./input` (examples: VCF, metadata CSV, etc.).  
Your `config.yaml` should reference container paths under `/data/...` so results persist (see below).

3) **CHARR Contamination Filtering Reference**

If you plan to run CHARR contamination filtering (recommended), you’ll need to download the reference database for the genome version of your data. Bear in mind that, if all your data belongs to the same version, you'll only need to download that reference. 

* **If your data was aligned to GRCh37**

To check the size and contents of the folder before downloading:

```
gsutil ls -l gs://gcp-public-data--gnomad/release/2.1.1/ht/genomes/gnomad.genomes.r2.1.1.sites.ht/
```

To download the reference data: 

```
gsutil cp -r gs://gcp-public-data--gnomad/release/2.1.1/ht/genomes/gnomad.genomes.r2.1.1.sites.ht/ .
```

* **If your data was aligned to GRCh38**

To check the size and contents of the folder before downloading:

```
gsutil ls -l gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/
```

To download the reference data: 

```
gsutil cp -r gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/ .
```

Once you have the reference table, move the files to the mounted directory: 

```
mv gnomad.genomes.vX.X.sites.ht input/
```

4) **Build the image**
```bash
docker compose build
```

5) **Run the pipeline**
```bash
docker compose up
```

Detached:
```bash
docker compose up -d
docker compose logs -f
```

Stop / cleanup:
```bash
docker compose down
```

## Volumes and paths 

The compose file mounts:

- `./input`  → `/data/input` (read-only)
- `./output` → `/data/output`
- `./work`   → `/data/work`

To ensure outputs appear on your host machine, configure `config.yaml` to use `/data/...` paths, e.g.:

```yaml
mt_from_vcf: "/data/work/synthetic.mt"
mt_afterQC: "/data/work/synthetic_afterQC.mt"
spark_local_dir: "/data/work/tmp"
tmp_dir: "/data/work/tmp"
local_tmpdir: "/data/work/tmp"
final_vcf_AF: "/data/output/final.vcf.bgz"
```

> `config.yaml` is mounted as a volume. Changing it doesn't require rebuild.  

## Running an interactive shell (debugging)

```bash
docker compose run --rm af_hail_pipeline bash
```

Inside the container you can test:

```bash
java -version
python -c "import hail as hl; hl.init(); print(hl.__version__)"
which grafanc
ldd "$(which grafanc)" | grep hts || true
```

## Troubleshooting

### `grafanc: error while loading shared libraries: libhts.so.X`
This indicates the dynamic linker cannot find HTSlib. The provided Dockerfile installs HTSlib and runs `ldconfig` to register `/usr/local/lib`. If you see this error, rebuild without cache:

```bash
docker compose build --no-cache
```

### "0 ancestry SNPs" during ancestry step
This can be expected depending on your test dataset / reference SNP set. In that case the ancestry annotation may be skipped.

### Where are my outputs?
If outputs are written to `/app/...` inside the container, they will **not** persist. Ensure paths in `config.yaml` point to `/data/output` or `/data/work`.

## License

Add your preferred license here (or keep repository default).


