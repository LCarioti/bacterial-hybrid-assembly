# Bacterial Hybrid Assembly Pipeline 

A robust Bash pipeline for automated hybrid de novo assembly of bacterial genomes, combining long, high-error reads (Oxford Nanopore) with short, high-accuracy reads (Illumina). 

The pipeline performs automated quality control, read filtering, long-read assembly, short-read polishing, evaluation, and structural re-mapping.

##  Workflow Overview
1. **Short-read QC**: Automated adapter trimming and quality filtering using `fastp`.
2. **Long-read Filtering**: Size and quality selection using `filtlong` (supports Standard and Aggressive modes).
3. **De Novo Assembly**: Long-read assembly using `flye` (`--nano-hq`).
4. **Polishing**: Hybrid error-correction using Illumina reads via `POLCA` (MaSuRCA).
5. **Evaluation**: Assembly metrics generation against a reference using `QUAST`.
6. **Remapping & Coverage**: Reads are mapped back to the polished assembly using `minimap2` (ONT) and `bwa mem` (Illumina) for coverage verification.

##  Tools & Installation

This pipeline orchestrates several widely used bioinformatic tools. Since all tools are available on GitHub, you can clone and compile them manually. Ensure that all binaries or wrappers are accessible in your system's `$PATH`.

### Dependency List & Sources

| Tool | Purpose | GitHub Repository / Source |
| :--- | :--- | :--- |
| **fastp** | Short-read QC and adapter trimming | [OpenGene/fastp](https://github.com/opengene/fastp) |
| **filtlong** | Long-read quality and length filtering | [rrwick/Filtlong](https://github.com/rrwick/Filtlong) |
| **Flye** | De novo long-read assembly | [mikolayenko/Flye](https://github.com/rrwick/Filtlong) |
| **MaSuRCA (POLCA)** | Polishing using short reads | [alekseyzimin/masurca](https://github.com/alekseyzimin/masurca) |
| **QUAST** | Quality assessment of the assembly | [ablab/quast](https://github.com/ablab/quast) |
| **minimap2** | Long-read mapping against assembly | [lh3/minimap2](https://github.com/lh3/minimap2) |
| **BWA** | Short-read indexing and mapping | [lh3/bwa](https://github.com/lh3/BWA) |
| **samtools** | BAM file sorting, indexing, and statistics | [samtools/samtools](https://github.com/samtools/samtools) |

---

### Manual Installation from GitHub (Example Guide)

If you wish to replicate the native installation, follow these steps for each tool.

#### 1. Core Tools (Compilation required)
For tools written in C/C++ like `fastp`, `minimap2`, `bwa`, and `samtools`, clone the repo and run `make`:
```bash
# Example for Minimap2
git clone https://github.com.git
cd minimap2 && make
# Repeat similar compilation steps for bwa, fastp, and samtools as per their GitHub readmes.
```

#### 2. Python & C++ Hybrid Tools
Tools like `Flye` and `Filtlong` require compilation and a proper Python environment:
```bash
# Example for Filtlong
git clone https://github.com.git
cd Filtlong && make bin/filtlong

# Example for Flye
git clone https://github.com.git
cd Flye && python setup.py install
```

#### 3. Polishing & Evaluation Tools
* **POLCA** is bundled inside **MaSuRCA**. Clone [alekseyzimin/masurca](https://github.com/alekseyzimin/masurca), follow their installer script, and ensure `polca.sh` (or `polca`) is extracted.
* **QUAST** requires Python dependencies. Clone [ablab/quast](https://github.com/ablab/quast) and install it via `python setup.py install`.

---

###  Making Tools Globally Accessible

To make the pipeline universal, the script searches for tools dynamically in your environment. You **must** add the paths of your compiled binaries to your user configuration (e.g., `~/.bashrc` or `~/.zshrc`):

```bash
# Open your bashrc
nano ~/.bashrc

# Add your GitHub installation folders to the PATH (Modify with your actual paths)
export PATH="/path/to/github/fastp:\$PATH"
export PATH="/path/to/github/minimap2:\$PATH"
export PATH="/path/to/github/bwa:\$PATH"
export PATH="/path/to/github/samtools:\$PATH"
export PATH="/path/to/github/Filtlong/bin:\$PATH"
export PATH="/path/to/github/Flye/bin:\$PATH"
export PATH="/path/to/github/masurca/bin:\$PATH"
export PATH="/path/to/github/quast:\$PATH"

# Save, exit, and reload the terminal
source ~/.bashrc
```

The pipeline will automatically run a pre-flight check to verify if all 8 tools are correctly exposed and executable before processing your sequencing data.

## 📖 Usage

```bash
./bacterial_pipeline.sh <R1.fq.gz> <R2.fq.gz> <nanopore.fq.gz> <reference.fasta> <threads> [mode]
```

### Arguments:
* `<R1>` / `<R2>`: Forward and reverse Illumina paired-end reads (FastQ compressed).
* `<nanopore.fastq.gz>`: Raw Oxford Nanopore reads.
* `<reference.fasta>`: Close reference genome for QUAST evaluation.
* `<threads>`: Number of CPU cores to allocate.
* `[mode]`: (*Optional*) Assembly mode:
  * `standard`: Keeps ONT reads > 2000 bp. Ideal for preserving small plasmids. (Default)
  * `aggressive`: Keeps ONT reads > 6000 bp (top 90%). Best for resolving complex repeats, rRNA operons, and closing chromosomes.

### Example:
```bash
./bacterial_pipeline.sh reads_R1.fastq.gz reads_R2.fastq.gz ont_raw.fastq.gz ref.fasta 16 aggressive
```

##  Output Structure
The pipeline organizes outputs into clean, structured directories:
* `raw/` & `reference/`: Symlinks to your input datasets.
* `qc/`: Trimmed Illumina reads and `fastp` HTML/JSON reports.
* `nanopore/`: Size-filtered ONT reads.
* `assembly/`: Raw Flye output, polished genome (`*.PolcaCorrected.fa`), and `QUAST` reports.
* `mapping/`: Sorted and indexed BAM files for both ONT and Illumina reads, including `flagstat` coverage statistics.
* `logs/`: Complete stdout/stderr console logs for reproducibility.

##  Next Steps & Downstream Analysis

Once the pipeline completes successfully, you will obtain a polished, high-quality draft genome (`assembly/assembly.fasta.PolcaCorrected.fa`). To fully characterize your bacterial isolate, we recommend the following downstream analyses:

### 1. Functional Annotation with Bakta
To predict genes, tRNA, rRNA, and functionally annotate your genome, **Bakta** is highly recommended due to its updated databases and antimicrobial resistance (AMR) gene detection. You can use it via the web interface or run it locally.
* **Web Interface:** Upload your `*.PolcaCorrected.fa` directly to the official [Bakta Web Server](https://bakta.computational.bio/).
* **Command-line Software:** Clone the repository from [OSF-Biolab/bakta](https://github.com/oschwengers/bakta) and run it locally:
  ```bash
  bakta --db /path/to/bakta_db --threads <threads> --output assembly/annotation/ assembly/assembly.fasta.PolcaCorrected.fa
  ```

### 2. Multi-Locus Sequence Typing (MLST) via Institut Pasteur
To identify the Sequence Type (ST) and clonal complex of your isolate for epidemiological tracking, you can use the official **Institut Pasteur MLST databases**.
* **Web Interface:** You can directly upload your `*.PolcaCorrected.fa` file to the [Institut Pasteur MLST Web Portal](https://pasteur.fr).
* **Command Line Alternative:** Alternatively, you can use the command-line tool `mlst` (by Torsten Seemann) which includes Pasteur schemes:
  ```bash
  mlst assembly/assembly.fasta.PolcaCorrected.fa
  ```
  
##  Citations
If you use this pipeline in your research, please remember to cite the individual tools utilized in the workflow (links to publications can be found in their respective GitHub repositories).

##  License
This project is licensed under the **GNU General Public License v3.0**. You are free to copy, modify, and distribute this software, provided that any derivative works are also licensed under the GPLv3. See the [LICENSE](LICENSE) file for details.

