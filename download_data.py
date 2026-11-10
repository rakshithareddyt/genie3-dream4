"""
download_data.py
Downloads the DREAM4 In Silico Size-100 Multifactorial data and gold standards.

Source: Bioconductor 'DREAM4' experiment-data package (public git repo).
It bundles the original DREAM4 challenge files as plain TSV.
Raw data is NOT committed to GitHub; run this script to fetch it.

Output: data/insilico_size100_{1..5}/{multifactorial,wildtype,knockouts,knockdowns,goldStandard}.tsv
"""
import subprocess, tarfile, shutil
from pathlib import Path

REPO = "https://git.bioconductor.org/packages/DREAM4"
TMP = Path("_dream4_tmp")
OUT = Path("data")

def main():
    if OUT.exists() and any(OUT.glob("insilico_size100_*/multifactorial.tsv")):
        print("Data already present in ./data - nothing to do.")
        return
    # 1) shallow-clone the package (about 3 MB)
    if TMP.exists():
        shutil.rmtree(TMP)
    subprocess.run(["git", "clone", "--depth", "1", REPO, str(TMP)], check=True)
    # 2) unpack the tarball holding the original challenge files
    tar_path = TMP / "inst/extdata/lightlyProcessedDownloadedData.tar.gz"
    with tarfile.open(tar_path) as tar:
        tar.extractall(TMP / "extracted")
    src = TMP / "extracted/lightlyProcessedDownloadedData"
    # 3) keep only what we need, skipping macOS '._' junk files
    OUT.mkdir(exist_ok=True)
    for i in range(1, 6):
        d = OUT / f"insilico_size100_{i}"
        d.mkdir(exist_ok=True)
        for fname in ["multifactorial.tsv", "goldStandard.tsv", "wildtype.tsv", "knockouts.tsv", "knockdowns.tsv"]:
            shutil.copy(src / f"insilico_size100_{i}" / fname, d / fname)
    shutil.rmtree(TMP)
    print("Done. Files in ./data")

if __name__ == "__main__":
    main()
