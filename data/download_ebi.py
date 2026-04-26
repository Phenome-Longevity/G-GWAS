#!/usr/bin/env python3
"""
Step 0a: Download harmonised GWAS summary statistics from EBI GWAS Catalog.

Downloads all available harmonised TSV files via FTP.
Each file: ~10-500MB compressed, contains per-variant summary statistics.

Usage:
    python data/download_ebi.py --output_dir /path/to/raw_gwas/ --max_files 50000

Output: directory of GCST*.h.tsv.gz files.
"""

import os, time, argparse, json
import urllib.request
import urllib.error


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


FTP_BASE = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"


def get_study_list():
    """Get list of available GWAS studies from EBI API."""
    log("Fetching study list from EBI GWAS Catalog API...")
    url = "https://www.ebi.ac.uk/gwas/rest/api/studies?size=1000"
    studies = []
    page = 0
    while url:
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read())
            for s in data.get('_embedded', {}).get('studies', []):
                acc = s.get('accessionId', '')
                if acc.startswith('GCST'):
                    studies.append(acc)
            links = data.get('_links', {})
            url = links.get('next', {}).get('href', None)
            page += 1
            if page % 10 == 0:
                log(f"  Page {page}, {len(studies)} studies so far")
        except Exception as e:
            log(f"  API error at page {page}: {e}")
            break
    log(f"Found {len(studies)} studies")
    return studies


def download_harmonised(gcst, output_dir):
    """Download harmonised summary statistics for a single study."""
    out_path = os.path.join(output_dir, f"{gcst}.h.tsv.gz")
    if os.path.exists(out_path):
        return 'skip'

    # EBI FTP structure: GCSTXXXXXX/ contains harmonised/ subdirectory
    prefix = gcst[:8]  # e.g., GCST0005
    ftp_dir = f"{FTP_BASE}{gcst}/"

    try:
        # Try direct harmonised file pattern
        url = f"{ftp_dir}harmonised/{gcst}.h.tsv.gz"
        urllib.request.urlretrieve(url, out_path)
        return 'ok'
    except urllib.error.URLError:
        pass

    try:
        # Alternative pattern
        url = f"{ftp_dir}harmonised/"
        with urllib.request.urlopen(url, timeout=30) as resp:
            listing = resp.read().decode()
        # Parse for .h.tsv.gz files
        for line in listing.split('"'):
            if line.endswith('.h.tsv.gz'):
                file_url = f"{ftp_dir}harmonised/{line}"
                urllib.request.urlretrieve(file_url, out_path)
                return 'ok'
    except:
        pass

    return 'fail'


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output_dir', required=True, help='Directory for raw GWAS files')
    p.add_argument('--max_files', type=int, default=50000)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log("=" * 70)
    log("DOWNLOAD EBI GWAS CATALOG")
    log("=" * 70)

    studies = get_study_list()
    stats = {'ok': 0, 'skip': 0, 'fail': 0}

    for i, gcst in enumerate(studies[:args.max_files]):
        result = download_harmonised(gcst, args.output_dir)
        stats[result] += 1
        if (i + 1) % 100 == 0:
            log(f"  {i+1}/{len(studies)} | ok={stats['ok']} skip={stats['skip']} fail={stats['fail']}")

    log(f"\nDone. Stats: {stats}")
    log(f"Output: {args.output_dir}")
