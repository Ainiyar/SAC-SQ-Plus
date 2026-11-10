#!/usr/bin/env python3
"""Download the three UCI source archives used by the fitted-cost experiments."""
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "occupancy": ("https://archive.ics.uci.edu/static/public/357/occupancy%2Bdetection.zip", "uci_occupancy.zip"),
    "eeg": ("https://archive.ics.uci.edu/static/public/264/eeg%2Beye%2Bstate.zip", "uci_eeg.zip"),
    "secom": ("https://archive.ics.uci.edu/static/public/179/secom.zip", "uci_secom.zip"),
}


def main():
    for name, (url, filename) in SOURCES.items():
        target_dir = ROOT / name
        target_dir.mkdir(parents=True, exist_ok=True)
        target = target_dir / filename
        if target.exists() and target.stat().st_size:
            print(f"EXISTS {target} ({target.stat().st_size} bytes)")
            continue
        print(f"DOWNLOAD {name}: {url}", flush=True)
        with urlopen(url, timeout=120) as response:
            payload = response.read()
        if not payload:
            raise RuntimeError(f"Empty UCI response for {name}")
        target.write_bytes(payload)
        print(f"SAVED {target} ({len(payload)} bytes)", flush=True)


if __name__ == "__main__":
    main()
