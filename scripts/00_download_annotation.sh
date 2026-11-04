#!/usr/bin/env bash
# BTx623 v5.1 gene annotation (JGI, Fort Lauderdale data-use terms) from the SorghumBase FTP, pinned to release-11.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/annotation
F=Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.gff3.gz
URL=https://ftp.sorghumbase.org/release-11/gff3/sorghum_bicolorv5/$F
MD5=9140e1d3787919306025a7f1d7679976
if [ ! -s data/annotation/$F ]; then
  curl -L --fail -o data/annotation/$F.part "$URL" && mv data/annotation/$F.part data/annotation/$F
fi
echo "$MD5  data/annotation/$F" | md5sum -c -
