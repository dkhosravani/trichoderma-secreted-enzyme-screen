"""Download an entire UniProtKB proteome with pagination and provenance."""
import json
from datetime import datetime, timezone
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parent.parent
QUERY = 'proteome:UP000008984'
URL = 'https://rest.uniprot.org/uniprotkb/search'


def main():
    session = requests.Session()
    retry = Retry(total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504],
                  respect_retry_after_header=True)
    session.mount('https://', HTTPAdapter(max_retries=retry))

    records, accessions = [], set()
    expected_total, release, release_date = None, None, None
    url = URL
    params = {'query': QUERY, 'format': 'json', 'size': 500}
    page = 0

    while url:
        response = session.get(url, params=params, timeout=90)
        response.raise_for_status()
        if expected_total is None:
            expected_total = int(response.headers['x-total-results'])
            release = response.headers.get('x-uniprot-release')
            release_date = response.headers.get('x-uniprot-release-date')

        batch = response.json().get('results', [])
        if not batch:
            raise RuntimeError(f'Empty UniProt page {page + 1}')
        for record in batch:
            accession = record['primaryAccession']
            if accession in accessions:
                raise RuntimeError(f'Duplicate accession: {accession}')
            accessions.add(accession)
            records.append(record)
        page += 1
        print(f'Page {page}: {len(records)} / {expected_total}')
        url = response.links.get('next', {}).get('url')
        params = None

    if len(records) != expected_total:
        raise RuntimeError(f'Incomplete download: {len(records)} of {expected_total}')

    (ROOT / 'data').mkdir(exist_ok=True)
    (ROOT / 'results').mkdir(exist_ok=True)
    with (ROOT / 'data/proteome_records.json').open('w', encoding='utf-8') as handle:
        json.dump(records, handle)
    provenance = {
        'source': 'UniProtKB REST API',
        'query': QUERY,
        'downloaded_at_utc': datetime.now(timezone.utc).isoformat(),
        'uniprot_release': release,
        'uniprot_release_date': release_date,
        'expected_total': expected_total,
        'downloaded_count': len(records),
        'pages': page,
    }
    with (ROOT / 'results/metadata.json').open('w', encoding='utf-8') as handle:
        json.dump(provenance, handle, indent=2)
    print(f'Download complete: {len(records)} unique records')


if __name__ == '__main__':
    main()
