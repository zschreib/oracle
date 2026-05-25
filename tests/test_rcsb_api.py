"""
Standalone test for RCSB PDB sequence similarity search API.
Tests both experimental PDB and AlphaFold CSM searches.

Usage:
    python tests/test_rcsb_api.py
"""

import sys
import requests

# Short dark matter-like sequence - unlikely to hit experimental PDB
# so should trigger CSM fallback
DARK_SEQUENCE = (
    "MSKKQTITGIVLGLIAGTTLSADNKTATNLNAKNDLINAAKSDLTNAKSDLTNAKNELT"
    "NAKADLTNAKSDLTDSKTDLTDAKTDLADAKSDLTDAKNDLTNAKADLT"
)

# T4 lysozyme - should hit experimental PDB strongly
STRONG_SEQUENCE = (
    "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITK"
    "DEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRM"
    "LQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL"
)

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def search(sequence: str, content_type: str = "experimental",
           identity_cutoff: float = 0.20, max_results: int = 3) -> dict:
    """Search RCSB with specified content type."""
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": sequence
            }
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "results_content_type": [content_type],
            "paginate": {"start": 0, "rows": max_results}
        },
        "return_type": "polymer_entity"
    }

    response = requests.post(
        RCSB_SEARCH_URL,
        json=query,
        headers={"Content-Type": "application/json"},
        timeout=30
    )
    print(f"    Status: {response.status_code}")

    if response.status_code == 204 or not response.text.strip():
        return {"total_count": 0, "result_set": []}

    response.raise_for_status()
    return response.json()


def main():
    print("=" * 60)
    print("RCSB PDB search API test - experimental + CSM")
    print("=" * 60)

    # Test 1: Strong hit against experimental PDB
    print("\nTest 1: Strong sequence against experimental PDB")
    print(f"  Sequence: T4 lysozyme ({len(STRONG_SEQUENCE)} aa)")
    try:
        results = search(STRONG_SEQUENCE, content_type="experimental")
        total = results.get("total_count", 0)
        hits = results.get("result_set", [])
        print(f"  Total matches: {total}")
        if hits:
            print(f"  Top hit: {hits[0].get('identifier')} score={hits[0].get('score')}")
            print("  PASS - experimental search works")
        else:
            print("  WARN - no hits returned")
    except Exception as e:
        print(f"  FAILED: {e}")

    # Test 2: Dark matter against experimental PDB (expect no/few hits)
    print("\nTest 2: Dark matter sequence against experimental PDB")
    print(f"  Sequence: BATS dark matter ({len(DARK_SEQUENCE)} aa)")
    try:
        results = search(DARK_SEQUENCE, content_type="experimental")
        total = results.get("total_count", 0)
        hits = results.get("result_set", [])
        print(f"  Total matches: {total}")
        if hits:
            print(f"  Top hit: {hits[0].get('identifier')} score={hits[0].get('score')}")
        else:
            print("  No experimental hits (expected for dark matter)")
        print("  PASS - experimental search completed")
    except Exception as e:
        print(f"  FAILED: {e}")

    # Test 3: Dark matter against CSM (AlphaFold) - the key test
    print("\nTest 3: Dark matter sequence against AlphaFold CSMs")
    print(f"  Sequence: BATS dark matter ({len(DARK_SEQUENCE)} aa)")
    print("  content_type: computational")
    try:
        results = search(DARK_SEQUENCE, content_type="computational",
                        identity_cutoff=0.20)
        total = results.get("total_count", 0)
        hits = results.get("result_set", [])
        print(f"  Total matches: {total}")
        if hits:
            print(f"  Top hit: {hits[0].get('identifier')} score={hits[0].get('score')}")
            print("  PASS - CSM search works and found hits")
        else:
            print("  No CSM hits either - genuine dark matter")
            print("  PASS - CSM search completed without error")
    except Exception as e:
        print(f"  FAILED: {e}")
        import traceback
        traceback.print_exc()

    # Test 4: Timeout behavior - increase timeout
    print("\nTest 4: CSM search with longer timeout (60s)")
    try:
        query = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 1,
                    "identity_cutoff": 0.20,
                    "sequence_type": "protein",
                    "value": DARK_SEQUENCE
                }
            },
            "request_options": {
                "scoring_strategy": "sequence",
                "results_content_type": ["computational"],
                "paginate": {"start": 0, "rows": 3}
            },
            "return_type": "polymer_entity"
        }
        response = requests.post(
            RCSB_SEARCH_URL,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=60  # longer timeout
        )
        print(f"    Status: {response.status_code}")
        print(f"    Response time: completed within 60s")
        if response.status_code in (200, 204):
            print("    PASS - CSM search works with 60s timeout")
        else:
            print(f"    Response: {response.text[:200]}")
    except Exception as e:
        print(f"    FAILED even with 60s timeout: {e}")

    print()
    print("=" * 60)


if __name__ == "__main__":
    main()


def search_sequence(sequence: str, identity_cutoff: float = 0.3,
                    max_results: int = 5) -> dict:
    """
    Search RCSB PDB for structures with similar sequences using MMseqs2.

    Args:
        sequence:         Amino acid sequence string
        identity_cutoff:  Minimum sequence identity (0.0-1.0). Default 0.3
        max_results:      Maximum hits to return

    Returns:
        Full JSON response from RCSB Search API
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": sequence
            }
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "paginate": {
                "start": 0,
                "rows": max_results
            }
        },
        "return_type": "polymer_entity"
    }

    response = requests.post(
        RCSB_SEARCH_URL,
        json=query,
        headers={"Content-Type": "application/json"},
        timeout=30
    )
    print(f"  Status code: {response.status_code}")
    if response.status_code != 200:
        print(f"  Response: {response.text[:300]}")
    response.raise_for_status()
    return response.json()


def fetch_entity_info(pdb_id: str, entity_id: str) -> dict:
    """
    Fetch functional annotation for a PDB polymer entity.

    Returns name, source organism, and UniProt annotations.
    """
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    response = requests.get(url, timeout=15)
    if response.status_code == 200:
        return response.json()
    return {}



if __name__ == "__main__":
    main()