"""
Standalone test for the EBI InterProScan REST API.

Run this directly to verify the endpoint works before running the full agent.
No Anthropic API key needed, no tokens spent.

Usage:
    python tests/test_hmmer_api.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oracle.tools.hmmer import HmmerTool

# T4 lysozyme - well characterized, should return Pfam Phage_lysozyme domain
TEST_SEQUENCE = "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL"
TEST_ID = "T4_lysozyme_test"

def main():
    print("=" * 60)
    print("EBI InterProScan API test")
    print("=" * 60)
    print(f"Sequence: {TEST_ID}")
    print(f"Length:   {len(TEST_SEQUENCE)} aa")
    print()

    tool = HmmerTool(max_wait=300, poll_interval=10)

    print("Step 1: Submitting job to EBI InterProScan...")
    try:
        job_id = tool._submit(TEST_SEQUENCE)
        print(f"  Job ID: {job_id}")
    except Exception as e:
        print(f"  FAILED: {e}")
        sys.exit(1)

    print("Step 2: Polling for results (max 300s, checking every 10s)...")
    try:
        results = tool._poll(job_id)
        sequences = results.get("results", [])
        total_matches = sum(
            len(s.get("matches", [])) for s in sequences
        )
        print(f"  Status: FINISHED")
        print(f"  Sequences in result: {len(sequences)}")
        print(f"  Total matches: {total_matches}")
        if sequences:
            dbs = set()
            for s in sequences:
                for m in s.get("matches", []):
                    db = m.get("signature", {}).get(
                        "signatureLibraryRelease", {}
                    ).get("library", "unknown")
                    dbs.add(db)
            print(f"  Databases with hits: {sorted(dbs)}")
            # Print first match keys to debug db name extraction
            first_match = sequences[0].get("matches", [{}])[0]
            sig = first_match.get("signature", {})
            print(f"  First match signature keys: {list(sig.keys())}")
            print(f"  signatureLibraryRelease: {sig.get('signatureLibraryRelease', 'NOT FOUND')}")
    except Exception as e:
        print(f"  FAILED: {e}")
        sys.exit(1)

    print("Step 3: Parsing best hit...")
    try:
        evidence = tool._parse(results, TEST_SEQUENCE)
        print(f"  Annotation: {evidence.annotation}")
        print(f"  Confidence: {evidence.confidence.name}")
        print(f"  E-value:    {evidence.score}")
        print(f"  Coverage:   {evidence.coverage:.0%}" if evidence.coverage else "  Coverage: N/A")
        print(f"  Hit ID:     {evidence.hit_id}")
        print(f"  Reasoning:  {evidence.reasoning}")
    except Exception as e:
        print(f"  FAILED: {e}")
        sys.exit(1)

    print()
    print("=" * 60)
    if evidence.confidence.name in ("HIGH", "MODERATE"):
        print("PASS - InterProScan API is working correctly")
    else:
        print("WARN - API responded but confidence is low, check results above")
    print("=" * 60)


if __name__ == "__main__":
    main()