"""
Standalone test for the ESM Atlas fold API.

Tests structure prediction from sequence, then prints the first 50 lines
of PDB output so we know what we're working with before building the tool.

Usage:
    python tests/test_esmfold_api.py
"""

import sys
import requests

# Same sequence from the ESM docs example
TEST_SEQUENCE = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINS"
    "RWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQ"
    "AWIRGCRL"
)

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

def main():
    print("=" * 60)
    print("ESMFold API test")
    print("=" * 60)
    print(f"Sequence length: {len(TEST_SEQUENCE)} aa")
    print(f"Endpoint: {ESMFOLD_URL}")
    print()

    print("Submitting sequence (this may take 10-30s)...")
    try:
        response = requests.post(
            ESMFOLD_URL,
            data=TEST_SEQUENCE,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=60
        )
        print(f"Status code: {response.status_code}")
        print(f"Content-Type: {response.headers.get('Content-Type', 'unknown')}")
        print(f"Response length: {len(response.text)} chars")
        print()

        if response.status_code == 200:
            lines = response.text.splitlines()
            print(f"First 20 lines of response:")
            print("-" * 40)
            for line in lines[:20]:
                print(line)
            print("-" * 40)
            print(f"Total lines: {len(lines)}")
            print()

            # Check if it looks like a PDB file
            is_pdb = any(
                line.startswith(("ATOM", "HETATM", "HEADER", "REMARK"))
                for line in lines[:10]
            )
            print(f"Looks like PDB format: {is_pdb}")

            if is_pdb:
                atom_lines = [l for l in lines if l.startswith("ATOM")]
                print(f"ATOM records: {len(atom_lines)}")
                print()
                print("PASS - ESMFold API is working and returns PDB format")
            else:
                print(f"Response preview: {response.text[:300]}")
                print("WARN - Response received but may not be PDB format")
        else:
            print(f"FAILED - Response: {response.text[:300]}")

    except requests.exceptions.Timeout:
        print("FAILED - Request timed out after 60s")
        sys.exit(1)
    except Exception as e:
        print(f"FAILED - {e}")
        sys.exit(1)

    print("=" * 60)


if __name__ == "__main__":
    main()