#Fetches structures from the AlphaFold database
#by UniProt accession

#API docs: https://alphafold.ebi.ac.uk/api-docs


import os
import requests

API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"
CACHE_DIR = os.path.join(os.path.dirname(__file__), "cache")

os.makedirs(CACHE_DIR, exist_ok=True)


class AlphaFoldError(Exception):
    pass


def get_prediction_metadata(uniprot_id: str) -> dict:
    """
    Looks up the AlphaFold prediction for a UniProt accession.
    Returns the metadata dict, which includes pdbUrl, cifUrl, uniprotSequence,
    and confidence info -- see the AlphaFold API docs for the full shape.
    """
    #looks up prediction of 3d structure, returns meta data dictionary with the pdbUrl,
    #ciUrl, uniprotSequence and confidence info
    url = f"{API_BASE}/{uniprot_id}"
    response = requests.get(url, timeout=15)
    if response.status_code != 200:
        raise AlphaFoldError(
            f"No AlphaFold prediction found for {uniprot_id} "
            f"(status {response.status_code})"
        )
    data = response.json()
    if not data:
        raise AlphaFoldError(f"Empty response for {uniprot_id}")
    #API returns a list
    return data[0]


def get_pdb_path(uniprot_id: str) -> str:
    #returns local file path to the pdb structure for the uniprot accession and 
    #downloading if we obvi dpnt have it and then cachine cuz duh
    cache_path = os.path.join(CACHE_DIR, f"{uniprot_id}.pdb")
    if os.path.exists(cache_path):
        return cache_path

    metadata = get_prediction_metadata(uniprot_id)
    pdb_url = metadata["pdbUrl"]

    response = requests.get(pdb_url, timeout=30)
    if response.status_code != 200:
        raise AlphaFoldError(f"Could not download structure from {pdb_url}")

    with open(cache_path, "w") as f:
        f.write(response.text)

    return cache_path


def get_sequence(uniprot_id: str) -> str:
    #returns the amino acid sequence that alphafold predicted the sturcture for
    metadata = get_prediction_metadata(uniprot_id)
    return metadata["uniprotSequence"]


if __name__ == "__main__":
    import sys
    uid = sys.argv[1] if len(sys.argv) > 1 else "P69905"  # human hemoglobin alpha
    print(f"Fetching metadata for {uid}...")
    meta = get_prediction_metadata(uid)
    print("Entry ID:", meta.get("entryId"))
    print("Sequence length:", len(meta.get("uniprotSequence", "")))
    path = get_pdb_path(uid)
    print("Downloaded to:", path)