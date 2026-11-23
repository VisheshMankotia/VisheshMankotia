"""Calculate basic RDKit descriptors for a curated SMILES table.

Input CSV columns: compound, smiles
Output: compound_descriptors.csv

This is a starter utility; inspect and curate structures before interpreting results.
"""
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

INPUT = Path("data/compounds.csv")
OUTPUT = Path("results/compound_descriptors.csv")


def calculate(smiles: str):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return {
        "canonical_smiles": Chem.MolToSmiles(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "h_bond_donors": rdMolDescriptors.CalcNumHBD(mol),
        "h_bond_acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }


def main():
    if not INPUT.exists():
        raise FileNotFoundError(
            f"Create {INPUT} with columns 'compound' and 'smiles'. "
            "Do not use unverified or fabricated structures."
        )

    data = pd.read_csv(INPUT)
    required = {"compound", "smiles"}
    if not required.issubset(data.columns):
        raise ValueError(f"Input must contain columns: {sorted(required)}")

    rows = []
    for _, row in data.iterrows():
        descriptors = calculate(row["smiles"])
        rows.append({
            "compound": row["compound"],
            "smiles": row["smiles"],
            **(descriptors or {"error": "Invalid SMILES"}),
        })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUTPUT, index=False)
    print(f"Saved {len(rows)} records to {OUTPUT}")


if __name__ == "__main__":
    main()
