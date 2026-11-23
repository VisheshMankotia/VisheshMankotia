# CDK2–Tivozanib Molecular Docking 🔬

A reproducible learning workflow for preparing a protein and ligand, performing molecular docking, and inspecting predicted protein–ligand interactions.

> **Status:** Computational study / work in progress. This documentation records the workflow and an example docking output previously observed during practical work. It is not evidence of experimental binding or biological efficacy.

## Study overview
- **Target:** Cyclin-dependent kinase 2 (CDK2), PDB entry **1AQ1**, chain A (confirm chain and biological assembly before use).
- **Ligand:** Tivozanib, PubChem CID **9911830**.
- **Docking software:** AutoDock 4 / AutoDockTools.
- **Supporting tools:** Open Babel, UCSF Chimera, PyRx (where applicable), GROMACS (for a separate molecular-dynamics learning workflow).

## Workflow
1. Obtain the target structure from the Protein Data Bank and verify the biological context, chain, missing residues and bound molecules.
2. Prepare the protein: remove only molecules not required for the chosen protocol, add hydrogens, assign charges, and save in PDBQT format.
3. Obtain the ligand structure from PubChem, inspect stereochemistry and protonation, minimize it, and prepare ligand PDBQT.
4. Define the docking grid around a justified binding site. Record centre coordinates, dimensions and spacing.
5. Run AutoDock 4 with documented genetic-algorithm settings.
6. Cluster poses, compare binding energies and inspect key contacts.
7. Record software versions, input structure identifiers and all parameters for reproducibility.

## Example practical settings
| Parameter | Recorded setting |
|---|---:|
| Grid points | 60 × 60 × 60 |
| Grid spacing | 0.375 Å |
| Genetic algorithm runs | 10 |

The grid centre must be determined from the selected binding site and recorded in the run configuration; do not guess coordinates.

## Example output noted during practical work
One reported pose had an estimated binding free energy of **−6.63 kcal/mol** and an estimated inhibition constant (Ki) of **13.80 µM**. Re-run the workflow from the documented inputs before treating this as a reproducible result.

## Results to document
- Protein and ligand source/version.
- Protein preparation decisions.
- Ligand protonation, stereochemistry and charge method.
- Grid centre and dimensions.
- Exhaustiveness / GA settings and random seeds where supported.
- Best pose, cluster population and energy distribution.
- Hydrogen bonds, hydrophobic contacts, salt bridges and interacting residues.
- Limitations and comparison with experimental evidence.

## Suggested folder structure
```
CDK2-Tivozanib-Molecular-Docking/
├── README.md
├── config/
│   └── docking_parameters_template.md
├── results/
│   └── results_template.csv
└── references/
    └── sources.md
```

## Interpretation and limitations
Docking scores are approximate scoring-function outputs, not measured binding affinities. They depend on structure preparation, protonation, receptor flexibility, water treatment, grid placement and scoring assumptions. A docking pose alone cannot establish inhibition, selectivity, anticancer activity or clinical relevance. Molecular dynamics simulations and experimental assays require separate validation.

## Author
**Vishesh Mankotia** · [GitHub](https://github.com/VisheshMankotia) · [LinkedIn](https://www.linkedin.com/in/vishesh-mankotia-b8569b292/)
