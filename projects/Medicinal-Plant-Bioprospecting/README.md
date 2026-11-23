# Medicinal Plant Bioprospecting for Anticancer Leads 🌿

A computational research-learning project exploring how phytochemical information, cheminformatics and molecular docking can help prioritize plant-derived compounds for further investigation.

> **Status:** Exploratory research workflow. Candidate compounds and computational predictions must be validated experimentally; no anticancer efficacy is claimed.

## Research question
How can publicly available phytochemical and protein information be combined to prioritize plant-derived compounds for follow-up studies against cancer-relevant molecular targets?

## Proposed workflow
1. **Literature and plant selection** — define the plant species, plant part, extraction context and evidence for reported constituents.
2. **Compound curation** — collect compound names, identifiers, structures and source references from traceable databases or publications.
3. **Structure standardization** — inspect salts, duplicates, stereochemistry, protonation and canonical identifiers.
4. **Cheminformatics** — calculate molecular descriptors and perform exploratory similarity or clustering analyses.
5. **Target selection** — justify cancer-relevant protein targets using literature and biological context.
6. **Molecular docking** — prepare protein and ligand structures, document parameters and inspect poses.
7. **Prioritization** — combine evidence transparently; distinguish literature evidence, computational predictions and experimental observations.
8. **Validation plan** — propose appropriate biochemical or cell-based assays for future work.

## Tools
- Python and pandas for data handling.
- RDKit for molecular structure handling and descriptors.
- AutoDock / AutoDockTools or PyRx for docking workflows.
- XGBoost only where a suitable, labelled dataset and evaluation strategy are available.

## Suggested structure
```
Medicinal-Plant-Bioprospecting/
├── README.md
├── data/
│   └── README.md
├── notebooks/
├── src/
├── results/
└── references/
    └── sources.md
```

## Data provenance
For every compound, preserve:
- Plant species and plant part.
- Reported extraction or isolation method.
- Compound name and database identifier.
- Source publication or database record.
- Whether identity was experimentally confirmed in the cited study.
- Structure file provenance and any standardization performed.

Do not invent missing compound concentrations, assay results or target interactions. Keep raw source data separate from processed data and document transformations.

## Evaluation principles
- Avoid training predictive models on tiny or biased datasets.
- Use suitable train/test separation and cross-validation.
- Prevent data leakage between highly similar compounds.
- Report baseline performance and uncertainty.
- Treat docking scores and descriptor-based predictions as prioritization signals, not proof of activity.

## Roadmap
- [ ] Select and justify a plant species and plant part.
- [ ] Build a source-linked phytochemical dataset.
- [ ] Standardize structures and remove duplicates.
- [ ] Perform exploratory descriptor analysis.
- [ ] Select targets with literature-based justification.
- [ ] Run and document docking experiments.
- [ ] Develop an evidence table and validation plan.

## Author
**Vishesh Mankotia** · [GitHub](https://github.com/VisheshMankotia) · [LinkedIn](https://www.linkedin.com/in/vishesh-mankotia-b8569b292/)
