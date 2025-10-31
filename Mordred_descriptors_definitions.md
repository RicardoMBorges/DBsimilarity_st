
## **1. Constitutional Descriptors (Simple Composition)**

Quantify the elemental and compositional makeup of a molecule.
**Nature:** Zero-order, independent of molecular connectivity.

**Examples:**

* `nAtom`, `nHeavyAtom`, `nC`, `nO`, `nN`, `nF`, `nCl`, `nBr`, `nI`
* `nH`, `nB`, `nSi`, `nP`, `nS`, `nSe`, `nAs`, `nSn`, `nPb`
* `nRing`, `nARing`, `nFaRing`, `nFHRing`, `n5aRing`, `n7Ring`, etc.
* `nRot`, `nRotBond`, `nBO`, `nBonds`, `nDouble`, `nTriple`, `nAromBond`
* `AMW`, `MW`, `MolecularWeight`, `ATSm`, `nHetero`
* `nHBAcc`, `nHBDon`
* `nX`, `nY`, `nZ` (total counts per element)
* `FractionCSP3`, `nBridgehead`, `nSpiro`, `nFused`, `nSubstituents`

*Represents molecular size, elemental composition, unsaturation, hydrogen bonding, and structural complexity.*

---

## **2. Topological Descriptors (Graph-based 2D Connectivity)**

Derived from the connectivity matrix or distance matrix of the molecular graph.
**Nature:** Depend on how atoms are connected.

**Examples:**

* **Wiener indices:** `WienerIndex`, `WPath`, `WPol`
* **Zagreb indices:** `Zagreb1`, `Zagreb2`, `mZagreb1`, `mZagreb2`
* **Information content indices:** `IC1`–`IC5`, `TIC1`–`TIC5`, `SIC1`–`SIC5`, `BIC1`–`BIC5`
* **Chi indices (Kier–Hall):** `Chi0`, `Chi1`, `Chi1n`, `Chi1v`, `Chi2n`, `Chi2v`, `Chi3v`
* **Balaban & Randic:** `BalabanJ`, `Randic`, `SpAbs_Dz`, `SpMax_Dz`, `SpDiam_Dz`
* **Connectivity indices:** `Connectivity0`, `Connectivity1`, etc.
* **Harary, Schultz:** `Harary`, `Schultz`
* **Petitjean:** `PetitjeanNumber`, `PetitjeanShape`
* **Topological charge & mean:** `GGI1`–`GGI10`, `JGI1`–`JGI10`, `ATS*`, `GATS*`, `MATS*`
* **Edge adjacency:** `EEig*`, `VR*`, `DetourIndex`, `PathIndex`, `MDE*`, `MDEC*`

*Encodes branching, cyclicity, and connectivity complexity.*

---

## **3. Atom-Type E-State Descriptors (Electronic Topological State)**

Hybrid descriptors capturing both **topology and electronic state** per atom type.
**Nature:** Local environment + intrinsic electronic effect.

**Examples:**

* `SaaCH`, `SaaNH`, `SssNH`, `SdsCH`, `SdO`, `SssO`, `SsOH`, `SaaC`
* `SsBr`, `SsCl`, `SsI`, `SssP`, `SdS`, `SssS`, `SaaSe`, etc.
* `MAXDN`, `MAXDP`, `MINDN`, `MINDP`
* `NssCH2`, `NdsCH`, `NssO`, `NssNH`, `NaaCH`, `NsssCH`, `NddC`, etc.

*Reflects atom-level hybridization, substitution, and electronic polarization effects.*

---

## **4. E-State / Partial Charge–Surface Area Hybrid Descriptors**

Combine **electronic charge or E-state values** with **surface area partitions (VSA)**.
**Nature:** Property bins combining polarity and accessible area.

**Sub-families:**

* **EState_VSA*** — bins of E-state vs solvent-accessible area (EState_VSA1–11)
* **VSA_EState*** — inverse mapping (VSA_EState1–10)
* **PEOE_VSA*** — partial charge bins via PEOE algorithm (PEOE_VSA1–14)
* **SlogP_VSA*** — bins by hydrophobicity contribution (SlogP_VSA1–12)
* **SMR_VSA*** — bins by molar refractivity (SMR_VSA1–10)

*Estimate the balance of polar vs. non-polar regions and their size — crucial for ADME.*

---

## **5. Physicochemical Property Descriptors**

Directly estimate measurable properties from structure.

**Examples:**

* `LogS_FilterItLogS` (solubility)
* `SLogP`, `ALogP`, `CrippenLogP`
* `SMR` (molar refractivity)
* `McGowanVolume`
* `TopologicalPolarSurfaceArea`, `TopoPSA`
* `LabuteASA` (approx. molecular surface area)
* `apol`, `bpol` (atomic/bond polarizability)
* `VDWVolume`, `VDWArea`
* `MolarRefractivity`, `ExactMolWt`

*Correlate with solubility, lipophilicity, permeability, and polarizability.*

---

## **6. Shape and Geometrical Descriptors (3D)**

Depend on the 3D coordinates of a molecule.

**Examples:**

* `RadiusOfGyration`, `InertiaX/Y/Z`, `MomentOfInertia`, `Asphericity`, `Eccentricity`
* `PrincipalMomentRatio`, `Spherocity`, `PBF`, `PetitjeanShapeIndex`
* `MOMI-X/Y/Z`, `GeomPetitjeanIndex`
* `Vabc`, `MaxDepth`, `MinDepth`
* `LabuteASA`, `MolecularVolume`, `MolecularArea`

*Capture molecular compactness, elongation, anisotropy, and 3D size.*

---

## **7. 3D Distance-Matrix Derived Descriptors**

Numerical transforms of interatomic distance matrices, weighted by atomic properties.

### **a. MoRSE Descriptors**

*(Molecular Representation of Structures based on Electron Diffraction)*

* `Mor01m`–`Mor32m` (mass-weighted)
* `Mor01v`–`Mor32v` (volume-weighted)
* `Mor01p`–`Mor32p` (polarizability-weighted)
* `Mor01e`–`Mor32e` (electronegativity-weighted)
* `Mor01se`–`Mor32se` (sigma-electronegativity-weighted)

### **b. WHIM (Weighted Holistic Invariant Molecular)**

* `WHIM1`–`WHIM99` (various weightings: mass, polarizability, charge, volume)

### **c. GETAWAY (GEometry, Topology, and Atom-Weight AssemblY)**

* `HATS*`, `RCON*`, `RARS*`, `RMAX*`, `RMIN*`, etc.

*Describe 3D spatial distributions and are widely used in CoMFA-like 3D-QSAR.*

---

## **8. Functional Group and Fragment Counts**

Quantify structural motifs and pharmacophoric fragments.

**Examples:**

* `nHeteroRing`, `nFusedRing`, `nHBondDonors/Acceptors`
* `nAromRings`, `nARing`, `nFaRing`
* `nOHs`, `nCOOH`, `nCOO`, `nNH2`, `nCN`, `nNO2`, `nHalogen`
* `nROH`, `nROR`, `nRNH2`, `nR2NH`, `nR3N`, `nC=O`, etc.

*Encode substructural features responsible for key interactions.*

---

## **9. Index-Based Families (Advanced Topological Matrices)**

Mathematical transformations of adjacency and distance matrices.

**Examples:**

* **Burden eigenvalues:** `BEHv1`–`BEHv8`, `BELv1`–`BELv8`, etc.
* **Edge adjacency:** `SpAbs_EA`, `SpMax_EA`, `SpDiam_EA`
* **Eigenvalue sums/products:** `EEig0`, `EEig05`, `VR1`, `VR2`, etc.
* **Topological distance matrix:** `ATS*`, `GATS*`, `MATS*`
* **RDF (Radial Distribution Function):** `RDF010m`–`RDF150m`, `RDF010v`, `RDF010p`, etc.

*Describe how atomic properties are distributed in space or along the graph.*

---

## **10. Information-Theoretic Descriptors**

Derived from graph-theoretic entropy and complexity.

**Examples:**

* `IC1`–`IC5` (information content based on atom degrees)
* `TIC1`–`TIC5`, `SIC1`–`SIC5`, `BIC1`–`BIC5` (topological, structural, bonding IC)
* `MIC`, `ComplementaryIC`, `ComplexityIndex`

*Reflect molecular symmetry, redundancy, and complexity.*

---

## **11. Pharmacokinetic / Rule-Based Filters**

Simple Boolean or integer filters encoding drug-likeness.

**Examples:**

* `Lipinski`, `GhoseFilter`, `Veber`, `Egan`, `MDDRlikeRule`
* `FilterItLogS`, `FilterItLogP`
* `QED` (Quantitative Estimation of Drug-likeness)
* `TPSA`, `HBA`, `HBD`, `RotatableBonds`

*Used for filtering compound libraries or modeling ADMET.*

---

## **12. Indices of Shape and Complexity (Kier & Hall)**

Dimensionless indices describing molecular form.

**Examples:**

* `Kier1`, `Kier2`, `Kier3`
* `ShapeIndex`, `Flexibility`, `Eccentricity`, `Asphericity`, `InertialShapeFactor`

*Relate to steric bulk, flexibility, and linearity.*

---

## **13. Advanced Hybrid Descriptor Sets (Mordred-specific families)**

| Family                               | Meaning                                                        | Examples                                                                   |
| ------------------------------------ | -------------------------------------------------------------- | -------------------------------------------------------------------------- |
| **ETA (Extended Topochemical Atom)** | Combines topology, electronegativity, and atomic valence       | `ETA_eta_F`, `ETA_beta`, `ETA_dEpsilon_C`, `AETA_eta_B`, `AETA_dEpsilon_B` |
| **ATS/GATS/MATS**                    | 2D autocorrelation (unweighted, weighted by atomic properties) | `ATS1m`, `GATS2v`, `MATS3p`, etc.                                          |
| **RDF**                              | Radial Distribution Function — 3D property distribution        | `RDF020m`, `RDF120p`, `RDF100v`                                            |
| **3D autocorrelation (3D-MATS)**     | 3D analogs of MATS with spatial weighting                      | `3DMATS1m`, `3DMATS2v`                                                     |
| **GETAWAY**                          | Geometry + topology + atomic weights                           | `HATS2m`, `RCON5v`, `RMAX6p`                                               |
| **WHIM**                             | Weighted holistic molecular descriptors                        | `WHIM1`, `WHIM24v`, `WHIM38m`                                              |
| **MoRSE**                            | Electron diffraction–like signal transforms                    | `Mor03p`, `Mor12se`, `Mor27v`                                              |

*These encode multidimensional molecular information for robust QSAR/QSPR.*

---

## **14. Special Statistical or Derived Groups**

Mathematical manipulations of other descriptors.

**Examples:**

* `MAXDP`, `MINDP`, `MAXDN`, `MINDN` (max/min electron densities)
* `MeanMW`, `MeanEEig`, `VarEEig`, `SkewEEig`
* `VarianceWHIM`, `SumMorSE`

*Provide summary statistics over descriptor distributions.*

---

## **15. Identification / Index Descriptors**

Categorical identifiers are sometimes output by Mordred.

**Examples:**

* `MolecularId`, `MID_B`, `MID_C`, `MID_N`
* `FragCpx`, `FragIndex`, `FunctionalGroupCount`

*For bookkeeping or fragment tracking.*

---

## **Summary of Major Classes**

| **Class**                       | **Descriptor Families (Mordred)**        | **Typical Count Range** |
| ------------------------------- | ---------------------------------------- | ----------------------- |
| Constitutional                  | n*, RingCount, AMW, MW                   | 50–100                  |
| Topological                     | Chi, Balaban, Zagreb, WPath, IC/TIC      | 150–250                 |
| Electronic (E-State / AtomType) | SssCH, SdO, MAXDP, MINDP                 | 100–150                 |
| Charge/Surface Hybrid           | PEOE_VSA, SlogP_VSA, SMR_VSA, EState_VSA | 50–80                   |
| Physicochemical                 | logP, PSA, ASA, apol/bpol                | 20–40                   |
| Shape/3D                        | MOMI, PBF, Geom*, Vabc                   | 30–50                   |
| MoRSE/3D-field                  | Mor*, RDF*, WHIM, GETAWAY                | 300–500                 |
| Information                     | IC, TIC, SIC, BIC                        | 20–30                   |
| Functional Group                | nCOOH, nOH, nNH2, nArom                  | 40–60                   |
| Drug-likeness/Rules             | Lipinski, Ghose, Veber, QED              | 10–20                   |
| Others (ETA, Hybrid, Summary)   | ETA_*, MAXDP, MINDP                      | 100–200                 |

