MCI_CFC Project
===============

Abstract
--------

Brain structural (SC) and functional connectivities (FC) are intrinsically coupled, yet this relationship declines in mild cognitive impairment (MCI) and its mechanisms remain unclear. White matter hyperintensities (WMH) may further disrupt or reshape this coupling. Leveraging brain network communication models, this project quantifies how neural signals propagate along SC to explain CM-FC decoupling in MCI.

Background
----------

- Five communication measures (CM) spanning diffusion, path-accessibility, and routing strategies are derived from SC to characterize coupling profiles.
- The cohort comprises 186 individuals with MCI and 171 cognitively normal controls; comparisons are performed across diagnostic groups and tractography types.
- Lesion mapping plus mediation analyses evaluate how WMH burden alters coupling and impacts cognition.

Methods
-------

- Compute CM-FC coupling by correlating FC matrices with CM-derived communication patterns per subject and tractography variant.
- Evaluate diffusion-based (communicability, flow graph), path-accessibility, and routing-based strategies to determine which pathways drive decoupling.
- Perform lesion mapping focused on WMH-impaired tracts (e.g., corona radiata) and conduct mediation analyses to test whether coupling mediates WMH–memory relationships.

Results Summary
---------------

- MCI participants show significant whole-brain CM-FC decoupling, especially within default mode, sensorimotor, and visual networks.
- Diffusion-based strategies exhibit the strongest decoupling, whereas routing-based coupling is relatively preserved.
- WMH-impaired networks have reduced diffusion-based coupling but relatively enhanced path-accessibility coupling, implying compensatory rerouting.
- WMH in the corona radiata disrupts CM-FC coupling; diffusion-based coupling partially mediates the link between WMH burden and memory deficits.

Conclusion
----------

Diffusion-based CM-FC decoupling characterizes MCI and offers a framework for understanding how WMH lesions alter large-scale communication. These findings highlight network-level vulnerabilities and compensatory pathways relevant to white matter injury.

Overview
--------

- End-to-end pipeline for lesion quantification, connectome construction, and coupling analyses.
- Contains preprocessing pipelines (`proc_mat/`), analysis notebooks/scripts (`analysis/`), and WMH segmentation utilities (`WMHseg/`).

Prerequisites
-------------

- Python 3.9+ with scientific stack (numpy, pandas, nibabel, brainspace, pingouin, neuroCombat, netneurotools, etc.).
- MATLAB/FSL/ANTs (optional) for specific preprocessing scripts referenced in `svd_*.sh`.
- SLURM-capable cluster if running the provided `sbatch-*.sh` batch scripts.

Quick Start
-----------

1. Create and activate a virtual environment.
2. Install Python dependencies (pip, conda, or requirements file once defined).
3. Set environment variables if needed:
   - `MCI_CFC_HOME` to override the default project root.
   - `MCI_CFC_N_JOBS` to control multiprocessing workers.
4. Prepare data folders under `data/` (subject CSVs, connectomes, lesion maps, etc.).
5. Run processing scripts, e.g.:
   ```bash
   python analysis/corr_glo_lm_cfc.py
   python analysis/corr_net_lm_cfc.py
   ```

Repository Highlights
---------------------

- `analysis/`: statistical models, plotting, and subject selection helpers.
- `proc_mat/`: structural/functional matrix preprocessing, harmonization, QC.
- `WMHseg/`: shell wrappers around LST / SynthSeg lesion workflows.
- `preprocessing/`: preprocessing for BOLD and DTI.
