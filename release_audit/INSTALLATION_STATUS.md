# Installation status

Local candidate checks passed on Ubuntu 20.04.6 x86_64, Python 3.11.15, with public PyPI dependencies. Core 0.3.0 and wrapper 0.1.1rc1 were built as separate wheels and source distributions and installed together in two fresh standard virtual environments. No editable install, project PYTHONPATH or system-site-packages were used. Imports resolved inside each virtual environment. `pip check`, CLI help, synthetic H5AD generation, Both/Original/Expanded analysis, seed details and score-run replay passed outside the checkout.

The original regression suites passed 189 tests. The portable installed-package suite passed 78 tests after correcting one new harness cache-interface mistake; the initial 77-pass/1-fail record is retained. Mouse exact-output/source-integrity parity passed 2059 checks. Human numeric/NA/coverage parity passed 13 checks while explicitly recording a changed provenance digest. Fourteen additional CLI checks passed, including 11 malformed-input refusals.

The 256-profile × 795-gene synthetic example selects 240 target profiles across four sample records (no confirmed animals). The initial Both runs took about 20.6–20.8 seconds internally; process high-water RSS observed in the per-version exports reached about 0.41 GiB. Timing covers only this small implementation example. Exact environment/performance and 36 observed output-table schemas are recorded alongside this file.

Initial network-restricted bootstrap/build attempts failed and remain in the retained logs. Approved public-network attempts succeeded. No offline dependency closure is claimed. Final artifact hashes and any final artifact reinstalls are recorded in the internal release report rather than inferred from these initial checks.

GitHub Actions is configured for Ubuntu 24.04/Python 3.11 but has not run. Public repository/tag/asset downloads and anonymous installation have not run. Windows, macOS and million-cell scalability are unverified. The candidate remains staged pending public-license and Reference/mapping redistribution approval.
