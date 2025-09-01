# Hashira Placements Assignment – Polynomial Secret Reconstruction

## Problem
Reconstruct the polynomial and find the secret `f(0)` from JSON input files containing base-encoded values.  
Method used: **Newton’s Divided Differences.**

## Files
- `solve.js` – Main program
- `sample.json` – First testcase
- `second.json` – Second testcase

## Requirements
- Node.js (v18 or higher)

## How to Run
Open a terminal inside the project folder:

```bash
# Run with the sample testcase
node solve.js sample.json

# Run with the second testcase
node solve.js second.json

# (Optional) Print full polynomial coefficients
node solve.js second.json --coeffs


Expected Output

Sample Testcase
Secret f(0): 3

Second Testcase
Secret f(0): -6290016743746469796