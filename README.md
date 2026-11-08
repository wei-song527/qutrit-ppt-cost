Title: Second-Level Termination of the PPT Entanglement Cost Hierarchy in Bipartite Three-Level Systems
Author: Wei Song，Xiao-Lan Zong, and Ming Yang
Material: Supplementary Material — Rational Certificates and Exact Verification

Description: This package contains 126 groups of rational certificates, as well as Python programs for exact verification of positive semidefinite lower bounds and parameter-domain coverage. In addition, it provides auxiliary verification programs for the boundary analytic proof, a verification guide, and the complete verification log.

1.Purpose and Relation to the Paper

Sections III–V of the paper give the shortening construction of the dual chain and the related matrix proof. Among them, Section IV introduces the finite certificates, and Section V gives the analytic proof of the boundary lemma. Appendix A explains the concrete steps of exact verification.

This package can be used to reproduce the rational computations in the finite-certificate part.
qutrit_ppt_cost_supplement.pdf explains the conventions used and the correspondence between the files and the content of the paper. The complete proof still needs to be combined with the assumptions and mathematical arguments in the paper. The programs passing verification does not mean that the entire theorem has been verified by a formal proof assistant.

2.Running Requirements

Python 3.9 or higher is required, and only the Python standard library is used.
No third-party packages, SDP solvers, or numerical computation libraries need to be installed, and no internet access or random seed is required.

Please extract the archive completely, and make sure that the data/ folder is in the same directory as the three Python files.

Do not use -O or -OO when running, and do not set PYTHONOPTIMIZE.
The assertions in the programs participate in verification, so the entry program refuses to run in optimized mode.

This version has been tested in a CPython 3.12.14 environment on Linux.
Other Python versions and operating systems have not been tested one by one.

3.Running Verification

Enter the extracted supplementary-material directory and run:
python3 run_verification.py

To perform the same complete audit as the attached log, run:
python3 run_verification.py --audit-all-units

On Windows, one may use:
py -3 run_verification.py
py -3 run_verification.py --audit-all-units

To save the log of this run separately while keeping the original reference log, run:
python3 run_verification.py --audit-all-units > verification_new.log 2>&1

It is also possible to run from another directory by using the full path to the script.
The programs locate input files relative to the script directory, not relative to the current working directory.

4.Expected Output and Error Handling

When verification succeeds, the program exits with code 0 and outputs:
HASHES 7 PASS
[item-by-item checks for the middle and dyadic regions]
[JSON summary of finite-region lower bounds and computation counts]
FINITE RATIONAL CERTIFICATES: PASS.
[boundary arithmetic, index, and sensitivity checks]
BOUNDARY SUPPORT AND SENSITIVITY CHECKS: PASS.
ALL EXACT CHECKS PASSED.

The complete audit checks 126 certificate groups, computes map values on 17010 upper-triangular matrix units and 13608 lower-triangular matrix units, totaling 30618 computations.
The default mode does not directly compute lower-triangular matrix units, but reconstructs the corresponding contributions using adjoint identities. Both modes yield the same residual bounds.

Two tests with intentionally introduced errors should output REJECTED.
This means that the program correctly identified the modified test data, not that the original certificates failed verification. These modifications are made only in memory and do not alter the files.

If an assertion failure, missing input file, checksum mismatch, or other exception occurs, the program terminates with a nonzero exit code, and this run cannot be regarded as passing verification.
If a missing file is reported, first check whether extraction was complete. If a hash value does not match, restore the original file; do not accept unverified changes by modifying the hash value.

5.File Description

README.TXT
Package description, running requirements, commands, and verification scope.

qutrit_ppt_cost_supplement.pdf
Verification guide, for reading and checking.

qutrit_ppt_cost_supplement.tex
Editable REVTeX 4.1 source file of the verification guide.

run_verification.py
Entry program, responsible for running the two verification modules.
If either module fails, the overall run also ends in failure.

verify_certificates.py
Rebuilds the rational Choi maps, computes coefficient errors and polynomial lower bounds,
and checks the norm conditions and the coverage of the finite parameter domain.

verify_boundary.py
Checks rational inequalities, partial-transpose index relations,
symmetric projector identities, and index permutations in the Frobenius norm estimates.
It separately checks whether the program correctly reports an error when a transpose is omitted or a covering leaf node is missing.

data/dyadic.json
Six certificates, covering 1/640 <= t <= 1/10.
The denominator is 10^10, and the required lower bound is 1/100000.

data/middle.json
119 certificates, with denominator 10^8, and the required lower bound is 1/5000.
The partition tree has 243 nodes in total, including 119 leaf nodes verified by certificates,
and three leaf nodes containing no feasible parameters.

data/upper.json
One certificate, covering [13/16,1]^3.
The denominator is 10^8, and the required lower bound is 1/100.

SHA256SUMS.json
Records the SHA-256 hash values of this README, the three programs, and the three data files.
The verification programs check the inputs against this fixed seven-file list.
The hash values are used only to check the byte content of the files, not to prove positive semidefiniteness.

verification.log
All standard output and standard error output produced by the complete audit of the released version.

The input file list does not include the list file itself, the PDF and source file of the guide, or verification.log.
Running verification does not require LaTeX. Only when compiling the guide source file are REVTeX 4.1, Latin Modern, and hyperref required.

6.Data Format and Mathematical Conventions

The JSON files store data as integer arrays and rational-number strings.
Certificate computations use fractions.Fraction, not binary floating point.

Each record in the middle and dyadic files contains the following fields:
id; box=[x_lower,x_upper,t_lower,t_upper]; beta; gram_denominator; maps;
scalar_gram_numerators; scalar_weight_numerators; sphere_numerators;
weight_denominator; mu.

Here alpha=1/2, s=(1+beta)/2.
The constraint on y is reconstructed from box and the parameter domain given in Section IV of the paper.

The middle file also contains root and tree.
Each tree node performs one of the following operations: splitting an x or t interval at the exact midpoint, specifying a certificate, or excluding an infeasible rectangular region.
The lower endpoint of the dyadic region is 1/640.

In the upper file, alpha=2/3 and beta=1/3, hence s=2/3.
Additional fields include amplitude_box and lower_polynomial_rows.
lower_polynomial_rows stores integer coefficients; dividing by 10000 gives the actual coefficients, in the order:
1, x^2, y^2, t^2, xy, xt, yt.

Each map stores, block by block, the Choi coordinate indices and the integer numerators of the Gram factor.
Dividing the integer numerators by gram_denominator gives L, and the corresponding Choi block is L L^T, with denominator the square of gram_denominator.
The index flattening is (p,b) -> 9*p+b.

All maps have output dimension nine.
The input dimension of S is three, and the input dimensions of the remaining maps are all nine.
AH,AHT,AE,AET,AR,ART in the files correspond respectively to A_H,A_HT,A_E,A_ET,A_R,A_RT in the paper.

All maps are extended complex-linearly to complex matrices.
The full transposes in AHT,AET,ART are retained.
When estimating errors, each pair of adjoint terms corresponding to an off-diagonal matrix unit is included, so the resulting bounds apply to complex Hermitian inputs.
Deviations between the rational coefficients and the target coefficients are explicitly estimated and subtracted from the positive lower bounds given by the certificates.

7.Coverage and Scope of Conclusions

The analytic proof covers 0 < t <= 1/640.
The finite certificates cover 1/640 <= t <= 1, and both parts include the common endpoint.

Section III of the paper directly handles the case ||C|| <= 2.
When ||C|| > 2, the remaining assumptions guarantee that t is nonzero.

The boundary module is used to check auxiliary arithmetic relations and identities.
Positive semidefiniteness in the general case, joint estimates for traces, and Schur complement arguments are all proved in the paper, and also cover the case of singular matrix blocks.

This paper proves that on C^3 ⊗ C^3, for all p >= 2, one has chi_p = chi_2.
The exact zero-error asymptotic PPT entanglement cost is log2(chi_2) = log2(kappa_3).

This conclusion does not determine the LOCC entanglement cost, nor does it prove that the cost equals log2(kappa_2).
It neither gives a uniform termination level applicable to all dimensions, nor states that every quantum state must be computed to the second level.

8.Version Record

The three certificate JSON files are byte-for-byte identical to the files attached to the previous manuscript version, and the mathematical verification algorithm is also unchanged.
This time, only the paper title and author information were added in the programs, and the explanatory file name in the input file list was unified as README.TXT.
There is no need to rerun numerical optimization or certificate search.

This version has re-extracted the seven input files and their list into a clean directory, and completed the full audit again.
The test environment is described in Section 2, and all run output is in verification.log.

This package provides the final rational certificates and their verification programs; it does not provide the complete workflow for rerunning the early numerical search.
