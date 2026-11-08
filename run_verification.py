# Supplemental Material for
# Second-level termination of the PPT entanglement-cost hierarchy for two qutrits
# Wei Song; release 2026-11-07.
"""Reproduce the exact checks accompanying qutrit_ppt_cost_pra.tex."""
import argparse
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Run without -O: assertions are part of the verifier.')

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--audit-all-units',action='store_true',
                    help='also compute every lower matrix unit and compare row bounds')
args=parser.parse_args()
root=Path(__file__).resolve().parent
commands=[['verify_certificates.py']+(['--audit-all-units'] if args.audit_all_units else []),
          ['verify_boundary.py']]
for command in commands:
    print('CHECKING',' '.join(command),flush=True)
    subprocess.run([sys.executable,str(root/command[0]),*command[1:]],check=True)
print('ALL EXACT CHECKS PASSED.',flush=True)
print('The universal conclusion also uses the mathematical arguments in the manuscript.')
