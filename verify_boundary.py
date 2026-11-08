# Supplemental Material for
# Second-level termination of the PPT entanglement-cost hierarchy for two qutrits
# Wei Song; release 2026-11-07.
"""Exact arithmetic, indexing, and sensitivity checks for the written proof.

These checks support, and do not replace, its universal matrix arguments.
Python 3.9+, standard library only; no numerical or SDP dependencies.
"""
from fractions import Fraction as F
from pathlib import Path
import json, copy
import verify_certificates as a

def arithmetic():
    tau = F(1, 640)


    def positive(label, number):
        assert number > 0, (label, number)
        print(label + ": " + str(number) + " > 0")


    positive("trace budget implies T < 4t^2", 4 - (3 + 2*tau))
    positive("symmetric denominator exceeds one", (2 - 4*tau) - 1)
    assert F(1, 24) + F(1, 24) + F(11, 12) == 1
    positive("shared arm loss", 50*(2 - 4*tau) - 96)
    positive("b loss at most t^2", (2 - 4*tau) - 384*tau)
    positive("antisymmetric pivot", 2 - F(21, 10)*tau - 6*tau**2)
    positive("inverse polynomial upper coefficient", 9 - F(4113, 500) - F(159, 25)*tau)
    positive("inverse polynomial remaining factor", F(1, 50) - 9*tau)
    positive("flag constant", F(3, 220) - F(12, 11)*tau)
    positive("flag trace coefficient", F(1, 22) - F(33, 40)*tau)
    positive("positive flag allowance on actual trace budget", F(33, 40) - (3+2*tau)/4)
    positive("e < 1/3", F(1, 3) - 203*tau)
    positive("contradiction p < 2t", 2 - F(191, 100) - tau)
    positive("sqrt(t) < 1/25", F(1, 625) - tau)
    positive("linear deficit", F(48, 25) - tau - 2*tau**2 - F(191, 100))
    assert F(191, 100) - F(53, 50) - F(33, 40) == F(1, 40)
    assert F(1, 40) - F(21, 2)*tau == F(11, 1280)
    assert F(1, 4)/tau - 150 == 10
    positive("final margin above t/128", F(11, 1280) - F(1, 128))
    assert F(1, 10*2**6) == tau
    assert 6 + 119 + 1 == 126
    print("Closed analytic and finite regions meet exactly at t=1/640.")
    print("ALL SHORT-PROOF ARITHMETIC CHECKS PASSED.")

def indexing():
    # Direct four-index construction, distinct from the sparse index formula.
    def pt_dense(m,d):
        out=[[F(0)]*(d*d) for _ in range(d*d)]
        for i in range(d):
            for j in range(d):
                for k in range(d):
                    for l in range(d):out[d*i+l][d*k+j]=m[d*i+j][d*k+l]
        return out
    for p in range(9):
        for q in range(9):
            unit=[[F(0)]*9 for _ in range(9)];unit[p][q]=1
            got=pt_dense(unit,3)
            assert {(i,j):v for i,row in enumerate(got) for j,v in enumerate(row) if v}==a.gamma({(p,q):F(1)})
    sym=[[F(int(i==j),2) for j in range(4)] for i in range(4)]
    for i in range(2):
        for j in range(2):sym[2*i+j][2*j+i]+=F(1,2)
    target=[[F(int(i==j),2)+F(int(i in [0,3] and j in [0,3]),2) for j in range(4)] for i in range(4)]
    assert pt_dense(sym,2)==target
    old=[0,1,3,4];arms=[2,5,6,7];images=[]
    for i in old:
        for j in arms:
            p=3*(i//3)+j%3;q=3*(j//3)+i%3
            images.append((p,q) if p in old else (q,p))
    assert sorted(images)==[(p,q) for p in old for q in arms]
    print('Direct four-index PT on 81 units; symmetric trace identity; Frobenius permutation: PASS')

def sensitivity(root):
    data=json.loads((root/'data/middle.json').read_text());r=data['records'][0]
    true=a.regional(r);assert true['margin']>=F(1,5000)
    original=a.variable_inputs
    def missing_transpose(which,u):
        terms=original(which,u)
        if which=='H':terms['AHT']=u
        return terms
    try:
        a.variable_inputs=missing_transpose
        bad=a.regional(r)
    finally:a.variable_inputs=original
    assert bad['margin']<0
    print('Deliberately omitted complex transpose: REJECTED; resulting margin',bad['margin'])
    rows={r['id']:{'box':list(map(F,r['box']))} for r in data['records']}
    a.tree_check(data,rows)
    corrupt=copy.deepcopy(data)
    leaf=next(k for k,v in corrupt['tree'].items() if 'certificate' in v)
    del corrupt['tree'][leaf]
    try:a.tree_check(corrupt,rows)
    except AssertionError:print('Deliberately missing coverage leaf: REJECTED')
    else:raise AssertionError('Coverage corruption was not detected')

def run(root):
    arithmetic()
    indexing()
    sensitivity(root)
    print('BOUNDARY SUPPORT AND SENSITIVITY CHECKS: PASS.',flush=True)

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Run without -O.')
    run(Path(__file__).resolve().parent)
