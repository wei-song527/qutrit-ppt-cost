# Supplemental Material for
# Second-level termination of the PPT entanglement-cost hierarchy for two qutrits
# Wei Song; release 2026-11-07.
"""Exact rational verifier for the 126 two-qutrit certificate collections.

Python 3.9+, standard library only. Every bound is reconstructed from the data.
Upper-triangular matrix units include both adjoint contributions. The optional
full-unit audit independently computes the lower units and compares all rows.
"""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict

AUDIT_ALL_UNITS=False
EVALUATIONS={'upper_units':0,'lower_audit_units':0}

def plus(*terms):
    d=defaultdict(Q)
    for factor,matrix in terms:
        for ij,v in matrix.items():d[ij]+=factor*v
    return {ij:v for ij,v in d.items() if v}
def transpose(m):return {(j,i):v for (i,j),v in m.items()}
def gamma(m):return {(3*(i//3)+j%3,3*(j//3)+i%3):v for (i,j),v in m.items()}
I={(i,i):Q(1) for i in range(9)}
cols=[{0:1},{4:1},{8:1}]
for sign in [1,-1]:
    for i,j in [(0,1),(0,2),(1,2)]:cols.append({3*i+j:1,3*j+i:sign})
def congr(m):
    return {(a,b):v for a in range(9) for b in range(9)
            if (v:=sum(s*t*m.get((i,j),0) for i,s in cols[a].items() for j,t in cols[b].items()))}
N=congr(I);norms=[N[i,i] for i in range(9)]
def skew(m):return {(i-6,j-6):v for (i,j),v in m.items() if i>=6 and j>=6}

def decode(record):
    den=record['gram_denominator'];assert den in [10**8,10**10]
    maps={}
    for name,d in record['maps'].items():
        dim=3 if name=='S' else 9
        assert d['input_dimension']==dim and d['output_dimension']==9
        op=defaultdict(dict);covered=set()
        for block in d['blocks']:
            ids=block['indices'];g=block['gram_numerators']
            assert ids and len(ids)==len(g) and len(set(ids))==len(ids)
            assert not covered.intersection(ids);covered.update(ids)
            assert all(type(i)==int and 0<=i<9*dim for i in ids)
            width=len(g[0]);assert width>0 and all(len(row)==width for row in g)
            assert all(type(x)==int for row in g for x in row)
            for a,ia in enumerate(ids):
                for b,ib in enumerate(ids):
                    value=sum(g[a][k]*g[b][k] for k in range(width))
                    if value:op[ia//9,ib//9][ia%9,ib%9]=Q(value,den**2)
        assert covered==set(range(9*dim));maps[name]=dict(op)
    assert set(maps)==set('AH AHT AE AET AR ART U Q B V S'.split())
    return maps
def apply(op,m):
    return plus(*[(v,op.get(ij,{})) for ij,v in m.items()])

def core(v):
    m={(4*i,4*i):v[i] for i in range(3)}
    for value,(i,j) in zip(v[3:],[(0,1),(0,2),(1,2)]):
        m[4*i,4*j]=m[4*j,4*i]=value
    return {ij:q for ij,q in m.items() if q}
def constant_inputs(v):
    c=core(v);g=congr(plus((1,I),(-1,gamma(c))))
    return {'U':plus((2,I)),'Q':g,'B':g,'V':plus((2,I),(-1,c)),
            'S':{(i,i):2*(1-v[3+i]) for i in range(3)}}
def variable_inputs(which,unit):
    z=congr(gamma(unit));zz=congr(unit)
    if which=='H':return {'AH':unit,'AHT':transpose(unit),'Q':plus((-1,z)),'B':z,'S':skew(z)}
    if which=='E':return {'AE':unit,'AET':transpose(unit),'Q':plus((-1,z)),'B':plus((-1,z)),'V':plus((-1,unit)),'S':skew(z)}
    return {'AR':unit,'ART':transpose(unit),'U':plus((-1,unit)),
            'Q':plus((-1,zz)),'B':plus((-1,zz)),'V':plus((-2,gamma(unit)))}
def mapped(maps,inputs):return plus(*[(1,apply(maps[k],v)) for k,v in inputs.items()])
def reconstruct(record,alpha):
    maps=decode(record);sums=[Q(0)]*9;full_sums=[Q(0)]*9
    for which in 'HER':
        factor=1 if which=='H' else (1-2*alpha if which=='E' else 0)
        def error_at(i,j):
            u={(i,j):Q(1)}
            return plus((1,mapped(maps,variable_inputs(which,u))),
                        (-factor,congr(gamma(u))))
        for i in range(9):
            for j in range(i,9):
                error=error_at(i,j);EVALUATIONS['upper_units']+=1
                if i==j:assert error==transpose(error)
                for (a,b),v in error.items():
                    sums[a]+=abs(v)
                    if i<j:sums[b]+=abs(v)
                    if AUDIT_ALL_UNITS:full_sums[a]+=abs(v)
                if AUDIT_ALL_UNITS and i<j:
                    lower=error_at(j,i);EVALUATIONS['lower_audit_units']+=1
                    assert lower==transpose(error),(which,i,j)
                    for (a,b),v in lower.items():full_sums[a]+=abs(v)
    if AUDIT_ALL_UNITS:assert sums==full_sums
    delta=max(sums)
    beta=Q(record['beta'])
    def residual(v):return plus((1,N),(-beta,congr(gamma(core(v)))),(-1,mapped(maps,constant_inputs(v))))
    base=residual([Q(0)]*6);coeff=[base]
    for k in range(6):
        v=[Q(0)]*6;v[k]=1;coeff.append(plus((1,residual(v)),(-1,base)))
    eta=max(sum(abs(m.get((i,j),0)) for m in coeff for j in range(9) if j!=i)/norms[i] for i in range(9))
    return coeff,delta,eta

# Polynomials are dictionaries of (x,y,t) exponent triples.
one={(0,0,0):Q(1)};variables=[{tuple(int(i==j) for i in range(3)):Q(1)} for j in range(3)]
def multiply(a,b):
    out=defaultdict(Q)
    for e,u in a.items():
        for f,v in b.items():out[tuple(i+j for i,j in zip(e,f))]+=u*v
    return dict(out)
features=[one]+[multiply(variables[i],variables[j]) for i,j in [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]]
def scalar_basis(box):
    xl,xu,tl,tu=box;yl=max(tl,1-tu*tu);g=[]
    for var,l,u in zip(variables,[xl,yl,tl],[xu,xu,tu]):g.extend([plus((1,var),(-l,one)),plus((u,one),(-1,var))])
    g.extend([plus((1,variables[0]),(-1,variables[1])),plus((1,variables[1]),(-1,variables[2]))])
    return [one]+g+[multiply(g[i],g[j]) for i in range(8) for j in range(i,8)]
def regional(record):
    coeff,delta,eta=reconstruct(record,Q(1,2));box=list(map(Q,record['box']));xl,xu,tl,tu=box
    assert 0<=xl<xu<=1 and 0<tl<tu<=1
    beta=Q(record['beta']);assert beta==1-(tl+tu)/2 and 0<=beta<=1
    assert (1+beta)*(2*xu*xu+tu*tu)<=4
    basis=scalar_basis(box);den=record['weight_denominator'];gd=record['gram_denominator']
    assert den==gd
    weights=record['scalar_weight_numerators'];sphere=record['sphere_numerators'];grams=record['scalar_gram_numerators']
    assert len(weights)==45 and all(len(r)==9 for r in weights) and len(sphere)==len(grams)==9
    assert all(type(v)==int and v>=0 for r in weights for v in r)
    assert all(type(v)==int and v>=0 for v in sphere)
    mu=Q(record['mu']);errs=[]
    spherepoly=plus((-2,one),*[(1,p) for p in features[1:4]])
    z=[one]+variables
    for i in range(9):
        target=plus(*[(c.get((i,i),0)/norms[i],f) for c,f in zip(coeff,features)])
        g=grams[i];assert len(g)==4 and all(len(row)==4 for row in g)
        assert all(type(v)==int for row in g for v in row)
        poly=plus((mu,one),(Q(sphere[i],den),spherepoly),*[(Q(weights[k][i],den),b) for k,b in enumerate(basis)])
        for k in range(4):
            factor=plus(*[(Q(g[j][k],gd),z[j]) for j in range(4)])
            poly=plus((1,poly),(1,multiply(factor,factor)))
        errs.append(sum(abs(v) for v in plus((1,target),(-1,poly)).values()))
    bound=mu-max(errs)-delta-eta
    return {'box':box,'margin':bound,'delta':delta,'eta':eta,'polynomial_error':max(errs)}

def tree_check(data,records):
    root=list(map(Q,data['root']));assert root==[Q(13,16),Q(1),Q(1,10),Q(13,16)]
    tree=data['tree'];seen=set();used=set();discarded=[]
    stack=[('',root)]
    while stack:
        key,box=stack.pop();assert key in tree and key not in seen;seen.add(key);node=tree[key]
        l,u,a,b=box
        if set(node)=={'certificate'}:
            k=node['certificate'];assert k not in used and records[k]['box']==box;used.add(k)
        elif set(node)=={'discard'}:
            assert node['discard'] is True and 2*u*u+b*b<2;discarded.append(key)
        else:
            assert set(node)=={'split'} and node['split'] in ['x','t']
            if node['split']=='x':m=(l+u)/2;children=[[l,m,a,b],[m,u,a,b]]
            else:m=(a+b)/2;children=[[l,u,a,m],[l,u,m,b]]
            stack.extend([(key+str(i),c) for i,c in enumerate(children)])
    assert seen==set(tree) and used==set(records) and 3*Q(13,16)**2<2
    return len(seen),len(used),len(discarded)

def balanced(record):
    assert Q(record['alpha'])==Q(2,3) and Q(record['beta'])==Q(1,3)
    c,delta,eta=reconstruct(record,Q(2,3))
    rows=record['lower_polynomial_rows'];assert len(rows)==9 and all(len(row)==7 for row in rows)
    vals=[]
    for i,row in enumerate(rows):
        assert all(type(v)==int for v in row) and all(v>=0 for v in row[1:])
        lo=[Q(v,10000) for v in row]
        assert all(cc.get((i,i),0)/norms[i]>=v for cc,v in zip(c,lo))
        vals.append(lo[0]+sum(lo[1:])*Q(13,16)**2)
    margin=min(vals)-delta-eta;assert margin>=Q(1,100)
    return {'margin':margin,'delta':delta,'eta':eta}

def run(root,audit_all_units=False):
    global AUDIT_ALL_UNITS
    AUDIT_ALL_UNITS=audit_all_units
    EVALUATIONS.update(upper_units=0,lower_audit_units=0)
    results={};collections=0
    manifest=root/'SHA256SUMS.json';d=json.loads(manifest.read_text())
    for name,want in d.items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==want,name
    assert set(d)=={'README.TXT','run_verification.py','verify_certificates.py',
                    'verify_boundary.py','data/middle.json','data/dyadic.json','data/upper.json'}
    print('HASHES',len(d),'PASS',flush=True)
    for name,margin in [('data/middle.json',Q(1,5000)),('data/dyadic.json',Q(1,100000))]:
        data=json.loads((root/name).read_text());rows={}
        assert Q(data['uniform_margin'])==margin
        for r in data['records']:
            k=r['id'];assert type(k)==int and k not in rows
            result=regional(r);assert result['margin']>=margin,(name,k,result)
            rows[k]=result;collections+=1
            print(name,'region',k,'PASS',flush=True)
        if name=='data/middle.json':
            assert data['schema']=='qutrit-general-region-v1'
            assert len(rows)==119;results['tree_counts']=tree_check(data,rows)
        else:
            assert data['schema']=='qutrit-boundary-dyadic-v2'
            assert set(rows)==set(range(6)) and Q(data['threshold'])==Q(1,640)
            for k,r in rows.items():
                b=Q(1,10*2**k);assert r['box']==[1-b*b/2,Q(1),b/2,b]
            assert rows[5]['box'][2]==Q(1,640)
            assert rows[0]['box'][3]==Q(1,10)
            assert all(rows[k+1]['box'][3]==rows[k]['box'][2] for k in range(5))
        results[name]={'minimum_margin':str(min(r['margin'] for r in rows.values())),
                      'maximum_defect':str(max(r['delta'] for r in rows.values()))}
    upper=json.loads((root/'data/upper.json').read_text())
    assert upper['schema']=='qutrit-proportional-split-v1'
    results['upper']={k:str(v) for k,v in balanced(upper).items()}
    collections+=1;assert collections==126
    assert EVALUATIONS['upper_units']==126*3*45
    assert EVALUATIONS['lower_audit_units']==(126*3*36 if audit_all_units else 0)
    results['collections']=collections
    results['matrix_unit_evaluations']=dict(EVALUATIONS)
    results['full_units_reconstructed']=126*3*81
    results['all_adjoint_identities_and_row_sums_audited']=audit_all_units
    print(json.dumps(results,indent=2),flush=True)
    print('FINITE RATIONAL CERTIFICATES: PASS.',flush=True)

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Run without -O.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-all-units',action='store_true')
    args=parser.parse_args()
    run(Path(__file__).resolve().parent,args.audit_all_units)
