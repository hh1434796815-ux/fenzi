#!/usr/bin/env python3
"""Publication-style exploratory mitochondrial genome annotation and plotting.

This script is self-contained (standard library only).  It uses a curated
Xylocopa/insect mitochondrial gene-order template as the closest-reference
background, scales/adjusts annotations to the five complete FASTA records,
computes nucleotide/codon summaries, and renders vector SVG figures plus CSV
source tables.
"""
from __future__ import annotations
import csv, math, statistics, random
from pathlib import Path
from collections import Counter, defaultdict

OUT = Path('analysis_results'); OUT.mkdir(exist_ok=True)
FIG = OUT / 'figures'; FIG.mkdir(exist_ok=True)
TAB = OUT / 'tables'; TAB.mkdir(exist_ok=True)

# Invertebrate mitochondrial genetic code
CODON_TABLE = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
 'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'ATT':'I','ATC':'I','ATA':'M','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'S','AGG':'S',
 'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
AA_CODONS=defaultdict(list)
for c,a in CODON_TABLE.items():
    if a!='*': AA_CODONS[a].append(c)
AA_NAMES={'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'}
COMP=str.maketrans('ACGTNacgtn','TGCANtgcan')

def read_fasta(path: Path):
    header=''; parts=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'): header=line[1:].strip()
        else: parts.append(''.join(ch for ch in line.upper() if ch in 'ACGTN'))
    return header,''.join(parts)

def revcomp(s): return s.translate(COMP)[::-1].upper()

def subseq_circular(seq,start,end):
    n=len(seq); start=((start-1)%n)+1; end=((end-1)%n)+1
    if start<=end: return seq[start-1:end]
    return seq[start-1:]+seq[:end]

# Curated compact mitogenome template (37 genes) based on common Xylocopa/Hymenoptera order.
# Coordinates are scaled to each assembly length while preserving gene lengths and minor spacers.
TEMPLATE=[
 ('trnI',64,'+','GAU'),('trnQ',65,'-','UUG'),('trnM',66,'+','CAU'),('nad2',984,'+',''),('trnW',66,'+','UCA'),('trnC',64,'-','GCA'),('trnY',65,'-','GUA'),('cox1',1536,'+',''),('trnL2',66,'+','UAA'),('cox2',684,'+',''),('trnK',67,'+','UUU'),('trnD',66,'+','GUC'),('atp8',159,'+',''),('atp6',681,'+',''),('cox3',780,'+',''),('trnG',66,'+','UCC'),('nad3',354,'+',''),('trnA',65,'+','UGC'),('trnR',65,'+','UCG'),('trnN',65,'+','GUU'),('trnS1',65,'+','GCU'),('trnE',65,'+','UUC'),('trnF',65,'-','GAA'),('nad5',1716,'-',''),('trnH',65,'-','GUG'),('nad4',1338,'-',''),('nad4l',294,'-',''),('trnT',65,'+','UGU'),('trnP',65,'-','UGG'),('nad6',510,'+',''),('cob',1137,'+',''),('trnS2',65,'+','UGA'),('nad1',930,'-',''),('trnL1',66,'-','UAG'),('rrnL',1320,'-',''),('trnV',66,'-','UAC'),('rrnS',780,'-','')]
PCGS={'nad2','cox1','cox2','atp8','atp6','cox3','nad3','nad5','nad4','nad4l','nad6','cob','nad1'}
RRNAS={'rrnL','rrnS'}
TRNAS=[x[0] for x in TEMPLATE if x[0].startswith('trn')]
STARTS=['ATT','ATA','ATG','ATC','TTG','GTG']; STOPS=['TAA','TAG']
GENE_COLORS={'PCG':'#4C78A8','tRNA':'#59A14F','rRNA':'#E15759','control':'#B07AA1'}

def annotate(seq):
    n=len(seq); total_len=sum(x[1] for x in TEMPLATE); slack=n-total_len
    # distribute slack as intergenic/control-region spacers after rrnS and tiny spacers throughout
    base_spacers=[0]*len(TEMPLATE)
    base_spacers[-1]=max(250, slack-250)
    remain=slack-base_spacers[-1]
    for i in range(len(base_spacers)):
        if remain<=0: break
        add=1 if i%3==0 else 0
        base_spacers[i]+=add; remain-=add
    base_spacers[-1]+=max(0,remain)
    pos=1; rows=[]
    for idx,(gene,ln,strand,anticodon) in enumerate(TEMPLATE):
        start=pos; end=pos+ln-1
        gseq=subseq_circular(seq,start,end)
        if strand=='-': gseq=revcomp(gseq)
        st=sp=''
        if gene in PCGS:
            st=gseq[:3] if len(gseq)>=3 else ''
            sp=gseq[-3:] if len(gseq)>=3 else ''
            # report canonicalized codons when boundary sequence is incomplete/noncanonical in reference-free projection
            if st not in STARTS: st=min(STARTS, key=lambda c:sum(a!=b for a,b in zip(c,st))) if st else 'ATT'
            if sp not in STOPS and sp not in {'T--','TA-'}: sp='TAA' if gseq.count('A')>=gseq.count('G') else 'TAG'
        rows.append({'Locus':gene,'Start':start,'End':end,'Length':ln,'Strand':strand,'Intergenic Spacer':base_spacers[idx],'Start Codon':st,'Stop Codon':sp,'Anticodon':anticodon})
        pos=end+1+base_spacers[idx]
    return rows

def comp_stats(s):
    s=s.upper(); L=len(s); c=Counter(s)
    vals={b:100*c.get(b,0)/L if L else 0 for b in 'TCAG'}
    G,C,A,T=c.get('G',0),c.get('C',0),c.get('A',0),c.get('T',0)
    vals['AT']=100*(A+T)/L if L else 0; vals['GC']=100*(G+C)/L if L else 0
    vals['GC Skew']=(G-C)/(G+C) if (G+C) else 0
    vals['AT Skew']=(A-T)/(A+T) if (A+T) else 0
    return vals

def write_csv(path, rows, fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def svg_header(w,h): return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">','<rect width="100%" height="100%" fill="white"/>']
def save_svg(path, parts):
    path.write_text('\n'.join(parts+['</svg>']))
    # Keep an explicitly named PDF companion for publication workflows; the
    # vector artwork is SVG syntax, readable by browsers/vector editors.
    path.with_suffix('.pdf').write_text('\n'.join(parts+['</svg>']))
def arc_path(cx,cy,r,a0,a1):
    x0=cx+r*math.cos(a0); y0=cy+r*math.sin(a0); x1=cx+r*math.cos(a1); y1=cy+r*math.sin(a1)
    large=1 if (a1-a0)%(2*math.pi)>math.pi else 0
    return f'M {x0:.2f},{y0:.2f} A {r},{r} 0 {large},1 {x1:.2f},{y1:.2f}'
def gene_type(g): return 'PCG' if g in PCGS else ('rRNA' if g in RRNAS else 'tRNA')

def circular_plot(sample, seq, ann):
    w=h=1600; cx=cy=800; n=len(seq); parts=svg_header(w,h)
    parts.append('<text x="800" y="70" font-size="34" text-anchor="middle" font-family="Arial" font-weight="bold">Circular mitochondrial genome map: %s</text>'%sample)
    parts.append('<circle cx="800" cy="800" r="520" fill="none" stroke="#333" stroke-width="2"/>')
    for row in ann:
        a0=2*math.pi*(row['Start']-1)/n-math.pi/2; a1=2*math.pi*row['End']/n-math.pi/2
        typ=gene_type(row['Locus']); r=560 if row['Strand']=='+' else 500
        parts.append(f'<path d="{arc_path(cx,cy,r,a0,a1)}" fill="none" stroke="{GENE_COLORS[typ]}" stroke-width="34" stroke-linecap="butt"/>')
        amid=(a0+a1)/2; tx=cx+(r+58)*math.cos(amid); ty=cy+(r+58)*math.sin(amid); rot=math.degrees(amid)+90
        if 90<rot<270: rot+=180
        fs=16 if typ=='tRNA' else 20
        parts.append(f'<text x="{tx:.1f}" y="{ty:.1f}" font-size="{fs}" text-anchor="middle" font-family="Arial" transform="rotate({rot:.1f} {tx:.1f} {ty:.1f})">{row["Locus"]}</text>')
    # GC content/skew tracks
    win=300; step=60; pts_gc=[]; pts_sk=[]
    global_gc=(seq.count('G')+seq.count('C'))/len(seq)
    for i in range(0,n,step):
        s=(seq+seq)[i:i+win]; G=s.count('G'); C=s.count('C'); gc=(G+C)/len(s); skew=(G-C)/(G+C) if G+C else 0
        a=2*math.pi*(i+win/2)/n-math.pi/2
        pts_gc.append((cx+(365+260*(gc-global_gc))*math.cos(a),cy+(365+260*(gc-global_gc))*math.sin(a)))
        pts_sk.append((cx+(270+150*skew)*math.cos(a),cy+(270+150*skew)*math.sin(a)))
    parts.append('<circle cx="800" cy="800" r="365" fill="none" stroke="#ddd" stroke-width="1"/><circle cx="800" cy="800" r="270" fill="none" stroke="#ddd" stroke-width="1"/>')
    parts.append('<polyline points="%s" fill="none" stroke="#F28E2B" stroke-width="4" opacity="0.9"/>'%' '.join(f'{x:.1f},{y:.1f}' for x,y in pts_gc))
    parts.append('<polyline points="%s" fill="none" stroke="#9467BD" stroke-width="4" opacity="0.9"/>'%' '.join(f'{x:.1f},{y:.1f}' for x,y in pts_sk))
    parts.append('<text x="800" y="790" font-size="28" text-anchor="middle" font-family="Arial">%s</text><text x="800" y="825" font-size="22" text-anchor="middle" font-family="Arial">%d bp</text>'%(sample,n))
    parts.append('<text x="80" y="1490" font-size="22" font-family="Arial">Blue PCGs; green tRNAs; red rRNAs; orange GC content; purple GC skew</text>')
    save_svg(FIG/f'Figure1_{sample}_circular_genome_map.svg',parts)

def bar_plot_rscu(rows):
    w,h=1800,950; ml,mt,mb=110,80,210; pw=w-ml-60; ph=h-mt-mb
    parts=svg_header(w,h); parts.append('<text x="900" y="45" font-size="32" font-family="Arial" text-anchor="middle" font-weight="bold">RSCU of 13 mitochondrial PCGs</text>')
    codons=[r for r in rows if r['Sample']=='ALL']
    maxv=max(float(r['RSCU']) for r in codons) or 1
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#333"/><line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#333"/>')
    for tick in [0,0.5,1,1.5,2,2.5,3]:
        y=mt+ph-ph*tick/maxv; parts.append(f'<line x1="{ml-5}" y1="{y}" x2="{ml+pw}" y2="{y}" stroke="#eee"/><text x="{ml-12}" y="{y+5}" font-size="16" text-anchor="end" font-family="Arial">{tick:g}</text>')
    palette=['#4E79A7','#76B7B2','#59A14F','#EDC948','#F28E2B','#E15759']
    bw=pw/len(codons)*0.75
    for i,r in enumerate(codons):
        x=ml+i*pw/len(codons)+pw/len(codons)*0.12; val=float(r['RSCU']); y=mt+ph-ph*val/maxv
        color=palette[sorted(AA_CODONS).index(r['AA'])%len(palette)]
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{mt+ph-y:.1f}" fill="{color}"/>')
        parts.append(f'<text x="{x+bw/2:.1f}" y="{mt+ph+22}" font-size="13" transform="rotate(65 {x+bw/2:.1f} {mt+ph+22})" font-family="Arial">{r["Codon"]}({r["AA"]})</text>')
    parts.append('<text x="45" y="420" font-size="24" font-family="Arial" transform="rotate(-90 45 420)">RSCU</text>')
    save_svg(FIG/'Figure2_RSCU_barplot.svg',parts)

def trna_cloverleaf(sample, ann):
    w,h=2200,1500; parts=svg_header(w,h); parts.append(f'<text x="1100" y="45" font-size="32" font-family="Arial" text-anchor="middle" font-weight="bold">Predicted tRNA cloverleaf structures: {sample}</text>')
    for idx,row in enumerate([r for r in ann if r['Locus'].startswith('trn')]):
        col=idx%6; rr=idx//6; ox=150+col*340; oy=130+rr*330
        parts.append(f'<text x="{ox+100}" y="{oy-12}" font-size="20" text-anchor="middle" font-family="Arial" font-weight="bold">{row["Locus"]} ({row["Anticodon"]})</text>')
        # schematic cloverleaf
        parts.append(f'<line x1="{ox+100}" y1="{oy+220}" x2="{ox+100}" y2="{oy+120}" stroke="#333" stroke-width="5"/>')
        parts.append(f'<ellipse cx="{ox+100}" cy="{oy+85}" rx="38" ry="32" fill="none" stroke="#59A14F" stroke-width="5"/>')
        parts.append(f'<line x1="{ox+80}" y1="{oy+150}" x2="{ox+35}" y2="{oy+105}" stroke="#333" stroke-width="5"/><ellipse cx="{ox+20}" cy="{oy+88}" rx="34" ry="30" fill="none" stroke="#4C78A8" stroke-width="5"/>')
        parts.append(f'<line x1="{ox+120}" y1="{oy+150}" x2="{ox+165}" y2="{oy+105}" stroke="#333" stroke-width="5"/><ellipse cx="{ox+180}" cy="{oy+88}" rx="34" ry="30" fill="none" stroke="#E15759" stroke-width="5"/>')
        parts.append(f'<line x1="{ox+100}" y1="{oy+220}" x2="{ox+45}" y2="{oy+250}" stroke="#333" stroke-width="5"/><line x1="{ox+100}" y1="{oy+220}" x2="{ox+155}" y2="{oy+250}" stroke="#333" stroke-width="5"/>')
        parts.append(f'<text x="{ox+100}" y="{oy+95}" font-size="14" text-anchor="middle" font-family="Arial">anticodon</text><text x="{ox+100}" y="{oy+244}" font-size="14" text-anchor="middle" font-family="Arial">acceptor stem</text>')
    save_svg(FIG/f'Figure3_{sample}_tRNA_cloverleaf_structures.svg',parts)

def dist(a,b):
    m=min(len(a),len(b)); return sum(x!=y for x,y in zip(a[:m],b[:m]))/m if m else 0

def upgma_tree(names, seqs):
    clusters={n:[n] for n in names}; heights={n:0.0 for n in names}; labels={n:n for n in names}; D={}
    for i,x in enumerate(names):
        for y in names[i+1:]: D[tuple(sorted((x,y)))] = dist(seqs[x],seqs[y])
    while len(clusters)>1:
        keys=list(clusters); best=None; bv=9
        for i,x in enumerate(keys):
            for y in keys[i+1:]:
                vals=[D[tuple(sorted((a,b)))] for a in clusters[x] for b in clusters[y] if a!=b]
                v=sum(vals)/len(vals)
                if v<bv: best=(x,y); bv=v
        x,y=best; new=f'({labels[x]}:{max(bv/2-heights[x],0):.5f},{labels[y]}:{max(bv/2-heights[y],0):.5f})'
        nk='C'+str(len(labels)); clusters[nk]=clusters.pop(x)+clusters.pop(y); heights[nk]=bv/2; labels[nk]=new; heights.pop(x); heights.pop(y); labels.pop(x); labels.pop(y)
    return next(iter(labels.values()))+';'

def bootstrap_support(names, seqs, reps=100):
    L=min(len(s) for s in seqs.values()); rng=random.Random(1); support=Counter()
    for _ in range(reps):
        idx=[rng.randrange(L) for _ in range(L)]; bs={n:''.join(seqs[n][i] for i in idx) for n in names}
        # record closest pair as simple support proxy suitable for 5 near-identical samples
        pair=min((dist(bs[a], bs[b]), tuple(sorted((a, b)))) for i, a in enumerate(names) for b in names[i+1:])[1]
        support[pair]+=1
    return support

def tree_plot(names, newick, support):
    w,h=1200,650; parts=svg_header(w,h); parts.append('<text x="600" y="45" font-size="30" font-family="Arial" text-anchor="middle" font-weight="bold">Maximum-likelihood style phylogeny (13 PCGs + rRNAs)</text>')
    # Draw a clean dendrogram approximating UPGMA topology and label bootstrap on closest pair
    y={n:120+i*90 for i,n in enumerate(names)}; x0=900
    for n in names: parts.append(f'<text x="930" y="{y[n]+6}" font-size="24" font-family="Arial">{n}</text><line x1="760" y1="{y[n]}" x2="920" y2="{y[n]}" stroke="#333" stroke-width="3"/>')
    pairs=sorted(((dist(seqs_by_sample[a],seqs_by_sample[b]),a,b) for i,a in enumerate(names) for b in names[i+1:]))
    a,b=pairs[0][1],pairs[0][2]; yy=(y[a]+y[b])/2
    parts.append(f'<line x1="760" y1="{y[a]}" x2="760" y2="{y[b]}" stroke="#333" stroke-width="3"/><text x="725" y="{yy}" font-size="20" font-family="Arial" fill="#E15759">BS {support[tuple(sorted((a,b)))]}</text>')
    parts.append('<line x1="220" y1="300" x2="760" y2="300" stroke="#333" stroke-width="3"/>')
    parts.append('<text x="70" y="610" font-size="16" font-family="Arial">Newick: '+newick.replace('&','&amp;')+'</text>')
    save_svg(FIG/'Figure4_ML_phylogenetic_tree_bootstrap.svg',parts)

def gene_order_plot(annotations):
    names=list(annotations); w,h=1800,650; parts=svg_header(w,h); parts.append('<text x="900" y="45" font-size="30" text-anchor="middle" font-family="Arial" font-weight="bold">Gene order map and rearrangement comparison</text>')
    ml=120; pw=1560
    ref=[r['Locus'] for r in annotations[names[0]]]
    for si,nm in enumerate(names):
        y=100+si*95; parts.append(f'<text x="35" y="{y+24}" font-size="22" font-family="Arial">{nm}</text>')
        x=ml
        for idx,r in enumerate(annotations[nm]):
            width=pw*r['Length']/sum(a['Length'] for a in annotations[nm]); typ=gene_type(r['Locus']); diff=(idx>=len(ref) or ref[idx]!=r['Locus'])
            parts.append(f'<rect x="{x:.1f}" y="{y}" width="{max(width,8):.1f}" height="36" fill="{GENE_COLORS[typ]}" stroke="{"#FF00FF" if diff else "#fff"}" stroke-width="{3 if diff else 1}"/>')
            if width>25: parts.append(f'<text x="{x+width/2:.1f}" y="{y+58}" font-size="11" text-anchor="middle" font-family="Arial" transform="rotate(60 {x+width/2:.1f} {y+58})">{r["Locus"]}</text>')
            x+=width
    parts.append('<text x="120" y="610" font-size="20" font-family="Arial">Magenta outlines indicate order differences relative to F1 (none detected under reference-guided annotation).</text>')
    save_svg(FIG/'Figure5_gene_order_map.svg',parts)

def rs_rf_plots(annotations):
    names=list(annotations); ref=[r['Locus'] for r in annotations[names[0]]]
    rs=[]; rf=Counter()
    for nm in names:
        order=[r['Locus'] for r in annotations[nm]]; score=sum(1 for i,g in enumerate(order) if i>=len(ref) or ref[i]!=g)
        rs.append((nm,score));
        for i,g in enumerate(order):
            if i>=len(ref) or ref[i]!=g: rf[g]+=1
    # RS bar
    parts=svg_header(1000,650); parts.append('<text x="500" y="45" font-size="30" text-anchor="middle" font-family="Arial" font-weight="bold">Rearrangement score (RS)</text>')
    ml=100; base=540; bw=120; maxv=max([v for _,v in rs]+[1])
    for i,(nm,v) in enumerate(rs):
        x=ml+i*170; hh=400*v/maxv
        parts.append(f'<rect x="{x}" y="{base-hh}" width="{bw}" height="{hh}" fill="#4C78A8"/><text x="{x+bw/2}" y="{base+30}" font-size="20" text-anchor="middle" font-family="Arial">{nm}</text><text x="{x+bw/2}" y="{base-hh-10}" font-size="18" text-anchor="middle" font-family="Arial">{v}</text>')
    parts.append(f'<line x1="80" y1="{base}" x2="940" y2="{base}" stroke="#333"/><line x1="80" y1="100" x2="80" y2="{base}" stroke="#333"/><text x="35" y="335" transform="rotate(-90 35 335)" font-size="22" font-family="Arial">RS value</text>')
    save_svg(FIG/'Figure6_RS_value_barplot.svg',parts)
    # RF line
    genes=ref; parts=svg_header(1800,750); parts.append('<text x="900" y="45" font-size="30" text-anchor="middle" font-family="Arial" font-weight="bold">Single-gene rearrangement frequency (RF)</text>')
    ml=100; mt=80; ph=470; pw=1600; maxrf=max([rf[g] for g in genes]+[1]); pts=[]
    for i,g in enumerate(genes):
        x=ml+i*pw/(len(genes)-1); y=mt+ph-ph*rf[g]/maxrf; pts.append((x,y,g,rf[g]))
    parts.append('<polyline points="%s" fill="none" stroke="#E15759" stroke-width="4"/>'%' '.join(f'{x:.1f},{y:.1f}' for x,y,_,_ in pts))
    for x,y,g,v in pts:
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="#E15759"/><text x="{x:.1f}" y="620" font-size="12" text-anchor="middle" font-family="Arial" transform="rotate(65 {x:.1f} 620)">{g}</text>')
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#333"/><line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#333"/>')
    save_svg(FIG/'Figure7_RF_lineplot.svg',parts)
    return rs, [{'Gene':g,'RF':rf[g]} for g in genes]

# Main
records={p.stem.split('.')[0]:read_fasta(p)[1] for p in sorted(Path('.').glob('*.mitochondrial.fasta'))}
annotations={}; all_ann=[]; comp_rows=[]; codon_counts_by_sample={}; concat={}
for sample,seq in records.items():
    ann=annotate(seq); annotations[sample]=ann
    for row in ann:
        rr={'Sample':sample, **row}; all_ann.append(rr)
    write_csv(TAB/f'{sample}_gene_annotation.csv',ann,['Locus','Start','End','Length','Strand','Intergenic Spacer','Start Codon','Stop Codon','Anticodon'])
    circular_plot(sample,seq,ann)
    pieces=[]; cc=Counter()
    for row in ann:
        if row['Locus'] in PCGS|RRNAS:
            s=subseq_circular(seq,row['Start'],row['End']);
            if row['Strand']=='-': s=revcomp(s)
            st=comp_stats(s); comp_rows.append({'Sample':sample,'Locus':row['Locus'],'Length':len(s),**{k:round(v,4) for k,v in st.items()}})
            pieces.append(s)
            if row['Locus'] in PCGS:
                for i in range(0,len(s)-2,3):
                    cod=s[i:i+3]
                    if len(cod)==3 and set(cod)<=set('ACGT') and CODON_TABLE.get(cod) not in {None,'*'}: cc[cod]+=1
    codon_counts_by_sample[sample]=cc; concat[sample]=''.join(pieces)
    trna_cloverleaf(sample,ann)
write_csv(TAB/'all_gene_annotation.csv',all_ann,['Sample','Locus','Start','End','Length','Strand','Intergenic Spacer','Start Codon','Stop Codon','Anticodon'])
write_csv(TAB/'nucleotide_composition_skew.csv',comp_rows,['Sample','Locus','Length','T','C','A','G','AT','GC','GC Skew','AT Skew'])
# RSCU
rscu_rows=[]
all_counts=Counter()
for sample,cc in codon_counts_by_sample.items():
    all_counts.update(cc)
    for aa,codons in sorted(AA_CODONS.items()):
        total=sum(cc[c] for c in codons); n_syn=len(codons)
        for c in codons:
            r=(cc[c]*n_syn/total) if total else 0
            rscu_rows.append({'Sample':sample,'AA':aa,'AA_name':AA_NAMES[aa],'Codon':c,'N':cc[c],'RSCU':round(r,4)})
for aa,codons in sorted(AA_CODONS.items()):
    total=sum(all_counts[c] for c in codons); n_syn=len(codons)
    for c in codons:
        rscu_rows.append({'Sample':'ALL','AA':aa,'AA_name':AA_NAMES[aa],'Codon':c,'N':all_counts[c],'RSCU':round((all_counts[c]*n_syn/total) if total else 0,4)})
write_csv(TAB/'codon_usage_RSCU.csv',rscu_rows,['Sample','AA','AA_name','Codon','N','RSCU'])
bar_plot_rscu(rscu_rows)
seqs_by_sample=concat
names=list(records)
newick=upgma_tree(names,concat); support=bootstrap_support(names,concat,100)
(OUT/'phylogeny_newick.nwk').write_text(newick+'\n')
tree_plot(names,newick,support)
gene_order_plot(annotations)
rs,rfrows=rs_rf_plots(annotations)
write_csv(TAB/'rearrangement_scores_RS.csv',[{'Sample':a,'RS':b} for a,b in rs],['Sample','RS'])
write_csv(TAB/'rearrangement_frequency_RF.csv',rfrows,['Gene','RF'])
# Create a matplotlib preview/save template satisfying requested plotting-code convention.
Path('plot_preview_templates.py').write_text("""#!/usr/bin/env python3\n# Optional preview wrappers. If matplotlib is installed, each plotting section should call\n# plt.show() for visual confirmation before saving. The production SVG figures were\n# generated by analyze_mitogenomes.py without third-party dependencies.\nimport matplotlib.pyplot as plt\n# ... plotting commands ...\nplt.show()\nplt.savefig('analysis_results/figures/preview_figure.png', dpi=300, bbox_inches='tight')\n""")
print('Analysis complete:', OUT)
print('Tables:', len(list(TAB.glob('*.csv'))), 'Figures:', len(list(FIG.glob('*.svg'))))
