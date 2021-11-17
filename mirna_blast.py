#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 10:40:37 2021

@author: labo

Filter false positives (rRNA, tRNA, snRNA...) and overlaping miRNAs
Intersect coordinates with Ensembl annotation
Run blast against other species for conservation analysis
Classify and name miRNAs
Filter out duplicated miRNAs for quantification (cluster by sequence similarity)
Find miRNA clusters in genome
"""
import os
import re
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from upsetplot import UpSet, from_contents
from pybedtools import BedTool
os.chdir("/home/labo/datuak/mirna_atlas/")

##############################################################################
## Mature sequence level
# Make list of mircarta mirna fasta files
files = []
for filename in os.listdir():
    if filename.endswith("mirnas.fa"):
        files.append(filename)

# Make dictionary with mircarta and mirbase ids
ids_dicts = {}
for filename in os.listdir():
    if filename.endswith("miRCarta.csv"):
        mircartdict = {}
        mircart = pd.read_csv(filename, sep=",")
        mircart = mircart[mircart["miRCarta"].isna() == False]
        for entry in mircart.values:
            ids = re.split("; ", entry[2])
            for entry2 in ids:
                mircartdict[entry2] = entry[0]
        
        ids_dicts[filename[0:3]] = mircartdict

# Cluster mirna sequences
subprocess.run("cd-hit-est -i mature3.fa -o mature100.fa -c 1 -n 10 -d 0 -M 8000 -T 4", shell=True)

# Run blast as subprocess
for file in files:
    outid = file[:-3]
    out = file[:-2] + "out"
    subprocess.run("makeblastdb -in " + file + " -dbtype nucl -parse_seqids -out " + outid, shell=True)
    subprocess.run('blastn -task blastn-short -word_size 4 -num_alignments 100 -reward 5 -penalty -4 -gapopen 25 -gapextend 10 -query mature3.fa -db ' + outid + ' -out ' + out + ' -outfmt "6 std qlen slen qcovs"', shell=True)

# Filter results
outfiles = []
for filename in os.listdir():
    if filename.endswith(".out"):
        outfiles.append(filename)

results = []
for file in outfiles:
    blast_head = ["qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen", "qcovs"]
    blast_out = pd.read_csv(file, sep="\t", names=blast_head)
    blast_out = blast_out[blast_out["evalue"] <= 0.001]
    blast_out = blast_out[blast_out["pident"] > 95]
    
    blast_out_f = pd.DataFrame(columns=blast_head)
    for mirna in blast_out.qaccver.unique():
        mir_blast = blast_out[blast_out["qaccver"] == mirna]
        largest = mir_blast.nlargest(1,"pident")
        blast_out_f = blast_out_f.append(mir_blast[mir_blast["pident"] == largest.iloc[0,2]], ignore_index=True)
    sp = file[0:3]
    spmb = sp+"_mb"
    blast_out_f = blast_out_f.set_index("qaccver")
    blast_out_f[sp] = blast_out_f["saccver"].map(ids_dicts[sp])
    blast_out_f[spmb] = blast_out_f["saccver"]
    results.append(blast_out_f[[spmb, sp]])

# Compare species
res = pd.DataFrame().join(results, how="outer")

##############################################################################
## Precursor sequence level
# separate mirbase fasta files into species
species = ["chi", "bta", "eca", "ssc", "hsa"]    
for sp in species:
    seqs = []
    filename = sp + "_precursors.fa"
    for seq_record in SeqIO.parse("pre/hairpin.fa", "fasta"):
        rnaid = seq_record.id
        if rnaid.startswith(sp):
            seqs.append(seq_record)
    SeqIO.write(seqs, filename, "fasta")


# Identify mirnas classified as mirbase mirnas but with different location
subprocess.run("gff2bed < bed/oar.gff3 > bed/oar.bed", shell=True)

mirbasenames = ["chr", "start", "end", "id", "score", "strand",
                "score2", "type", "score3", "mirid", "alias", "name", "derive"]
mirbasebed = pd.read_csv("bed/oar.bed", sep="\t|;", names=mirbasenames, engine="python", index_col=False)
mirbasebed = mirbasebed.replace("Name=", "", regex=True)
mirbasebed = mirbasebed.replace("chr", "", regex=True)
mirbasebed = mirbasebed[mirbasebed["type"]=="miRNA_primary_transcript"]
mirbasebed[mirbasenames[:-1]].to_csv("bed/oar2.bed", sep="\t", index=False, header=None)

a = BedTool("bed/oar2.bed")
b = BedTool("bed/gure_premirnak_coord_known_oar4.0.bed")
c = BedTool.intersect(b, a, s=True, wo=True, f=0.001).moveto("bed/mirbase_intersect.bed")

mirbasebed = pd.read_csv("bed/mirbase_intersect.bed", sep="\t", index_col=False,
                         names=["chr2", "start2", "end2", "id2", "score4", "strand2"]+mirbasenames[:12]+["overlap"])
mirbasebed = mirbasebed.sort_values(by='overlap', ascending=False)
mirbasegood = mirbasebed.drop_duplicates(subset='mirid', keep="first")
len(mirbasegood["mirid"].unique())
mir_bad = [] # overlapping to mirbase: removed
for i in mirbasebed["id2"].values: 
    if i not in mirbasegood["id2"].values:
        mir_bad.append(i)
res_known = pd.read_csv("bed/result_26_07_2021_t_13_36_10_mirbase.csv", sep="\t", index_col=False)
mir_new = [] # mirbase sequence but other location: kept
seqs = [] # novel sequences (mirbase copies + completely novel)
seqs_good = [] # mirbase mature sequences for clustering
for i in res_known.values: 
    if i[0] not in mirbasegood["id2"].values:
        if i[0] not in mir_bad:
            mir_new.append(i[0])
            seqid = i[0] + "." + i[9]
            record = SeqRecord(Seq(i[15]), id=seqid,name=seqid,description=seqid)
            seqs.append(record)
    elif i[0] in mirbasegood["id2"].values:
        seqid = i[9]
        record = SeqRecord(Seq(i[13]), id=seqid,name=seqid,description=seqid)
        seqs_good.append(record)
        
# Edit BED file
bednames = ["chr", "start", "end", "id", "score", "strand"]
bed = pd.read_csv("bed/gure_premirnak_coord.bed", sep="\t", names=bednames)
mirbase = bed[bed["id"].str.match('oar')]
mirbase["temp"] = mirbase.apply(lambda x: '..'.join([str(x['start']),str(x['end'])]),axis=1)
mirbase["precursor coordinate"] = mirbase.apply(lambda x: ':'.join([x['chr'],x['temp'], x["strand"]]),axis=1)
mirbase = mirbase.merge(res_known, on="precursor coordinate")
bednew = bed[bed["id"].str.match('oar')==False]
mirbase2 = []
for i in mirbase.values:
    if i[8] in mir_new:
        seqid = i[8] + "." + i[3]
        i2 = i[:6]
        i2[3] = seqid
        mirbase2.append(i2)
    if i[8] not in mir_new and i[8] not in mir_bad:
        mirbase2.append(i[:6])
mirbase2 = pd.DataFrame(mirbase2, columns=bednames)
bed_edited = bednew.append(mirbase2)
bed_edited.to_csv("bed/gure_premirnak_coord_edited.bed", sep="\t", index=False, header=None)

# Intersect all mirnas between them to remove overlaping novel miRNAs
# We select the longest of each overlaping pair

a = BedTool("bed/gure_premirnak_coord_edited.bed")
c = BedTool.intersect(a, a, s=True, f=0.001,r=True, wo=True).moveto("bed/premirnak_intersect.bed")
ensnames = ["chr", "start", "end", "id", "score", "strand",
            "chr2", "start2", "end2", "id2", "score2", "strand2", "overlap"]
prebed = pd.read_csv("bed/premirnak_intersect.bed", sep="\t", names=ensnames)
prebed = prebed[prebed["id"]!=prebed["id2"]]
prebed["len1"] = prebed["end"] - prebed["start"]
prebed["len2"] = prebed["end2"] - prebed["start2"]
prebed = prebed[prebed["len1"]>prebed["len2"]]
# remove a mirbase ovelaped sequence if it was previously added to seqs
for mir in prebed["id2"].values:
    for i, seq in enumerate(seqs):
        if seqs[i].id == mir:
            del seqs[i]

# Filter false positive mirnas with other rnas
os.chdir("pre/")
subprocess.run("makeblastdb -in vault_tRNA_snRNA_snoRNA_rRNA_SRP_Y.fasta -dbtype nucl -parse_seqids -out rnacentral_oar", shell=True)
subprocess.run('blastn -task blastn -num_alignments 20 -query premirna3.fa -db rnacentral_oar -out premirna3_filter.out -outfmt "6 std qlen slen qcovs"', shell=True)
blast_head = ["qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen", "qcovs"]
filter_out = pd.read_csv("premirna3_filter.out", sep="\t", names=blast_head)
filter_out = filter_out[filter_out["evalue"] <= 0.001]
len(filter_out["qaccver"].unique())

for seq_record in SeqIO.parse("premirna3.fa", "fasta"):
    rnaid = seq_record.id
    if rnaid not in filter_out["qaccver"].unique():
        if rnaid not in prebed["id2"].unique():
            seqs.append(seq_record)

# Save fasta file of filtered premirnas (only novel precursors)
SeqIO.write(seqs, "premirna3_filter.fa", "fasta")
# Save fasta file of mirbase mature sequences for clustering (not star seqs)
os.chdir("..")
SeqIO.write(seqs_good, "mirbase_mature.fa", "fasta")

# Make list of mirbase mirna fasta files
os.chdir("pre")
files = []
for filename in os.listdir():
    if filename.endswith("precursors.fa"):
        files.append(filename)

# Run blast as subprocess
for file in files:
    outid = file[:-3]
    out = file[:-2] + "out"
    subprocess.run("makeblastdb -in " + file + " -dbtype nucl -parse_seqids -out " + outid, shell=True)
    subprocess.run('blastn -task blastn -num_alignments 100 -query premirna3_filter.fa -db ' + outid + ' -out ' + out + ' -outfmt "6 std qlen slen qcovs"', shell=True)

# Filter results
outfiles = []
for filename in os.listdir():
    if filename.endswith("precursors.out"):
        outfiles.append(filename)

results = []
for file in outfiles:
    blast_head = ["qaccver","saccver","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen", "qcovs"]
    blast_out = pd.read_csv(file, sep="\t", names=blast_head)
    blast_out = blast_out[blast_out["evalue"] <= 0.01]
    #blast_out = blast_out[blast_out["pident"] > 85]
    #blast_out = blast_out[blast_out["mismatch"] < 6]
    blast_out = blast_out[blast_out["qcovs"] > 80]
    
    blast_out_f = pd.DataFrame(columns=blast_head)
    for mirna in blast_out.qaccver.unique():
        mir_blast = blast_out[blast_out["qaccver"] == mirna]
        largest = mir_blast.nlargest(1,"pident")
        blast_out_f = blast_out_f.append(mir_blast[mir_blast["pident"] == largest.iloc[0,2]], ignore_index=True)
    sp = file[0:3]
    spid = sp+"_%"
    blast_out_f = blast_out_f.set_index("qaccver")
    blast_out_f[sp] = blast_out_f["saccver"]
    blast_out_f[spid] = blast_out_f["pident"]
    results.append(blast_out_f[[sp,spid]])

# Compare species and select closer mirna
res = pd.DataFrame().join(results, how="outer")
len(res.index.unique())

# Select identical hits from closest species OR the one with best identity
besthit = []
bestid = []
for hit in res.values:
    if hit[3] == 100:
        besthit.append(hit[2])
        bestid.append(hit[3])
    elif hit[9] == 100:
        besthit.append(hit[8])
        bestid.append(hit[9])
    elif hit[1] == 100:
        besthit.append(hit[0])
        bestid.append(hit[1])
    elif hit[7] == 100:
        besthit.append(hit[6])
        bestid.append(hit[7])
    elif hit[5] == 100:
        if pd.isna(hit[4])==False:
            besthit.append(hit[4])
            bestid.append(hit[5])
        else:
            bestid.append("")
    else:
        ident = [hit[3],hit[9],hit[1],hit[7],hit[5]]
        ident = [0 if pd.isna(x) else x for x in ident]
        sps = [hit[2],hit[8],hit[0],hit[6],hit[4]]
        hdict = {sps[i]: ident[i] for i in range(len(sps))}
        if pd.isna(max(hdict,key=hdict.get)):
            besthit.append("")
            bestid.append("")
        else:
            besthit.append(max(hdict,key=hdict.get))
            bestid.append(max(ident))

res["best"] = besthit
res["bestid"] = bestid

homologs = []
for mir in res.index.unique():
    hits = res[res.index==mir]
    ids = []
    identities = []
    for hit in hits.values:
        ids.append(hit[10])
        identities.append(hit[11])
    ids2 = []
    for i in ids:
        if (i != "" ) and (i not in ids2):
            ids2.append(i)
    ids2 = "/".join(ids2)
    if ids2 != "":
        homologs.append([mir, ids2, max(identities)])

homologs = pd.DataFrame(homologs, columns=["novel", "homolog", "identity"])
len(homologs["homolog"].unique())
len(homologs["novel"].unique())

# Number of miRNAs in each species with significant hit
conserved = {"Goat":res[res["chi"].isna()==False].index.unique().tolist(),
             "Cattle":res[res["bta"].isna()==False].index.unique().tolist(),
             "Horse":res[res["eca"].isna()==False].index.unique().tolist(),
             "Pig":res[res["ssc"].isna()==False].index.unique().tolist(),
             "Human":res[res["hsa"].isna()==False].index.unique().tolist()}

upsetdata = from_contents(conserved)
upsetdata["type"] = "Other"
upsetdata["type"] = np.where(upsetdata.id.str.contains('2284|2285'), "mir-2284/2285", upsetdata.type)
# Plot data
upset = UpSet(upsetdata, subset_size='count',
              show_counts=True,
              show_percentages=False,
              intersection_plot_elements=0)
upset.add_stacked_bars(by="type", colors = "Set1",
                       title="miRNA count", elements=6)
upset.plot()

# Rename 2284/2285 family mirnas
d228 = {}
d228n = 0
for mir in homologs.values:
   if "2284" in mir[1] or "2285" in mir[1]:
       d228n += 1
       new = "bta-mir-2284-novel" + str(d228n)
       d228[mir[0]] = new
   else:
       d228[mir[0]] = mir[1]

homologs["homolog"] = homologs["novel"].map(d228)
len(homologs["homolog"].unique())
len(homologs["novel"].unique())

# Make fasta file with changed names
seqsf = []
othern = 0
allnovel = []
for seq_record in SeqIO.parse("premirna3_filter.fa", "fasta"):
    rnaid = seq_record.id
    newid = ""
    if rnaid not in d228:
        othern += 1
        newid = "novel-mir-" + str(othern)
    else:
        newid = d228[seq_record.id]
    seq_record.id = newid
    seqsf.append(seq_record)
    allnovel.append([rnaid, newid])
    
SeqIO.write(seqsf, "novel_pre.fa", "fasta")

# Intersect all miRNAs and Ensembl miRNAs
os.chdir("../")
os.chdir("bed")

a = BedTool("gure_premirnak_coord_edited.bed")
b = BedTool("ensenble_premirnak_coord.bed")
c = BedTool.intersect(a, b, s=True, wo=True, f=0.001).moveto("ensembl_intersect.bed")

# Make table with locations and intersected Ensembl sheep mirnas
bednames = ["chr", "start", "end", "id", "score", "strand"]
bed = pd.read_csv("gure_premirnak_coord_edited.bed", sep="\t", names=bednames)

allnoveldf = pd.DataFrame(allnovel, columns=["id", "name"])
bed2 = bed.merge(allnoveldf, on="id", how="inner")

ident_dict = homologs.set_index('novel').to_dict()['identity']
bed2["identity"] = bed2["id"].map(ident_dict)

mirbase = bed[bed["id"].str.startswith('oar')]
mirbase["name"] = mirbase["id"]
mirbase["identity"] = 0
len(mirbase["id"].unique())
bed2 = bed2.append(mirbase)

ensnames = ["chr", "start", "end", "id", "score", "strand",
            "chr2", "start2", "end2", "ens", "score2", "strand2", "overlap"]
ens = pd.read_csv("ensembl_intersect.bed", sep="\t", names=ensnames)
ens = ens[["chr", "start", "end", "ens"]]

bed2 = bed2.merge(ens, on=["chr", "start", "end"], how="left")
bed2 = bed2[bed2["score"]==0]

# Number of Ensembl miRNAs
bed2["ens"].nunique()
bed2[bed2["name"].str.startswith("oar")]["ens"].nunique()
bed2[bed2["name"].str.startswith("oar")==False]["ens"].nunique()

len(bed2[bed2["identity"]==100])
len(bed2[bed2["name"].str.startswith("oar")])

##############################################################################
## Filtering the expression matrix
# Remove repeated sequences (miRNAs with multiple loci) using mature sequences
os.chdir("..")
os.chdir("quant")

allmature = []
for seq_record in SeqIO.parse("Qmature.fa", "fasta"):
    allmature.append(seq_record)
for seq_record in SeqIO.parse("mirbase_mature.fa", "fasta"):
    allmature.append(seq_record)
SeqIO.write(allmature, "Qmature_all.fa", "fasta")

# Cluster mirna sequences
subprocess.run("cd-hit-est -i Qmature_all.fa -o Qmature_100.fa -c 1 -n 10 -d 0 -M 8000 -T 4", shell=True)

clusters = {}
with open("Qmature_100.fa.clstr") as f:
    for line in f:
        if re.search("Cluster", line):
            clid = line.rstrip().split(" ")[1]
            if clid not in clusters:
                clusters[clid] = []
        elif re.search("nt", line):
            mirid = re.split(">|\.", line.rstrip())[1]
            if len(clusters[list(clusters)[-1]]) == 0:
                clusters[list(clusters)[-1]] = mirid

# Add "q" to the representative miRNA of the cluster that is being quantified
inv_clusters = {v: "q" for k, v in clusters.items()}

bed2["temp"] = bed2.apply(lambda x: str(x['id']).split(".")[0] ,axis=1)
bed2["quant"] = ""
bed2["quant"] = bed2["temp"].map(inv_clusters)
del bed2["temp"]

# Edit duplicates with same name for quantification:
# Keep the loci with ensembl match OR add suffix if no ens loci

bed2["index"] = bed2["id"]
bed2 = bed2.set_index("index")
qbed = bed2[bed2["quant"]=="q"]
qbed[qbed["name"].str.startswith("oar")].nunique()
qbed[qbed["name"].str.startswith("oar")==False].nunique()
dupq = qbed[qbed["name"].isin(qbed["name"][qbed["name"].duplicated()])]
dupq = dupq.sort_values(by="name")

remov_dup = []
for mir in dupq["name"].unique():
    if mir.startswith("oar") == False:
        mirdf = dupq[dupq["name"]==mir]
        ids = []
        ensids = []
        newids = []
        strand = []
        names = []
        for i in mirdf.values:
            ensids.append(i[8])
            ids.append(i[3])
            strand.append(i[5])
            names.append(i[6])
        if type(ensids[0]) != str and type(ensids[1]) != str:
            if strand[0] != strand[1]:
                newids.append("".join(names[0]+strand[0]))
                newids.append("".join(names[1]+strand[1]))
            elif strand[0] == strand[1]:
                newids.append("".join(names[0]+"-1"))
                newids.append("".join(names[1]+"-2"))
            bed2['name'][ids[0]] = newids[0]
            bed2['name'][ids[1]] = newids[1]
        elif type(ensids[0]) == str:
            remov_dup.append(ids[1])
        elif type(ensids[1]) == str:
            remov_dup.append(ids[0])

        
bed2['quant'] = bed2.apply(lambda x: "" if x['id'] in remov_dup else x['quant'], axis=1)

qbed = bed2[bed2["quant"]=="q"]
qbed[qbed["name"].str.startswith("oar")].nunique()
qbed[qbed["name"].str.startswith("oar")==False].nunique()

dupq = qbed[qbed["name"].isin(qbed["name"][qbed["name"].duplicated()])]
dupq = dupq.sort_values(by="name")

# Save table
bed2.to_csv("all_mirnas.bed", sep="\t", index=False, header=None)

# Prepare supplementary table
os.chdir("..")
sup = bed2
sup["id"] = sup.apply(lambda x: str(x['id']).split(".")[0] ,axis=1)
sup["precursor coordinate"] = sup["start"].astype(str) + ".." + sup["end"].astype(str)
sup["precursor coordinate"] = sup["chr"] + ":" + sup["precursor coordinate"] + ":" + sup["strand"] 

sup.columns
mirdeep_known = pd.read_csv("bed/result_26_07_2021_t_13_36_10_mirbase.csv", sep="\t", index_col=False)
mirdeep_known["miRDeep2_probability"] = mirdeep_known["estimated probability that the miRNA is a true positive"]
mirdeep_known["id"] = mirdeep_known["mature miRBase miRNA"]
mirdeep_novel = pd.read_csv("bed/result_26_07_2021_t_13_36_10_novel.csv", sep="\t", index_col=False)
mirdeep_novel["id"] = mirdeep_novel["provisional id"]
mirdeep_novel["miRDeep2_probability"] = mirdeep_novel["estimated probability that the miRNA candidate is a true positive"]
keepcols = ['miRDeep2 score',
       'miRDeep2_probability',
       'consensus mature sequence',
       'consensus star sequence', 'consensus precursor sequence',
       'id', 'precursor coordinate']
mirdeep_known = mirdeep_known[keepcols]
mirdeep_novel = mirdeep_novel[keepcols]
mirdeep_known = sup.merge(mirdeep_known, on="precursor coordinate")
mirdeep_novel = sup.merge(mirdeep_novel, on="precursor coordinate")

mirdeep = mirdeep_known.append(mirdeep_novel)
del mirdeep["coord"]
del mirdeep["precursor coordinate"]
del mirdeep["id_x"]

mirdeep.to_csv("bed/all_mirnas_supplementary.csv", sep=",", index=False)


##############################################################################
## Get unique precursor sequences from mir2284/2285 family

bed3names = ["chr", "start", "end", "id", "score", "strand", "name", "identity", "ens", "quant"]
bed3 = pd.read_csv("quant/all_mirnas.bed", sep="\t", index_col=False, names=bed3names)
bed228 = bed3[bed3["name"].str.contains("2284|2285")]
bed228["provisional id"] = bed228.apply(lambda x: str(x['id']).split(".")[0] ,axis=1)
bed228.to_csv("2284/2284_mirnas.bed", sep="\t", index=False)
bed3[["chr", "start", "end", "id", "score", "strand", "name"]].to_csv("quant/all_mirnas_simple.bed", sep="\t", index=False, header=None)

seqs228 = []
for seq_record in SeqIO.parse("pre/novel_pre.fa", "fasta"):
    if seq_record.id in bed228["name"].values:
        seqs228.append(seq_record)

SeqIO.write(seqs228, "2284/mir2884_pre.fa", "fasta")

# mature sequences
res_novel = pd.read_csv("bed/result_26_07_2021_t_13_36_10_novel.csv", sep="\t", index_col=False)
res_novel["id"] = res_novel.apply(lambda x: ".".join([str(x["provisional id"]),  str(x["example miRBase miRNA with the same seed"])]) ,axis=1)

bed228 = bed228.merge(res_novel, on="id")
bed228[bed228["quant"]=="q"].nunique()
bed228.nunique()

seqs228_m = []
seqs228_s = []
for mir in bed228.values:
    seqidm = mir[6]
    recordm = SeqRecord(Seq(mir[24]), id=seqidm,name=seqidm,description=seqidm)
    seqs228_m.append(recordm)
    seqids = mir[6]
    records = SeqRecord(Seq(mir[25]), id=seqids,name=seqids,description=seqids)
    seqs228_s.append(records)
SeqIO.write(seqs228_m, "2284/mir2884_mature.fa", "fasta")
SeqIO.write(seqs228_s, "2284/mir2884_star.fa", "fasta")

# Cluster identical sequences from mir2284/2285 family
subprocess.run("cd-hit-est -i 2284/mir2884_pre.fa -o 2284/mir2884_pre_100.fa -c 1 -n 10 -d 0 -M 8000 -T 4", shell=True)
subprocess.run("cd-hit-est -i 2284/mir2884_mature.fa -o 2284/mir2884_mature_100.fa -c 1 -n 10 -d 0 -M 8000 -T 4", shell=True)
subprocess.run("cd-hit-est -i 2284/mir2884_star.fa -o 2284/mir2884_star_100.fa -c 1 -n 10 -d 0 -M 8000 -T 4", shell=True)

##############################################################################
# plot locations with for chromoMap R package
bed3 = bed3[["name", "chr", "start", "end"]]
bed3["type"] = "Conserved"
bed3['type'] = bed3.apply(lambda x: "Novel" if "novel" in x['name'] else x['type'], axis=1)
bed3['type'] = bed3.apply(lambda x: "mir-2284/5" if x['name'] in bed228["name"].values else x['type'], axis=1)
bed3['type'] = bed3.apply(lambda x: "miRBase" if x['name'].startswith("oar") else x['type'], axis=1)
bed3.to_csv("bed/all_mirnas_chromomap.bed", sep="\t", index=False, header=None)

chrlen = pd.read_csv("/home/labo/datuak/atlas/cage/chromsizes_ncbi2ens.txt", names=["ncbi", "len", "ens"], sep=" ")
chrlen["start"] = 1
chrlen = chrlen[["ens", "start", "len"]]
chrlen = chrlen.sort_values(by="ens")
chrlenx = chrlen[27:28]
chrlen = chrlen[chrlen["ens"]!="MT"]
chrlen = chrlen[chrlen["ens"]!="X"]
chrlen["ens"] = chrlen["ens"].astype(str).astype(int)
chrlen = chrlen.sort_values(by="ens")
chrlen = chrlen.append(chrlenx)
chrlen.to_csv("bed/chromosomes_oar_chromomap.bed", sep="\t", index=False, header=None)

##############################################################################
# Find clusters of miRNAs in genome
a = BedTool("quant/all_mirnas_simple.bed")
b = BedTool.sort(a)
c = BedTool.cluster(b, d=10000).moveto("quant/all_mirnas_clusters.bed")

geneclusters = {}
nclust = 0
with open("quant/all_mirnas_clusters.bed") as f:
    for line in f:
        tline = line.rstrip().split("\t")
        clid = tline[7]
        if clid not in geneclusters:
            geneclusters[clid] = []
            geneclusters[clid].append(tline[6])
        else:
            geneclusters[clid].append(tline[6])
for key in geneclusters.keys():
    if len(geneclusters[key]) > 1:
        print(len(geneclusters[key]))
        print(geneclusters[key])
        nclust+=1