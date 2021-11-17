#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 12:37:51 2021

@author: labo
"""
import pandas as pd
import seaborn as sns
import math
import os
import scipy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rcParams
import matplotlib
rcParams.update({'figure.autolayout': True})
os.chdir("/")

### Mean expression by miRNA type (novel, conserved, 2284 family)
# Open tissue info
coldata = pd.read_csv("expr/sample_info.csv", sep=",", index_col=0)
# Open TPM expression
expr = pd.read_table("expr/normalizedEM.csv", delimiter=",", index_col=0)
# Open miRNA classification
types = pd.read_table("bed/all_mirnas_chromomap.bed", delimiter="\t", index_col=False, names=["id", "chr", "start", "end", "type"])
typedict = types.set_index('id').to_dict()['type']

# miRNA custom colormap
cmap1=plt.cm.get_cmap('tab20b', 20)
cmap2=plt.cm.get_cmap('tab20c', 20)

mycols = []
for i in range(cmap1.N):
   rgba = cmap1(i)
   mycols.append(matplotlib.colors.rgb2hex(rgba))
for i in range(cmap2.N):
   rgba = cmap2(i)
   mycols.append(matplotlib.colors.rgb2hex(rgba))

mycmap = mcolors.ListedColormap(mycols[0:30])

# Tissue colormap
c25 = ["#1C86EE", "#E31A1C", "#008B00", "#6A3D9A", "#FF7F00", "#000000",
       "#FFD700", "#7EC0EE", "#FB9A99", "#90EE90", "#CAB2D6", "#FDBF6F",
       "#B3B3B3", "#EEE685", "#B03060", "#FF83FA", "#FF1493", "#0000FF",
       "#36648B", "#00CED1", "#00FF00", "#8B8B00", "#CDCD00", "#8B4500",
       "#A52A2A"]

mycols_t = c25[0:21]
mycmap_t = mcolors.ListedColormap(mycols)

# calculate means
def calc_means(counts, typedict):
    means_d = []
    for index, row in counts.iterrows():
        geneid = index.rsplit('-', 1)[0]
        means = row.mean()
        logmean = math.log2(means)
        if geneid in types["id"].values:
            mirtype = typedict[geneid]
            means_d.append([index, means, logmean, mirtype])
        elif index in types["id"].values:
            mirtype = typedict[index]
            means_d.append([index, means, logmean, mirtype])

    means_df = pd.DataFrame(means_d)
    means_df = means_df.set_index(0)
    means_df.columns = ["CPM", "log2 CPM", "miRNA type"]
    return means_df

means = calc_means(expr, typedict)

# Separate means by tissue
means_tissue_list = []
for tis in coldata["Tissue"].unique():
    expr_t = expr[coldata[coldata["Tissue"]==tis].index]
    expr_t = expr_t.loc[(expr_t!=0).any(axis=1)]
    means_t = calc_means(expr_t, typedict)
    means_t["tissue"] = tis
    means_tissue_list.append(means_t)

means_tissue = pd.concat(means_tissue_list)
means_2284 = means_tissue[means_tissue["miRNA type"]=="mir-2284/5"]

# Plot violin plot by miRNA type expression
colors = sns.color_palette("Set2", 4)
fig = plt.figure(figsize=(5,4), dpi=300)
bplot = sns.violinplot(y='log2 CPM', x='miRNA type', 
                 data=means, 
                 width=1,
                 palette=colors,
                 saturation=0.75)
plt.savefig("expr/average_means.svg")

# Mean values
means_type = means.groupby(['miRNA type'])['CPM'].mean()
means_error = means.groupby(['miRNA type'])['CPM'].sem()

# Plot tissues mir2284/5
fig = plt.figure(figsize=(9,4), dpi=300)
bplot = sns.violinplot(y='log2 CPM', x='tissue', 
                 data=means_2284, 
                 width=1,
                 palette=mycols_t,
                 saturation=0.75)
plt.xticks(rotation=90)
plt.title("mir-2284/2285 family expression")
plt.savefig("expr/average_means_2284_tissues.svg")

means_2284.groupby(['tissue'])['CPM'].median()
means_2284.groupby(['tissue'])['CPM'].sem()

# Perform mannwhitney U test
for tis in coldata["Tissue"].unique():
    st = scipy.stats.mannwhitneyu(means_2284["CPM"],
                                  means_2284[means_2284["tissue"]==tis]["CPM"],
                                  alternative="less")
    print(tis, st[1])


# Number of 2284 expressed in each tissue
meanexpr = pd.read_table("expr/normalized_expresion_per_tissue.csv", delimiter=",", index_col=0)

means_2284_1cpm = means_2284[means_2284["CPM"]>1]
means_2284["tissue"].value_counts(sort=False)
means_2284_1cpm = means_2284_1cpm["tissue"].value_counts(sort=False)
means_2284_1cpm = means_2284_1cpm.to_frame()
means_2284_1cpm["id"] = means_2284_1cpm.index
means_2284_1cpm = means_2284_1cpm.reindex(meanexpr.columns)

fig = plt.figure(figsize=(9,4), dpi=300)
bplot = sns.barplot(x= "id", y='tissue', 
                 data=means_2284_1cpm,
                 palette=mycols_t,
                 saturation=0.75)
plt.xticks(rotation=90)
plt.title("Number of mir-2284/2285 family miRNAs > 1 CPM")
plt.savefig("expr/number_2284_tissues.svg")

##############################################################################
### Most expessed miRNAs by tissue

topmirs = pd.DataFrame()
for tis in coldata["Tissue"].unique():
    sortedm = meanexpr[tis].sort_values(ascending=False)
    top10 = sortedm.nlargest(5)
    other = sortedm[sortedm.index.isin(top10.index)==False]
    otherdf = pd.DataFrame({top10.to_frame().columns[0]: [sum(other)]}, index=["other"])
    top11 = top10.to_frame().append(otherdf)
    topmirs = topmirs.join(top11, how="outer")

topmirs = topmirs.transpose()
topmirs = topmirs / 1000000
topmirs["Tissue"] = topmirs.index

# Plot
topmirs.plot(
        x = 'Tissue',
        kind = 'bar',
        stacked = True,
        width=0.9,
        rot=30,
        ylabel="miRNA expresion fraction",
        colormap=mycmap,
        figsize=(14,7)).legend(loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig("expr/top_mirnas.svg")

##############################################################################
# Plot most specific miRNA expression across tissues
expr2 = meanexpr.transpose()
spec = pd.read_table("expr/Tissuespecifity.csv", delimiter=",", names=["id","tsi","tissue"])
spec = spec.sort_values(by="tissue")
spec["type"] = "Conserved"
spec['type'] = spec.apply(lambda x: "Novel" if "novel" in x['id'] else x['type'], axis=1)
spec['type'] = spec.apply(lambda x: "mir-2284/5" if "2284-novel" in x['id'] else x['type'], axis=1)
spec['type'] = spec.apply(lambda x: "mir-2284/5" if "2285-novel" in x['id'] else x['type'], axis=1)
spec['type'] = spec.apply(lambda x: "miRBase" if x['id'].startswith("oar") else x['type'], axis=1)
spec.groupby("type").count()
spec.groupby("tissue").count()

# example single mir
expr2["oar-miR-370-5p"].plot(
        x = 'Tissue',
        kind = 'bar',
        width=0.9,
        rot=30,
        ylabel="miRNA expresion (CPM)",
        colormap=mycmap,
        figsize=(14,7))

# mirbase specific
spec_mirbase = spec[spec["type"]=="miRBase"]
expr2_mirbase = expr2[spec_mirbase["id"]]

expr2_mirbase.plot(
        subplots=True,
        grid=False,
        layout=(6,3),
        kind = 'bar',
        width=0.9,
        rot=90,
        ylabel="CPM",
        #colormap=mycmap,
        color=[mycols_t],
        figsize=(14,10),
        sharex=True, sharey=False, legend=False)
plt.savefig("expr/specific_mirbase.svg")

# top 1 per tissue
grouped_spec = spec.groupby("tissue")
maxspec = grouped_spec.max()
maxspec = maxspec.reset_index()
expr2_max = expr2[maxspec["id"]]

expr2_max.plot(
        subplots=True,
        grid=False,
        layout=(6,3),
        kind = 'bar',
        width=0.9,
        rot=90,
        ylabel="CPM",
        #colormap=mycmap,
        color=[mycols_t],
        figsize=(14,10),
        sharex=True, sharey=False, legend=False)
plt.savefig("expr/specific_maxpertissue.svg")

# selected mirnas for article figure
#spec_sorted = spec.sort_values(by="tissue")
selected = expr2[["bta-mir-208a-3p","eca-mir-8908b-1-3p","chi-mir-122-3p","hsa-mir-373--3p","oar-miR-433-5p","oar-miR-654-3p","chi-mir-133b-5p","bta-mir-6715-3p"]]
selected.plot(
        subplots=True,
        grid=False,
        layout=(8,1),
        kind = 'bar',
        width=0.9,
        rot=90,
        ylabel="CPM",
        #colormap=mycmap,
        color=[mycols_t],
        figsize=(5,12),
        sharex=True, sharey=False, legend=False)
plt.savefig("expr/specific_selected.svg")

# 2284 specific
spec_2284 = spec[spec["type"]=="mir-2284/5"]
expr2_2284 = expr2[spec_2284["id"]]

expr2_2284.plot(
        subplots=True,
        grid=False,
        layout=(4,4),
        kind = 'bar',
        width=0.9,
        rot=90,
        ylabel="CPM",
        #colormap=mycmap,
        color=[mycols_t],
        figsize=(14,7),
        sharex=True, sharey=False, legend=False)
plt.savefig("expr/specific_2284.svg")

##############################################################################
# Matrix for heatmap
expr_2284 = expr[expr.index.isin(means_2284.index)]
expr_2284.to_csv("2284/normalized_2284.csv")

##############################################################################
# T Specific miRNAs in testis and brain
expr_spec = meanexpr[meanexpr.index.isin(spec["id"])]
expr_spec = expr_spec.transpose()

brain_testis_all = []
brain_testis_filter = []
for mir in expr_spec.columns:
    max2 = expr_spec.nlargest(2, mir)
    if "Brain" in max2.index:
        if "Testis" in max2.index or "Epididymis" in max2.index:
            ratio = max2[mir][0] / max2[mir][1]
            brain_testis_all.append(mir)
            if 0.20 < ratio < 5: 
                brain_testis_filter.append(mir)

expr2_br_tes = expr2[brain_testis_filter]

expr2_br_tes.plot(
        subplots=True,
        grid=False,
        layout=(4,4),
        kind = 'bar',
        width=0.9,
        rot=90,
        ylabel="CPM",
        #colormap=mycmap,
        color=[mycols_t],
        figsize=(14,7),
        sharex=True, sharey=False, legend=False)
plt.savefig("expr/specific_brain_testis.svg")
