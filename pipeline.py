import os
import glob
from functools import partial
import time
from collections import Counter
import itertools
import multiprocessing
import gzip
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, precision_score, recall_score, accuracy_score
from sklearn.model_selection import StratifiedKFold
from concurrent.futures import ProcessPoolExecutor

##############################
# Configuration
##############################
FASTQ_DIR     = '/home/b3nkaza/dataset-2'            # directory containing all .fastq files
METADATA_CSV  = 'multivirus_labels.csv'     # must have columns Run, Library Name, Virus, HPI

KMER_RANGE    = range(1, 21)                 # e.g. unit‐test with range(3,5) before full run
N_BOOT        = 50                          # bootstrap replicates for clustering
N_SPLITS      = 4                           # k‑fold CV folds
N_WORKERS     = 64                           # parallel workers for k‑mer counting
RF_KWARGS     = {'n_estimators':200,        # RF hyperparameters
                 'n_jobs':1,
                 'random_state':0}

##############################
# 1) ACGT‑only normalized k‑mer counting
##############################
def count_kmers(fq, k):
    proc = multiprocessing.current_process().name
    name = os.path.basename(fq)
    print(f"[{proc}] Counting {k}-mers in {name}")
    counts = Counter()
    valid = set('ACGT')
    open_fn = gzip.open if fq.endswith(".gz") else open
    with open_fn(fq, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seq = line.strip().upper()
                for j in range(len(seq)-k+1):
                    kmer = seq[j:j+k]
                    if all(c in valid for c in kmer):
                        counts[kmer] += 1
    total = sum(counts.values())
    if total>0:
        for kmer in counts:
            counts[kmer] /= total
    return counts

##############################
# 2) Merge into DataFrame
##############################
def merge_counts(dicts, names):
    all_kmers = sorted(set(itertools.chain.from_iterable(d.keys() for d in dicts)))
    df = pd.DataFrame(0.0, index=names, columns=all_kmers)
    for n, d in zip(names, dicts):
        for kmer, freq in d.items():
            df.at[n, kmer] = freq
    return df

##############################
# 3) Cluster purities
##############################
def cluster_rep_purities(X, y, n_boot=N_BOOT, method='average'):
    out = []
    classes = np.unique(y)
    for rep in range(n_boot):
        idx = np.random.choice(len(X), len(X), replace=True)
        Xb, yb = X[idx], y[idx]
        Z = linkage(Xb, method=method)
        cl = fcluster(Z, t=len(classes), criterion='maxclust')
        df = pd.DataFrame({'true':yb,'cl':cl})
        purity = (df.groupby('cl')['true']
                    .agg(lambda x: x.value_counts().iloc[0]).sum())/len(df)
        rec = {'replicate':rep,'purity':purity}
        for c in classes:
            sub = df[df.true==c]
            if len(sub):
                maj = sub.cl.value_counts().idxmax()
                rec[f'pc_{c}'] = (sub.cl==maj).mean()
            else:
                rec[f'pc_{c}'] = np.nan
        out.append(rec)
    return out

##############################
# 4) RF metrics
##############################
def rf_fold_metrics(X, y, n_splits=N_SPLITS, rf_kwargs=None):
    out = []
    classes = np.unique(y)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
    for fold, (tr, te) in enumerate(skf.split(X, y), 1):
        clf = RandomForestClassifier(**(rf_kwargs or {}))
        clf.fit(X[tr], y[tr])
        pred = clf.predict(X[te])
        m = {
            'fold': fold,
            'f1':        f1_score(y[te], pred, average='weighted'),
            'precision': precision_score(y[te], pred, average='weighted'),
            'recall':    recall_score(y[te], pred, average='weighted'),
            'accuracy':  accuracy_score(y[te], pred)
        }
        for c in classes:
            true_c = (y[te]==c)
            pred_c = (pred==c)
            m[f'pc_{c}_f1']        = f1_score(true_c, pred_c, zero_division=0)
            m[f'pc_{c}_precision'] = precision_score(true_c, pred_c, zero_division=0)
            m[f'pc_{c}_recall']    = recall_score(true_c, pred_c, zero_division=0)
        out.append(m)
    return out

##############################
# 5) Single‐run wrapper
##############################
def run_pipeline(label_array):
    """Runs the full pipeline on labels=label_array (same order as fastq_files)."""
    cluster_recs = []
    rf_recs      = []
    summary_recs = []

    for k in KMER_RANGE:
        print(f"\n>>> k={k} <<<")
        # k‐mer counting
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=N_WORKERS) as exe:
            counts = list(exe.map(partial(count_kmers, k=k), fastq_files))
        t_kmer = time.time()-t0

        # merge
        names = [os.path.basename(f).split('_')[0] for f in fastq_files]
        X_df = merge_counts(counts, names)
        nfeat = X_df.shape[1]
        X = X_df.values
        y = label_array  # already aligned by index
        
        # clustering
        t1 = time.time()
        creps = cluster_rep_purities(X, y, n_boot=N_BOOT)
        t_h = time.time()-t1
        for r in creps:
            r.update({'k':k})
        cluster_recs += creps

        # RF
        t2 = time.time()
        rfolds = rf_fold_metrics(X, y, n_splits=N_SPLITS, rf_kwargs=RF_KWARGS)
        t_r = time.time()-t2
        for r in rfolds:
            r.update({'k':k})
        rf_recs += rfolds

        # summary
        summary_recs.append({
            'k':k, 'n_features':nfeat,
            'time_kmer':t_kmer,
            'time_hclust':t_h,
            'time_rf':t_r
        })

    return cluster_recs, rf_recs, summary_recs

##############################
# 6) Main
##############################
if __name__ == '__main__':
    import time

    # load + filter metadata
    md = pd.read_csv(METADATA_CSV, usecols=['Run','Virus','HPI'])
    md = md[pd.to_numeric(md['HPI'],errors='coerce').isin([12,24])].set_index('Run')
    print("Runs:", md.index.tolist())

    # find FASTQs
    all_fqs = sum((glob.glob(os.path.join(FASTQ_DIR,p)) for p in ['*.fastq','*.fastq.gz']),[])
    fastq_files = [fq for fq in sorted(all_fqs)
                   if os.path.basename(fq).split('.')[0] in md.index]
    if not fastq_files:
        raise SystemExit("No FASTQs after HPI filter")
    print(f"Using {len(fastq_files)} FASTQs")

    # True‐label run
    labels_true = md.loc[[os.path.basename(f).split('.')[0] for f in fastq_files],'Virus'].values
    (cl_true, rf_true, sum_true) = run_pipeline(labels_true)
    pd.DataFrame(cl_true).to_csv('cluster_replicates.csv', index=False)
    pd.DataFrame(rf_true).to_csv('rf_folds.csv',      index=False)
    pd.DataFrame(sum_true).to_csv('k_timing_summary.csv', index=False)
    print("→ Wrote original CSVs")

    # Permuted‐label run
    np.random.seed(0)
    labels_perm = np.random.permutation(labels_true)
    (cl_perm, rf_perm, sum_perm) = run_pipeline(labels_perm)
    pd.DataFrame(cl_perm).to_csv('cluster_replicates_perm.csv', index=False)
    pd.DataFrame(rf_perm).to_csv('rf_folds_perm.csv',      index=False)
    pd.DataFrame(sum_perm).to_csv('k_timing_summary_perm.csv', index=False)
    print("→ Wrote permuted CSVs")