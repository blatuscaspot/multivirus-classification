import os
import glob
import time
import gzip
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from scipy.cluster.hierarchy import linkage, to_tree
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, precision_score, recall_score, accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# PARAMETERS
ALL_FASTQ_DIR = '/home/b3nkaza/dataset-2/datasets-unused/'  # Directory containing all processed .fastq files
METADATA_CSV  = '/home/b3nkaza/dataset-2/multivirus_labels.csv'
KMER_RANGE    = range(1, 16)
N_BOOT        = 50
BOOTSTRAP_FRACTION = 0.8
N_SPLITS      = 4
RF_KWARGS     = {'n_estimators': 200, 'n_jobs': 1, 'random_state': 0}

def count_kmers(fq, k):
    counts_by_acc = defaultdict(Counter)
    valid = set('ACGT')
    open_fn = gzip.open if fq.endswith(".gz") else open
    with open_fn(fq, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                acc = line.strip().split()[0][1:]
                run_id = acc.split('.')[0]
            elif i % 4 == 1:
                seq = line.strip().upper()
                for j in range(len(seq) - k + 1):
                    kmer = seq[j:j + k]
                    if all(c in valid for c in kmer):
                        counts_by_acc[run_id][kmer] += 1
    for run_id, counter in counts_by_acc.items():
        total = sum(counter.values())
        if total > 0:
            for kmer in counter:
                counter[kmer] /= total
    return counts_by_acc

def merge_counts(dicts, names):
    all_kmers = sorted(set(k for d in dicts for k in d))
    df = pd.DataFrame(0.0, index=names, columns=all_kmers)
    for n, d in zip(names, dicts):
        for kmer, freq in d.items():
            df.at[n, kmer] = freq
    return df

def compute_purity(linkage_matrix, labels):
    tree, _ = to_tree(linkage_matrix, rd=True)
    id_to_label = {i: labels[i] for i in range(len(labels))}
    label_to_indices = defaultdict(list)
    for i, label in enumerate(labels):
        label_to_indices[label].append(i)
    class_purities = {}
    for label, indices in label_to_indices.items():
        indices_set = set(indices)
        def find_subtree(node):
            if node.is_leaf():
                return node if node.id in indices_set else None
            left = find_subtree(node.left)
            right = find_subtree(node.right)
            if left and right:
                return node
            return left or right
        subtree = find_subtree(tree)
        if subtree:
            leaf_ids = subtree.pre_order()
            all_labels = [id_to_label[i] for i in leaf_ids]
            count = Counter(all_labels)
            purity = count[label] / len(all_labels)
            class_purities[label] = round(purity, 4)
        else:
            class_purities[label] = 0.0
    return np.mean(list(class_purities.values()))

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
            'f1': f1_score(y[te], pred, average='weighted'),
            'precision': precision_score(y[te], pred, average='weighted'),
            'recall': recall_score(y[te], pred, average='weighted'),
            'accuracy': accuracy_score(y[te], pred)
        }
        for c in classes:
            true_c = (y[te] == c)
            pred_c = (pred == c)
            m[f'pc_{c}_f1'] = f1_score(true_c, pred_c, zero_division=0)
            m[f'pc_{c}_precision'] = precision_score(true_c, pred_c, zero_division=0)
            m[f'pc_{c}_recall'] = recall_score(true_c, pred_c, zero_division=0)
        out.append(m)
    return out

def run_pipeline(label_array, accession_list, fq_path):
    cluster_recs, rf_recs, summary_recs = [], [], []
    for k in KMER_RANGE:
        print(f"\n>>> k={k} <<<")
        t0 = time.time()
        all_counts = count_kmers(fq_path, k)
        t_kmer = time.time() - t0
        X_df = merge_counts([all_counts[acc] for acc in accession_list], accession_list)
        nfeat = X_df.shape[1]
        X = X_df.values
        y = label_array
        mock_idx = np.where(y == "Mock")[0]
        selected_mock = np.random.choice(mock_idx, size=min(4, len(mock_idx)), replace=False)
        keep_idx = np.hstack([selected_mock] + [np.where(y == cls)[0] for cls in np.unique(y) if cls != "Mock"])
        X = X[keep_idx]
        y = y[keep_idx]
        cluster_reps = []
        for i in range(N_BOOT):
            idx = np.random.choice(len(X), int(len(X) * BOOTSTRAP_FRACTION), replace=False)
            X_sample = X[idx]
            y_sample = y[idx]
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_sample)
            linkage_matrix = linkage(X_scaled, method="ward")
            purity = compute_purity(linkage_matrix, y_sample)
            cluster_reps.append({'replicate': i, 'k': k, 'purity': purity})
        cluster_recs += cluster_reps
        t1 = time.time()
        rfolds = rf_fold_metrics(X, y, n_splits=N_SPLITS, rf_kwargs=RF_KWARGS)
        for r in rfolds:
            r.update({'k': k})
        rf_recs += rfolds
        t2 = time.time()
        summary_recs.append({
            'k': k,
            'n_features': nfeat,
            'time_kmer': t_kmer,
            'time_hclust': t1 - t0,
            'time_rf': t2 - t1
        })
    return cluster_recs, rf_recs, summary_recs

if __name__ == '__main__':
    md = pd.read_csv(METADATA_CSV, usecols=['Run', 'Virus']).set_index('Run')
    md.index = md.index.astype(str)
    fastq_files = sorted(glob.glob(os.path.join(ALL_FASTQ_DIR, '*.fastq')) +
                         glob.glob(os.path.join(ALL_FASTQ_DIR, '*.fastq.gz')))
    if not fastq_files:
        raise SystemExit("No FASTQ files found")

    def extract_accession_counts(fq_path):
        accession_counts = Counter()
        open_fn = gzip.open if fq_path.endswith(".gz") else open
        with open_fn(fq_path, 'rt') as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    acc = line.strip().split()[0][1:]
                    run_id = acc.split('.')[0]
                    accession_counts[run_id] += 1
        return accession_counts

    for fq in fastq_files:
        dataset_name = os.path.basename(fq).replace('.fastq', '').replace('.gz', '')
        print(f"\n=== Processing dataset: {dataset_name} ===")
        acc_counts = extract_accession_counts(fq)
        read_counts_df = pd.DataFrame.from_dict(acc_counts, orient='index', columns=['ReadCount'])
        read_counts_df.index.name = 'Run'
        read_counts_df.index = read_counts_df.index.astype(str)
        merged_df = md.join(read_counts_df, how='inner')
        if len(merged_df) == 0:
            print(f"No metadata match for dataset {dataset_name}, skipping.")
            continue
        labels_true = merged_df['Virus'].values
        accession_order = merged_df.index.tolist()
        cl_out, rf_out, sum_out = run_pipeline(labels_true, accession_order, fq)
        pd.DataFrame(cl_out).to_csv(f'cluster_replicates_{dataset_name}.csv', index=False)
        pd.DataFrame(rf_out).to_csv(f'rf_folds_{dataset_name}.csv', index=False)
        pd.DataFrame(sum_out).to_csv(f'k_timing_summary_{dataset_name}.csv', index=False)
        print(f"→ Wrote outputs for {dataset_name}")
