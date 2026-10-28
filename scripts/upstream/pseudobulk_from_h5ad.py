"""
Build GSE290679 donor × lineage pseudobulk counts and metadata from the author-annotated h5ad.

This is an upstream preprocessing step. The reproducibility package begins at the
resulting pseudobulk matrix; it does not rerun this script automatically.

Usage:
    python pseudobulk_from_h5ad.py --input-h5ad path/to/GSE290679_raw_labeled_filtered.h5ad --output-dir path/to/output

Output files (names match config.example.json and Required input inventory.tsv):
    - GSE290679_pseudobulk_counts.tsv.gz
    - GSE290679_pseudobulk_metadata.csv

Requirements: h5py, numpy, scipy
"""
from pathlib import Path
import argparse, csv, gzip, json, collections
import numpy as np, h5py, scipy.sparse as sp


def decode(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def vector(g):
    if isinstance(g, h5py.Dataset):
        return decode(g[:])
    cats = decode(g['categories'][:])
    return [cats[c] if c >= 0 else '' for c in g['codes'][:]]


def lineage(s):
    if s in ['CD14 Mono', 'CD16 Mono']:
        return 'Monocyte'
    if s.startswith('CD4 T') or s == 'CD4 CTL':
        return 'CD4 T'
    if s.startswith('CD8 T'):
        return 'CD8 T'
    # Original GSE279452 kept Plasma separate from B: do not fold plasmablasts into B.
    if s in ['B Naive', 'B Transitional', 'B Memory']:
        return 'B'
    if s in ['NK', 'NK CD56 Bright']:
        return 'NK'
    return ''


def main():
    parser = argparse.ArgumentParser(description='Build GSE290679 pseudobulk from h5ad')
    parser.add_argument('--input-h5ad', required=True, help='Path to author-annotated GSE290679 h5ad file')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    args = parser.parse_args()

    h5ad_path = Path(args.input_h5ad)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not h5ad_path.exists():
        raise FileNotFoundError(f'Input h5ad not found: {h5ad_path}')

    with h5py.File(h5ad_path, 'r') as f:
        obs = f['obs']
        sample = vector(obs['orig.ident'])
        group = vector(obs['Group'])
        fine = vector(obs['celltype_fine'])
        barcodes = vector(obs['_index'])
        lin = [lineage(s) for s in fine]

        assert len(set(barcodes)) == len(barcodes)

        dg = {}
        for s, g in zip(sample, group):
            assert s not in dg or dg[s] == g
            dg[s] = g
        assert collections.Counter(dg.values()) == {'CI-Sep': 19, 'CI-NS': 19, 'NHC': 9}

        clinical = {k: vector(obs[k]) for k in ['SCARAB_ID', 'Age', 'Biologic_Sex', 'Source', 'ICU', 'Batch'] if k in obs}

        donor_rows = []
        for s in sorted(dg):
            ii = [i for i, x in enumerate(sample) if x == s]
            row = {'sample_id': s, 'group': dg[s]}
            for k, vv in clinical.items():
                values = sorted(set(vv[i] for i in ii))
                assert len(values) == 1, (s, k, values)
                row[k] = values[0]
            donor_rows.append(row)

        scarab = [r.get('SCARAB_ID', '') for r in donor_rows if r['group'] != 'NHC']
        assert all(x not in ['', 'nan', 'None', 'NA'] for x in scarab) and len(set(scarab)) == 38, \
            ('ICU independent donor identity gate failed', scarab)
        assert len(set(r['SCARAB_ID'] for r in donor_rows)) == 47

        ages = [float(r['Age']) for r in donor_rows]
        assert all(np.isfinite(a) and a >= 18 for a in ages), ('Adult audit gate', ages)

        keys = sorted(set((s, l) for s, l in zip(sample, lin) if l))
        lookup = {k: i for i, k in enumerate(keys)}
        ids = np.array([lookup.get((s, l), -1) for s, l in zip(sample, lin)])
        counts = np.bincount(ids[ids >= 0], minlength=len(keys))

        v = f['var']
        vidx = vector(v['_index'])
        features = vidx
        for candidate in ['gene_symbols', 'gene_symbol', 'feature_name', 'gene_name', 'Gene']:
            if candidate in v:
                features = vector(v[candidate])
                break

        shape = tuple(f['X'].attrs['shape'])
        assert shape == (len(sample), len(features))

        X = f['X']
        encoding = X.attrs['encoding-type']
        encoding = encoding.decode() if isinstance(encoding, bytes) else encoding
        assert encoding == 'csr_matrix', ('Unsupported format, inspect rather than transpose blindly', encoding)

        indptr = X['indptr'][:]
        result = np.zeros((len(features), len(keys)), dtype=np.float64)

        for start in range(0, len(sample), 10000):
            end = min(start + 10000, len(sample))
            a = int(indptr[start])
            b = int(indptr[end])
            data = X['data'][a:b]
            idx = X['indices'][a:b]
            assert np.isfinite(data).all() and (data >= 0).all() and np.max(np.abs(data - np.round(data))) < 1e-6
            mat = sp.csr_matrix((data, idx, indptr[start:end+1] - a), shape=(end - start, len(features)))
            ii = ids[start:end]
            valid = np.where(ii >= 0)[0]
            indicator = sp.csr_matrix(
                (np.ones(len(valid)), (valid, ii[valid])),
                shape=(end - start, len(keys))
            )
            result += (mat.T @ indicator).toarray()
            if start % 100000 == 0:
                print('aggregated', end, len(sample), flush=True)

        if len(set(features)) < len(features):
            unique = list(dict.fromkeys(features))
            mapping = {x: i for i, x in enumerate(unique)}
            new = np.zeros((len(unique), len(keys)))
            np.add.at(new, np.array([mapping[x] for x in features]), result)
            result = new
            features = unique

        counts_path = out_dir / 'GSE290679_pseudobulk_counts.tsv.gz'
        with gzip.open(counts_path, 'wt', encoding='utf-8', newline='') as out:
            w = csv.writer(out, delimiter='\t')
            w.writerow(['gene_symbol'] + ['||'.join(k) for k in keys])
            w.writerows([g] + [int(x) for x in result[i]] for i, g in enumerate(features))

        metadata_path = out_dir / 'GSE290679_pseudobulk_metadata.csv'
        clinical_by_sample = {r['sample_id']: r for r in donor_rows}
        fields = ['pseudobulk_id', 'sample_id', 'donor_id', 'cell_type', 'group', 'cell_count',
                  'Age', 'Biologic_Sex', 'Source', 'ICU', 'Batch']
        with open(metadata_path, 'w', encoding='utf-8', newline='') as out:
            w = csv.DictWriter(out, fieldnames=fields)
            w.writeheader()
            for i, k in enumerate(keys):
                dr = clinical_by_sample[k[0]]
                row = dict(
                    pseudobulk_id='||'.join(k),
                    sample_id=k[0],
                    donor_id=dr['SCARAB_ID'],
                    cell_type=k[1],
                    group=dg[k[0]],
                    cell_count=int(counts[i])
                )
                row.update({name: dr.get(name, '') for name in fields[6:]})
                w.writerow(row)

    print('PSEUDOBULK_COMPLETE', result.shape, flush=True)
    print('Outputs:', counts_path, metadata_path)


if __name__ == '__main__':
    main()
