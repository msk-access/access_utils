import sys
import json
import argparse
import yaml
from pathlib import Path
from collections import OrderedDict

import pandas as pd


cols_not_rename = [
    'Sample', 'cmoSampleName', 'cmoSampleClass', 'cmoPatientId', 'investigatorSampleId',
    'barcodeId', 'oncoTreeCode', 'tumorOrNormal', 'tissueLocation', 'specimenType', 'sampleOrigin', 'sex',
    "status"]

SAMPLE_META_COLS = [
    "cmoSampleName",
    "igoId",
    "cmoSampleClass",
    "cmoPatientId",
    "investigatorSampleId",
    "barcodeId",
    "oncoTreeCode",
    "tumorOrNormal",
    "tissueLocation",
    "specimenType",
    "sampleOrigin",
    "sex",
    "dnaInputNg",
    "captureInputNg",
    "captureConcentrationNm",
    "status",
    "minor_contamination",
    "major_contamination",
    "sex_mismatch",
    "fingerprint",
    "raw_coverage_a",
    "raw_coverage_b",
    "duplex_target_coverage",
    "MODE_INSERT_SIZE",
    "TOTAL_READS",
    "PCT_PF_UQ_READS_ALIGNED",
    "noise_percentage",
    "noise_n_sites"
]


def get_args():
    parser = argparse.ArgumentParser(description='Prepare data for multiqc.')
    parser.add_argument('--dir', type=str, default=".", help='Directory containing results.')
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    samples = create_sample_meta_info(args)
    samples = parse_sequence_qc(args, samples)
    samples = parse_biometrics(args, samples)
    samples = parse_picard(args, samples)

    save_sample_meta(samples)


def save_sample_meta(samples):

    samples = pd.DataFrame(samples).T
    sample_col = samples.pop('cmoSampleName')
    samples.insert(0, 'cmoSampleName', sample_col)

    tumors = samples[samples['tumorOrNormal']=='Tumor']
    if len(tumors) > 0:
        tumors.columns = tumors.columns.map(lambda x: x + '_plasma' if x not in cols_not_rename else x)

        tumors.to_csv('genstats_qc_status_plasma.csv', index=False)

    normals = samples[samples['tumorOrNormal']!='Tumor']
    if len(normals) > 0:
        normals.columns = normals.columns.map(lambda x: x + '_buffy' if x not in cols_not_rename else x)
        normals.to_csv('genstats_qc_status_buffy.csv', index=False)


def create_sample_meta_info(args):

    files = list(Path(args.dir).rglob("sample_info.json"))
    if len(files) == 0:
        print('could not find sample_info.json file!')
        sys.exit(1)

    sample_info = json.load(open(files[0].resolve()))

    for i, sample_data in enumerate(sample_info):
        cols = set(SAMPLE_META_COLS).difference(set(sample_data.keys()))

        for col in cols:
            sample_info[i][col] = ""

    sample_info = {i['cmoSampleName']: i for i in sample_info}

    return sample_info


def write_html(html_path_in, html_path_out, section_name, description, parent_id, parent_name):

    html_file = open(html_path_in, 'r').read()
    header = "<!--\n" + \
        "section_name: '{}'\n".format(section_name) + \
        "description: '{}'\n".format(description) + \
        "parent_id: '{}'\n".format(parent_id) + \
        "parent_name: '{}'\n".format(parent_name) + \
        "-->\n"

    html_file = header + html_file

    fout = open(html_path_out, 'w')
    fout.write(html_file)
    fout.close()


def get_matching_sample(samples, text):
    for sample in samples:
        if text.find(sample) != -1:
            return sample

    return text


def parse_picard_file(fpath):

    parsed_data = dict()
    keys = None
    commadecimal = None
    filehandle = open(fpath)

    for l in filehandle:

        if "HsMetrics" in l and "## METRICS CLASS" in l:
            keys = filehandle.readline().strip("\n").split("\t")
        elif keys:
            vals = l.strip("\n").split("\t")
            if len(vals) == len(keys):
                j = "NA"
                if keys[0] == "BAIT_SET":
                    j = vals[0]
                parsed_data[j] = dict()
                # Check that we're not using commas for decimal places
                if commadecimal is None:
                    for i, k in enumerate(keys):
                        if k.startswith("PCT_"):
                            if "," in vals[i]:
                                commadecimal = True
                            else:
                                commadecimal = False
                for i, k in enumerate(keys):
                    try:
                        if commadecimal:
                            vals[i] = vals[i].replace(".", "")
                            vals[i] = vals[i].replace(",", ".")
                        parsed_data[j][k] = float(vals[i])
                    except ValueError:
                        parsed_data[j][k] = vals[i]

    return parsed_data


def parse_picard_insert_size_file(fpath):
    parsed_data = dict()
    keys = None
    filehandle = open(fpath)

    for l in filehandle:
        if "InsertSizeMetrics" in l and "## METRICS CLASS" in l:

            keys = filehandle.readline().strip("\n").split("\t")
            vals = filehandle.readline().strip("\n").split("\t")

            orientation_idx = keys.index("PAIR_ORIENTATION")
            while len(vals) == len(keys):
                pair_orientation = vals[orientation_idx]
                parsed_data[pair_orientation] = OrderedDict()

                for i, k in enumerate(keys):
                    try:
                        parsed_data[pair_orientation][k] = float(vals[i])
                    except ValueError:
                        try:
                            parsed_data[pair_orientation][k] = float(vals[i].replace(",", "."))
                            print(
                                "Switching commas for points in '{}': {} - {}".format(
                                    fpath, vals[i], vals[i].replace(",", ".")
                                )
                            )
                        except ValueError:
                            parsed_data[pair_orientation][k] = vals[i]
                    except IndexError:
                        pass  # missing data

                vals = filehandle.readline().strip("\n").split("\t")

            # Skip lines on to histogram
            l = filehandle.readline().strip("\n")
            l = filehandle.readline().strip("\n")

    return parsed_data


def parse_picard(args, samples):

    # parse uncollapsed BAM pool A picard metrics

    for fpath in Path(args.dir).rglob("uncollapsed_bam_stats_pool_a*/*/*hs_metrics.txt"):
        sample = get_matching_sample(list(samples.keys()), str(fpath))
        parsed_data = parse_picard_file(fpath)

        if sample in samples and parsed_data:
            samples[sample]['raw_coverage_a'] = \
                parsed_data['MSK-ACCESS-v1_0-probe-A_baits']['MEAN_TARGET_COVERAGE']

            samples[sample]['PCT_PF_UQ_READS_ALIGNED'] = \
                parsed_data['MSK-ACCESS-v1_0-probe-A_baits']['PCT_PF_UQ_READS_ALIGNED']
            samples[sample]['TOTAL_READS'] = \
                parsed_data['MSK-ACCESS-v1_0-probe-A_baits']['TOTAL_READS'] / 2

    # parse uncollapsed BAM pool B picard metrics

    for fpath in Path(args.dir).rglob("uncollapsed_bam_stats_pool_b*/*/*hs_metrics.txt"):
        sample = get_matching_sample(list(samples.keys()), str(fpath))
        parsed_data = parse_picard_file(fpath)

        if sample in samples and parsed_data:
            samples[sample]['raw_coverage_b'] = \
                parsed_data['MSK-ACCESS-v1_0-probe-B_baits']['MEAN_TARGET_COVERAGE']

    # parse duplex BAM pool A picard metrics

    for fpath in Path(args.dir).rglob("duplex_bam_stats_pool_a*/*/*hs_metrics.txt"):
        sample = get_matching_sample(list(samples.keys()), str(fpath))
        parsed_data = parse_picard_file(fpath)

        if sample in samples and parsed_data:
            samples[sample]['duplex_target_coverage'] = \
                parsed_data['MSK-ACCESS-v1_0-probe-A_baits'].get('MEAN_TARGET_COVERAGE')


    # parse picard insert metrics

    for fpath in Path(args.dir).rglob("duplex_bam_stats_pool_a*/*/*insert_size_metrics.txt"):
        sample = get_matching_sample(list(samples.keys()), str(fpath))
        parsed_data = parse_picard_insert_size_file(fpath)

        if sample in samples and parsed_data:
            if 'FR' in parsed_data:
                samples[sample]['MODE_INSERT_SIZE'] = \
                    parsed_data['FR'].get('MODE_INSERT_SIZE')

    return samples


def parse_biometrics(args, samples):

    # get minor/major contamination

    fpath = list(Path(args.dir).rglob("collapsed*/*major_contamination.csv"))

    if fpath:
        data = pd.read_csv(fpath[0])

        for i in data.index:
            sample = data.at[i, 'sample_name']
            val = data.at[i, 'major_contamination']

            if sample in samples:
                samples[sample]['major_contamination'] = val
            else:
                print('{} found in major contamination file, but not in sample meta info!'.format(sample))


    # get minor contamination

    fpath = list(Path(args.dir).rglob("collapsed*/*minor_contamination.csv"))

    if fpath:
        data = pd.read_csv(fpath[0])

        for i in data.index:
            sample = data.at[i, 'sample_name']
            val = data.at[i, 'minor_contamination']

            if sample in samples:
                samples[sample]['minor_contamination'] = val
            else:
                print('{} found in minor contamination file, but not in sample meta info!'.format(sample))

    # get sex mismatch

    fpath = list(Path(args.dir).rglob("collapsed*/*sex_mismatch.csv"))

    if fpath:
        data = pd.read_csv(fpath[0])

        for i in data.index:
            sample = data.at[i, 'sample']

            if sample in samples:
                samples[sample]['sex_mismatch'] = data.at[i, 'sex_mismatch']
            else:
                print('{} found in sex mismatch file, but not in sample meta info!'.format(sample))

    # get sex genotype

    fpath = list(Path(args.dir).rglob("collapsed*/*genotype_comparison.csv"))

    if fpath:
        data = pd.read_csv(fpath[0])

        for sample in samples.keys():
            data_sample = data[data['ReferenceSample']==sample]

            if len(data_sample) == 0:
                print('{} not found in biometrics genotype file!'.format(sample))
                sys.exit(1)

            sample_group = data_sample['ReferenceSampleGroup'].tolist()[0]
            statuses = data_sample['Status'].value_counts().to_dict()

            # check if there are any other samples from same patient

            if len(data_sample[data_sample['QuerySampleGroup']==sample_group]) == 1:
                samples[sample]['fingerprint'] = 'NA'
            elif (statuses.get('Unexpected Match', 0) + statuses.get('Unexpected Mismatch', 0)) > 0:
                samples[sample]['fingerprint'] = 'fail'
            else:
                samples[sample]['fingerprint'] = 'pass'

    # format html file

    files = list(Path(args.dir).rglob("collapsed*/*minor_contamination_sites.html"))

    if not files:
        return samples

    write_html(
        files[0].resolve(),
        'minor_contamination_sites_mqc.html',
        'Contributing sites',
        'Charts the sites that contribute to minor contamination.',
        'contamination',
        'Contamination'
    )

    return samples


def parse_sequence_qc(args, samples):

    sequence_qc_data = {}
    sequence_qc_substitution_data = {}

    for path in Path(args.dir).rglob("*noise_acgt.tsv"):
        f_data = pd.read_csv(path.resolve(), sep='\t')
        f_data = f_data.set_index('sample_id').to_dict()
        s_name = path.name.replace('noise_acgt.tsv', '')

        sequence_qc_data[s_name] = {}

        sequence_qc_data[s_name].update({
            'N_minor_alleles': f_data['minor_allele_count'].get(s_name),
            'N_major_allele': f_data['major_allele_count'].get(s_name),
            'noise_percentage': float(f_data['noise_fraction'].get(s_name)) * 100,
            'contributing_sites': f_data['contributing_sites'].get(s_name)
        })

        if s_name in samples:
            samples[s_name]['noise_percentage'] = sequence_qc_data[s_name]['noise_percentage']
            samples[s_name]['contributing_sites'] = sequence_qc_data[s_name]['contributing_sites']

    for path in Path(args.dir).rglob("*noise_del.tsv"):
        f_data = pd.read_csv(path.resolve(), sep='\t')
        f_data = f_data.set_index('sample_id').to_dict()
        s_name = path.name.replace('noise_del.tsv', '')

        sequence_qc_data[s_name].update({
            'n_deletions': f_data['del_count'].get(s_name),
            'total_base_count_del': f_data['total_base_count'].get(s_name),
            'noise_percentage_del': float(f_data['noise_fraction'].get(s_name))*100,
            'contributing_sites_del': f_data['contributing_sites'].get(s_name)
        })

    for path in Path(args.dir).rglob("*noise_n.tsv"):
        f_data = pd.read_csv(path.resolve(), sep='\t')
        f_data = f_data.set_index('sample_id').to_dict()
        s_name = path.name.replace('noise_n.tsv', '')

        sequence_qc_data[s_name].update({
            'ns': f_data['n_count'].get(s_name),
            'total_base_count_ns': f_data['total_base_count'].get(s_name),
            'noise_percentage_ns': float(f_data['noise_fraction'].get(s_name))*100,
            'contributing_sites_ns': f_data['contributing_sites'].get(s_name)
        })

    for path in Path(args.dir).rglob("*noise_positions.tsv"):
        f_data = pd.read_csv(path.resolve(), sep='\t')
        f_data['A'] = f_data['A']/f_data['total_acgt']
        f_data['C'] = f_data['C']/f_data['total_acgt']
        f_data['T'] = f_data['T']/f_data['total_acgt']
        f_data['G'] = f_data['G']/f_data['total_acgt']

        maf_c_a = pd.concat([f_data.loc[f_data['ref']=='C', 'A'], f_data.loc[f_data['ref']=='G', 'T']]).values
        maf_c_g = pd.concat([f_data.loc[f_data['ref']=='C', 'G'], f_data.loc[f_data['ref']=='G', 'C']]).values
        maf_c_t = pd.concat([f_data.loc[f_data['ref']=='C', 'T'], f_data.loc[f_data['ref']=='G', 'A']]).values
        maf_t_a = pd.concat([f_data.loc[f_data['ref']=='T', 'A'], f_data.loc[f_data['ref']=='A', 'T']]).values
        maf_t_c = pd.concat([f_data.loc[f_data['ref']=='T', 'C'], f_data.loc[f_data['ref']=='A', 'G']]).values
        maf_t_g = pd.concat([f_data.loc[f_data['ref']=='T', 'G'], f_data.loc[f_data['ref']=='A', 'C']]).values

        s_name = path.name.replace('noise_positions.tsv', '')
        sequence_qc_substitution_data[s_name] = {}

        sequence_qc_substitution_data[s_name].update({
            'C>A': float(maf_c_a[maf_c_a>0].mean()),
            'C>G': float(maf_c_g[maf_c_g>0].mean()),
            'C>T': float(maf_c_t[maf_c_t>0].mean()),
            'T>A': float(maf_t_a[maf_t_a>0].mean()),
            'T>C': float(maf_t_c[maf_t_c>0].mean()),
            'T>G': float(maf_t_g[maf_t_g>0].mean())
        })

        for sub in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']:
            if pd.isna(sequence_qc_substitution_data[s_name][sub]):
                sequence_qc_substitution_data[s_name][sub] = 0

    sequence_qc_output = {
        'id': 'sequence_qc_table',
        'section_name': 'Noise metrics',
        'description': 'Noise metrics calculated using sequence_qc for duplex BAM file.',
        'plot_type': 'table',
        'parent_id': 'sequence_qc',
        'parent_name': 'Duplex Noise Metrics',
        'pconfig': {
            'id': 'sequence_qc_table',
            'title': 'Per-sample noise metrics'
        },
        'data': sequence_qc_data,
        'headers': {
            'noise_percentage': {'format': '{:,.2}%', 'title': '% Noise (ACGT)'},
            'noise_percentage_del': {'format': '{:,.2}%', 'title': '% Noise (ACGT + Del)'},
            'noise_percentage_ns': {'format': '{:,.2}%', 'title': '% Noise (ACGT + Ns)'},
            'contributing_sites': {'format': '{:,.0f}', 'title': 'N contibuting sites (ACGT)'},
            'contributing_sites_ns': {'format': '{:,.0f}', 'title': 'N contibuting sites (ACGT + Ns)'},
            'contributing_sites_del': {'format': '{:,.0f}', 'title': 'N contibuting sites (ACGT + Del)'}
        }
    }

    fout = open('sequence_qc_mqc.yaml', 'w')
    yaml.dump(sequence_qc_output, fout)
    fout.close()

    sequence_qc_substitution_output = {
        'id': 'sequence_qc_substitution',
        'section_name': 'Duplex noise subsitution',
        'description': 'Noise metrics for each type of substitution. Calculated using sequence_qc from duplex BAM file(s).',
        'plot_type': 'bargraph',
        'parent_id': 'sequence_qc',
        'parent_name': 'Duplex Noise Metrics',
        'pconfig': {
            'id': 'sequence_qc_table_substitution',
            'title': 'Noise by substitution type',
            'stacking': None,
            "hide_zero_cats": False,
        },
        'data': sequence_qc_substitution_data
    }

    paths = list(Path(args.dir).rglob("*_noise.html"))

    if len(sequence_qc_substitution_data) > 1:
        sequence_qc_substitution_output['plot_type'] = 'table'
        sequence_qc_substitution_output['headers'] = {
            'C>A': {'format': '{:,.3}', 'color': '#3977af'},
            'C>G': {'format': '{:,.3}', 'color': '#3977af'},
            'C>T': {'format': '{:,.3}', 'color': '#3977af'},
            'T>A': {'format': '{:,.3}', 'color': '#3977af'},
            'T>C': {'format': '{:,.3}', 'color': '#3977af'},
            'T>G': {'format': '{:,.3}', 'color': '#3977af'}
        }

        fout = open('sequence_qc_substitution_mqc.yaml', 'w')
        yaml.dump(sequence_qc_substitution_output, fout)
        fout.close()
    elif len(paths) > 0:
        write_html(
            paths[0].resolve(),
            'sequence_qc_mqc.html',
            'Duplex noise figures',
            'Charts showing additional duplex noise information for a single sample.',
            'sequence_qc',
            'Duplex Noise'
        )

    return samples


if __name__ == "__main__":
    main()