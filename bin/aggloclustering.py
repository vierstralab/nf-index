#!/bin/usr/python

#Import Libraries
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, pair_confusion_matrix, homogeneity_score, completeness_score
from sklearn.cluster import AgglomerativeClustering
import json
from sklearn.preprocessing import LabelEncoder
import pickle


def get_clustering_metrics(labels_pred, labels_true):
    labels_pred, labels_true = zip(*[(x, y) for x, y in zip(labels_pred, labels_true) if not pd.isna(y)])
    labels1 = np.array(labels_pred)
    labels2 = np.array(labels_true)
    
    assert len(labels1) == len(labels2)

    #Pair confusion  matrix
    N = pair_confusion_matrix(labels2, labels1) // 2
    
    #Adjusted Rand Index
    ari = 2 * (N[0, 0]*N[1, 1] - N[0, 1]*N[1, 0]) / ((N[0, 0] + N[0, 1])*(N[0, 1] + N[1, 1]) + (N[0, 0] + N[1, 0])*(N[1, 0] + N[1, 1]))
    #assert ari == adjusted_rand_score(labels2, labels1)
    
    #Fowlkes-Mallows score
    fm = np.sqrt((N[1, 1] / (N[1, 1] + N[0, 1]))*(N[1, 1] / (N[1, 1] + N[1, 0])))
    #assert fm == fowlkes_mallows_score(labels2, labels1)
    
    #Jaccard simmilarity
    jaccard = N[1, 1] / (N[1, 0] + N[0, 1] + N[1, 1])
    
    #Adjusted Mutual Information
    ami = adjusted_mutual_info_score(labels2, labels1)
    
    #Homogeniety and Completeness (not adjusted to by-chance values)
    hom = homogeneity_score(labels2, labels1)
    comp = completeness_score(labels2, labels1)
    
    return ami, ari, fm, hom, comp, jaccard, N

def main(json_object, meta, embedding):
    #Run Clustering Algorithm
    #Methods = single, average, weighted, centroid, median, and ward
    #Metrics = euclidian, manhattan, hamming, minkowski
    clustering_model = AgglomerativeClustering(*json_object)
    clustering_model.fit(embedding)
    clustered_labels =  clustering_model.labels_


    #Encode Ontology_Term
    label_encoder = LabelEncoder()
    true_labels = label_encoder.fit_transform(meta['ontology_term'])

    #Get Clustering Metrics
    metrics = get_clustering_metrics(clustered_labels, true_labels)

    #Save Metrics
    metrics_df = pd.DataFrame([metrics[0], metrics[1], metrics[2], metrics[5]]).T
    
    #Save Cluster Annotations
    annotations = pd.DataFrame(clustered_labels)
    return metrics_df, annotations, clustering_model


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as openfile:
        params = json.load(openfile)
    embedding = np.load(sys.argv[2])
    #Load Metadata
    meta_columns = ['tc_number', 'ds_number', 'ln_number', 'id', 
                'passed_agg_qc', 'hotspot1_spot', 'hotspot2_spot', 'nuclear_percent_duplication', 'taxonomy_name', 
                'ontology_id', 'ontology_term', 'system', 'subsystem', 'organ', 'organ_region', 'side_position', 
                'cell_type', 'cell_markers', 'disease_pathology', 'sampling_origin', 'state', 'germ_layer', 
                'growth_stage', 'description', 'donor_id', 'donor_age', 'donor_age_units', 'donor_sex', 
                'donor_ethnicity', 'donor_health_condition', 'vendor_name', 'perturbation', 'treatment_name', 
                'treatment_details', 'dose', 'dose_units', 'filtered_alignments_bam', 'bam_index', 'hotspot_peaks', 
                'hotspot_peaks_1per', 'hotspot_peaks_point1per', 'adapter', 'total_alignments', 'all_aligned', 
                'all_mapq_filter', 'autosomal_align', 'chr_m', 'hotspot1_num_bases', 'hotspot1_num_spots', 
                'hotspot2_num_bases', 'hotspot2_num_spots', 'insert_ft_eleven', 'insert_ls_ratio', 'paired_aligned', 
                'paired_autosomal_align', 'paired_nuclear_align', 'picard_median_insert_size', 'picard_percent_duplication',
                'raj_dataset', 'tag', 'footprintdir', 'num_footprints']
    
    meta = pd.read_table(sys.argv[3], header=None, names=meta_columns)
    prefix=sys.argv[4]

    metrics_df, annotations, clustering_model = main(params, meta, embedding)

    annotations.to_csv(f"{prefix}.annotations.txt", header=None, index=False)

    #Save Model
    with open(f"{prefix}.model.pkl") as out:
        pickle.dump(clustering_model, out)

    metrics_df.to_csv(f'{prefix}.metrics.tsv', sep="\t", header=None, index=False)
