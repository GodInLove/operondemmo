import datetime
import os
import matplotlib.pyplot as plt
import numpy
from operondemmo.hierarchical_cluster.gamma_domain import get_result_by_clustering, get_result_by_clustering2

from operondemmo.input_file_handle.handle_gff import get_gene_pos_strand, sorted_gene, \
    from_simple_gff_information_to_get
from operondemmo.operon import from_depth_file_to_get_co_matrix_co_expression, load_from_input_files


def compute_n(list_sort, p_or_n_predict):
    p_or_n_predict_fp = open(p_or_n_predict, 'r')
    p_or_n_predict_out = p_or_n_predict_fp.read().strip()
    p_or_n_predict_out = p_or_n_predict_out.split("\n")
    count_ = 0
    # print(len(p_or_n_predict_out))
    for item in p_or_n_predict_out:
        item = item.split(";")
        item = sorted(item)
        item = ";".join(item)
        if item in list_sort:
            count_ = count_ + 1
    return count_


def generate_output_files(simple_gff_file, depth_files, co_method, out_path, static_threshold_list):
    _gene_pos, _gene_strand = get_gene_pos_strand(simple_gff_file)
    gene_strand_m, gene_index_m, sorted_gene_m = from_simple_gff_information_to_get(_gene_pos, _gene_strand)
    matrix_co = from_depth_file_to_get_co_matrix_co_expression(depth_files, _gene_pos, co_method)
    i_iter = 0
    for each_t in static_threshold_list:
        if i_iter < 10:
            tmp_name = "0" + str(i_iter)
        else:
            tmp_name = str(i_iter)
        result_file = out_path + tmp_name + ".out"
        get_result_by_clustering2(result_file, gene_strand_m, gene_index_m, sorted_gene_m, matrix_co, each_t)
        i_iter = i_iter + 1


def get_list_from_file(file_path):
    f_p_fp = open(file_path)
    f_p_content = f_p_fp.read().strip()
    f_p_fp.close()
    f_p_out = f_p_content.split("\n")
    f_p_out_sort = []
    for item in f_p_out:
        item = item.split(";")
        item = sorted(item)
        item = ";".join(item)
        f_p_out_sort.append(item)
    return f_p_out_sort


def roc_curve(simple_gff, depth_files, method, out_path, p_d_list, n_d_list, tp_fp_file, roc_file_path, threshold_list):
    generate_output_files(simple_gff, depth_files, method, out_path, threshold_list)
    list_files = os.listdir(out_path)
    list_files = sorted(list_files)
    tp_fp_file_fp = open(tp_fp_file, 'w')
    tpr_fpr_file_fp = open(roc_file_path, 'w')
    p = len(p_d_list)
    n = len(n_d_list)
    for _file in list_files:
        print(_file)
        tp = compute_n(p_d_sort_list, out_path + _file)
        fp = compute_n(n_d_sort_list, out_path + _file)
        tp_fp_file_fp.write(str(tp) + "\t" + str(fp) + "\n")
        print(tp, fp)
        tpr = tp / p
        fpr = fp / n
        print(tpr, fpr)
        tpr_fpr_file_fp.write(str(tpr) + "\t" + str(fpr) + "\n")
    tp_fp_file_fp.close()
    tpr_fpr_file_fp.close()
    matrix_a = numpy.loadtxt(roc_file_path)
    tpr = matrix_a[..., 0].tolist()
    fpr = matrix_a[..., 1].tolist()
    return tpr, fpr


def get_t_list(start, stop, n):
    sep = (stop - start) / n
    t_list = []
    i = float(start)
    while i <= stop:
        t_list.append(i)
        i = i + sep
    return t_list


if __name__ == "__main__":
    path = "/home/lyd/document/2018.1/gamma_domain/"
    eco_simple_gff = path + "simple_eco.gff_3"
    eco_depth_files = load_from_input_files(path + "eco_count/")
    tmp_out_path = path + "g_d_out/"
    co_expression_method = 0
    roc_path = "/home/lyd/document/2018.1/roc/"
    positive_data = roc_path + "positiveData.txt"
    negative_data = roc_path + "negativeData.txt"
    tmp_tp_fp_file = path + "tp_fp.txt"
    tmp_roc_file_path = path + "tpr_fpr.txt"
    threshold_list = get_t_list(-1, 1, 50)
    p_d_sort_list = get_list_from_file(positive_data)
    n_d_sort_list = get_list_from_file(negative_data)
    tpr_0, fpr_0 = roc_curve(eco_simple_gff, eco_depth_files, co_expression_method, tmp_out_path,
                             p_d_sort_list, n_d_sort_list, tmp_tp_fp_file, tmp_roc_file_path, threshold_list)

    tmp_out_path = path + "g_d_out_1/"
    co_expression_method = 1
    tmp_tp_fp_file = path + "tp_fp_1.txt"
    tmp_roc_file_path = path + "tpr_fpr_1.txt"
    tpr_1, fpr_1 = roc_curve(eco_simple_gff, eco_depth_files, co_expression_method, tmp_out_path,
                             p_d_sort_list, n_d_sort_list, tmp_tp_fp_file, tmp_roc_file_path, threshold_list)

    tmp_out_path = path + "g_d_out_2/"
    co_expression_method = 2
    tmp_tp_fp_file = path + "tp_fp_2.txt"
    tmp_roc_file_path = path + "tpr_fpr_2.txt"
    tpr_2, fpr_2 = roc_curve(eco_simple_gff, eco_depth_files, co_expression_method, tmp_out_path,
                             p_d_sort_list, n_d_sort_list, tmp_tp_fp_file, tmp_roc_file_path, threshold_list)
    plt.title("ROC curve of gamma_domain algorithm\n(AUC="")")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.xlim(0)
    plt.ylim(0)
    plt.plot(fpr_0, tpr_0)
    plt.plot(fpr_1, tpr_1)
    plt.plot(fpr_2, tpr_2)
    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
    plt.show()
