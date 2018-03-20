import itertools

from operondemmo.input_file_handle.handle_gff import from_simple_gff_information_to_get, get_gene_pos_strand, \
    sorted_gene


def get_true_tu(p_d_f):
    p_d_fp = open(p_d_f, 'r')
    p_d_data = p_d_fp.read().strip()
    p_d_data = p_d_data.split("\n")
    t_tu = []
    for each in p_d_data:
        each_list = each.split(";")
        each_list = sorted(each_list)
        each_list = tuple(each_list)
        t_tu.append(each_list)
    t_tu = set(t_tu)
    return t_tu


def combination_tu(locus_list, c_num, tu_len):
    len_locus_list = len(locus_list)
    _i = 0
    tmp_list = set()
    while _i < len_locus_list - c_num + 1:
        splice_list = locus_list[_i:_i + c_num]
        c_list = itertools.combinations(splice_list, tu_len)
        c_list = set(c_list)
        cc_list = []
        for each_c in c_list:
            each_c = list(each_c)
            each_c = sorted(each_c)
            each_c = tuple(each_c)
            cc_list.append(each_c)
        cc_set = set(cc_list)
        tmp_list = tmp_list.union(cc_set)
        _i = _i + 1
    return tmp_list


def add_name(_gene_name, result_):
    name_result = []
    for each_result in result_:
        name_result.append(_gene_name[each_result])
    name_result = tuple(name_result)
    return name_result


if __name__ == "__main__":

    # test_lcs_list = ['b0001', 'b0009', 'b0003', 'b0008', 'b0005', 'b0006', 'b0007']
    # test_tu_set = combination_tu(test_lcs_list, 9, 8)
    # print(test_tu_set)
    test_gff_path = "/home/lyd/document/2018.1/gamma_domain/simple_eco.gff_3"
    test_gene_pos, test_gene_strand = get_gene_pos_strand(test_gff_path)
    gene_name = sorted_gene(test_gene_pos)
    print(gene_name, len(gene_name))
    _gene_strand, _gene_index, _gene_sort = \
        from_simple_gff_information_to_get(test_gene_pos, test_gene_strand)
    positive_data_file = "/home/lyd/document/2018.1/roc/positiveData.txt"
    true_tu_set = get_true_tu(positive_data_file)
    # print(true_tu)
    f_tu_set = set()
    # print(_gene_index)
    j = 0
    while j < len(_gene_index):
        each_ = _gene_index[j]
        tmp_len = len(each_)
        # print(each_, tmp_len)
        i = 1
        while i <= tmp_len and i <= 10:
            tmp_tu_set = combination_tu(each_, i + 1, i)
            f_tu_set = f_tu_set.union(tmp_tu_set)
            i = i + 1
        j = j + 1
    name_set = set()
    for each_tu in f_tu_set:
        tmp_tu = add_name(gene_name, each_tu)
        name_set.add(tmp_tu)
    name_new_set = name_set.difference(true_tu_set)
    negative_data_file = "/home/lyd/document/2018.1/roc/negativeData.txt"
    n_d_fp = open(negative_data_file, 'w')
    for each_tu in name_new_set:
        list_each_tu = list(each_tu)
        n_d_fp.write(";".join(list_each_tu) + "\n")
    n_d_fp.close()
    print(len(name_new_set))
    # 1403311
