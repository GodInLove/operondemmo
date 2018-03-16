import sys
import argparse

import os

from operondemmo.input_file_handle.handle_gff import auto_download, generate_simple_gff, sorted_gene, \
    from_simple_gff_to_get_list_gene_strand, get_gene_pos_strand, partition_gene

self_version = "0.0"

APP_VERSION = (
        '''
    ----------------------------------------------------------------------
    
    operondemmo-(%s) - an independent demmo of KNOWN operon predict method
    
    ----------------------------------------------------------------------
    ''' % self_version
)


def main():
    starting(prepare(sys.argv))


def prepare(argv):
    parser = argparse.ArgumentParser(description=APP_VERSION,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    advanced_argv = parser.add_argument_group("ADVANCED OPTIONS")
    parser.add_argument("-i", action="store", dest="input_files", default="null",
                        help="A directory to store a group of result files "
                             "through [samtools depth XXX > xxx.txt] command")
    parser.add_argument("-o", action="store", dest="output_path", default="OUT",
                        help="A directory include output data(operon file).default:OUT")
    parser.add_argument("-g", action="store", dest="gff_file", default="null",
                        help="The gff file of the prokaryote")
    parser.add_argument("-p", action="store", dest="process_thread", default=1, type=int,
                        help="Specify the number of processing threads (CPUs).default:1")
    parser.add_argument("-t", action="store", dest="threshold", default=0.6, type=int,
                        help="")
    advanced_argv.add_argument("-k", action="store", dest="kegg_id", default="null",
                               help="The kegg id of the prokaryote")
    advanced_argv.add_argument("--auto_gff", action="store_true", dest="auto_gff", default=False,
                               help="Auto download gff_file from NCBI Database")
    advanced_argv.add_argument("--person", action="store_true", dest="person", default=False,
                               help="Build co-expression matrix with person correlation")
    advanced_argv.add_argument("--spearman", action="store_true", dest="spearman", default=False,
                               help="Build co-expression matrix with spearman correlation")
    advanced_argv.add_argument("-v", "--version", action="version", version="operondemmo-" + self_version)
    if len(argv) == 1:
        print(parser.print_help())
        sys.exit(0)
    else:
        args = parser.parse_args(argv[1:])
        return args


def starting(args):
    """

    :type args: parser.parse_args()
    """
    if args.auto_gff:
        if args.kegg_id != "null":
            gff_file_path = auto_download(args.kegg_id)
        else:
            print("NEEDED KEGG ID.PLEASE check your input with option '-k'")
            return
    else:
        if args.kegg_id != "null":
            gff_file_path = args.gff_file
        else:
            print("NEEDED GFF_FILE.PLEASE check your input with option '-g'")
            return
    if args.person:
        co_expression_method = 1
    else:
        if args.spearman:
            co_expression_method = 2
        else:
            co_expression_method = 0
    if args.input_files != "null":
        input_files = os.listdir(args.input_files)
    else:
        print("NEEDED INPUT_FILES.PLEASE check your input with option '-i'")
        return
    if args.threshold > 1 or args.threshold < -1:
        print("IT CANNOT BE:", args.threshold, "PLEASE check your input with option '-t'")
        return
    if args.output_path[-1] != "/":
        output_path = args.output_path + "/"
    else:
        output_path = args.output_path
    operon_predict(args.threshold, input_files, output_path, gff_file_path, args.process_thread,
                   co_expression_method)


def operon_predict(threshold, input_files, output_path, gff_file_path, p, co_expression_method):
    gff_file = gff_file_path.split("/")
    gff_file = gff_file[-1]
    tmp_path = output_path + "tmp/"
    simple_gff_path = tmp_path + "simple_" + gff_file
    generate_simple_gff(gff_file_path, simple_gff_path)
    depth_files = []
    for each_file in os.listdir(input_files):
        depth_files.append(input_files + each_file)
    matrix_i_j = from_depth_file_to_get_co_matrix_co_expression(depth_files, simple_gff_path, co_expression_method)
    gene_sort = sorted_gene(simple_gff_path)
    list_strand = from_simple_gff_to_get_list_gene_strand(simple_gff_path)
    gene_pos, gene_strand = get_gene_pos_strand(simple_gff_path)
    final_gene_strand, final_gene_index, final_gene_sort = partition_gene(gene_sort, list_strand)
    result_file = output_path + "operon.txt"
    get_result_file(result_file, final_gene_strand, final_gene_index, final_gene_sort, matrix_i_j, threshold)


def from_depth_file_to_get_co_matrix_co_expression(depth_files, simple_gff_file, method):
    matrix_a = 0
    return matrix_a


def get_result_file(result_file, final_gene_strand, final_gene_index, final_gene_sort, matrix_i_j, threshold):
    pass
