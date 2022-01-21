import os, gzip, re, itertools, pickle, bisect, array, sys, glob
import concurrent.futures
import argparse
from signal import signal, SIGPIPE, SIG_DFL

sys.path.append(os.path.join(os.path.dirname(__file__), '.'))
import rw_table

global verbose
verbose = False

def read_gff(file):
    if re.search(r'.(gff|gff3)$', file):
        with open(file) as f:
            data = f.read()
    elif re.search(r'.(gff|gff3).gz$', file):
        with gzip.open(file, 'rt') as f:
            data = f.read()
    else:
        raise Exception("GFF extension is not *.gff/.gff3(.gz)")
    
    feature_list = []
    line_cnt = 0
    for line in data.splitlines():
        line_cnt += 1
        if line == '' or line[0] == '#':
            continue

        seqid, source, type, start, end, score, strand, phase, attributes = line.split('\t')

        id = parent = None
        result = re.findall(r'ID=(.+?)(;|$)', attributes)
        if result:
            id = result[0][0]
        
        result = re.findall(r'Parent=(.+?)(;|$)', attributes)
        if result:
            parent = result[0][0]
        
        if id == None and parent == None:
            if verbose:
                print(f'Warning: ID/Parent tag not found. Line {line_cnt}\nSkipping ...\t{line}', file=sys.stderr)
            continue
        
        feature = {
            'seqid': seqid,
            'source': source,
            'type': type,
            'start': int(start),
            'end': int(end),
            'score': score,
            'strand': strand,
            'phase': phase,
            #'attributes': attributes,
            'id': id,
            'parent': parent
        }
        feature_list.append(feature)

    return feature_list


def read_fasta(file):
    if re.search(r'\.(fa|fas|fasta)$', file):
        with open(file) as f:
            data = f.read()
    elif re.search(r'\.(fa|fas|fasta).gz$', file):
        with gzip.open(file, 'rt') as f:
            data = f.read()
    else:
        raise Exception("FASTA extension is not *.fa/.fas/.fasta(.gz)")

    seqid = ''
    sequence = ''
    sequence_table = {}
    for line in data.splitlines():
        if line == '':
            continue
        elif line[0] == '>':
            if seqid != '':
                sequence_table[seqid] = sequence
                sequence = ''
            result = re.findall(r'>(\S+)', line)
            seqid = result[0]
        else:
            sequence += line
    if seqid != '':
        sequence_table[seqid] = sequence

    return sequence_table


def create_kmer_list(kmer_size):
    bases = ['A', 'T', 'G', 'C']
    return [''.join(p) for p in itertools.product(bases, repeat=kmer_size)]


def scan_kmer_coord(sequence, kmer_size):
    sequence = sequence.translate(sequence.maketrans('atgcnuU', 'ATGCNTT'))
    kmer_coord_list = {}
    for i in range(len(sequence)-kmer_size+1):
        kmer = sequence[i:i+kmer_size]
        if 'N' in kmer:
            continue
        if kmer in kmer_coord_list:
            kmer_coord_list[kmer].append(i+1)   # 1-based coord
        else:
            kmer_coord_list[kmer] = array.array('I', [i+1])   # 1-based coord
    return kmer_coord_list


def reverse_complement(sequence):
    complement = {'A':'T', 'C':'G','G': 'C','T':'A', 'N':'N'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
    return reverse_complement


def write_kmer_table_by_sequence(table_dir_path, seqid, kmer, coord_list, format):
    if not format:
        format = 'csv'
    if format == 'pkl':
        if not os.path.isfile(f"{table_dir_path}/{kmer}.pkl.gz"):
            with gzip.open(f"{table_dir_path}/{kmer}.pkl.gz", mode='wb') as f:
                pickle.dump({}, f)
        with gzip.open(f"{table_dir_path}/{kmer}.pkl.gz", mode='rb') as f:
            current_kmer_seqid_coord_list = pickle.load(f)
        current_kmer_seqid_coord_list[seqid] = coord_list
        with gzip.open(f"{table_dir_path}/{kmer}.pkl.gz", mode='wb') as f:
            pickle.dump(current_kmer_seqid_coord_list, f)
    elif format == 'csv':
        coord_list_str = ','.join([str(c) for c in coord_list])
        with gzip.open(f"{table_dir_path}/{kmer}.csv.gz", mode='at') as f:
            f.write(f'{seqid},{coord_list_str}\n')


def read_kmer_table(table_dir_path, kmer):
    if os.path.isfile(f"{table_dir_path}/{kmer}.pkl.gz"):
        with gzip.open(f"{table_dir_path}/{kmer}.pkl.gz", mode='rb') as f:
            target_kmer_coord_table = pickle.load(f)
        return target_kmer_coord_table
    elif os.path.isfile(f"{table_dir_path}/{kmer}.csv.gz"):
        target_kmer_coord_table = {}
        with gzip.open(f"{table_dir_path}/{kmer}.csv.gz", mode='rt') as f:
            for line in f.readlines():
                tmp_list = line.split(',')
                seqid = tmp_list.pop(0)
                if tmp_list[0] == '\n':
                    coord_list = array.array('I', [])
                else:
                    coord_list = array.array('I', [int(x) for x in tmp_list])
                target_kmer_coord_table[seqid] = coord_list
        return target_kmer_coord_table
    raise Exception(f"Not found table: {table_dir_path}/{kmer}")


def generate_kmer_coord_table_by_sequence(fasta_file, kmer_size, out_dir, format):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    sequence_table = read_fasta(fasta_file)
    for seqid, sequence in sequence_table.items():
        kmer_coord_list = scan_kmer_coord(sequence, kmer_size)
        kmer_list = [kmer for kmer in kmer_coord_list.keys()]
        coord_list_list = [coord_list for coord_list in kmer_coord_list.values()]
        seqid_list = [seqid for i in range(len(kmer_list))]
        out_dir_list = [out_dir for i in range(len(kmer_list))]
        format_list = [format for i in range(len(kmer_list))]
        with concurrent.futures.ProcessPoolExecutor() as executer:
            executer.map(write_kmer_table_by_sequence, out_dir_list, seqid_list, kmer_list, coord_list_list, format_list)


def write_kmer_table(table_dir_path, kmer, seqs_coord_list, seqid_list, format):
    if not format:
        format = 'csv'
    if format == 'pkl':
        tmp_seqs_coord_list = {}
        for seqid in seqid_list:
            if seqid in seqs_coord_list:
                tmp_seqs_coord_list[seqid] = seqs_coord_list[seqid]
            else:
                tmp_seqs_coord_list[seqid] = array.array('I', [])
        seqs_coord_list = tmp_seqs_coord_list
        if not os.path.isfile(f"{table_dir_path}/{kmer}.pkl.gz"):
            with gzip.open(f"{table_dir_path}/{kmer}.pkl.gz", mode='wb') as f:
                pickle.dump({}, f)
        with gzip.open(f"{table_dir_path}/{kmer}.pkl.gz", mode='rb') as f:
            current_seqs_coord_list = pickle.load(f)
        new_seqs_coord_list = {**current_seqs_coord_list, **seqs_coord_list}
        with gzip.open(f"{table_dir_path}/{kmer}.pkl.gz", mode='wb') as f:
            pickle.dump(new_seqs_coord_list, f)
    elif format == 'csv':
        with gzip.open(f"{table_dir_path}/{kmer}.csv.gz", mode='at') as f:
            for seqid in seqid_list:
                if seqid in seqs_coord_list:
                    coord_list_str = ','.join([str(c) for c in seqs_coord_list[seqid]])
                    f.write(f'{seqid},{coord_list_str}\n')
                else:
                    f.write(f'{seqid},\n')


def generate_kmer_coord_table(fasta_file, kmer_size, out_dir, format):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    kmer_list = create_kmer_list(kmer_size)
    sequence_table = read_fasta(fasta_file)
    all_kmer_list = {v: {} for v in kmer_list}

    if len(sequence_table.keys()) > 10000:
        for seqid, sequence in sequence_table.items():
            kmer_coord_list = scan_kmer_coord(sequence, kmer_size)
            for kmer, coord in kmer_coord_list.items():
                all_kmer_list[kmer][seqid] = coord
    else:
        sequence_list = sequence_table.values()
        kmer_size_list = [kmer_size for i in range(len(sequence_list))]
        with concurrent.futures.ProcessPoolExecutor() as executer:
            kmer_coord_list_by_sequence = executer.map(scan_kmer_coord, sequence_list, kmer_size_list)

        for kmer_coord_list, seqid in zip(kmer_coord_list_by_sequence, sequence_table.keys()):
            for kmer, coord in kmer_coord_list.items():
                all_kmer_list[kmer][seqid] = coord
    
    out_dir_list = [out_dir for i in range(len(kmer_list))]
    format_list = [format for i in range(len(kmer_list))]
    seqid_list = [seqid for seqid in sequence_table.keys()]
    seqid_list_list = [seqid_list for i in range(len(kmer_list))]
    with concurrent.futures.ProcessPoolExecutor() as executer:
        executer.map(write_kmer_table, out_dir_list, kmer_list, all_kmer_list.values(), seqid_list_list, format_list)


def get_kmer_coord_list(source_coord_list, start, end, kmer_size):
    list_start = bisect.bisect_left(source_coord_list, start)
    list_end = bisect.bisect_right(source_coord_list, end - kmer_size + 1)
    result_list = source_coord_list[list_start:list_end]
    return result_list


def get_kmer_coord_list_revcomp(source_coord_list, start, end, kmer_size):
    list_start = bisect.bisect_left(source_coord_list, start)
    list_end = bisect.bisect_right(source_coord_list, end - kmer_size + 1)
    tmp_list = source_coord_list[list_start:list_end]
    result_list = [c + kmer_size - 1 for c in tmp_list]
    return result_list


def output_kmer_coord_on_feature(kmer_table_dir, feature_list, kmer_size=6, kmer='AAAAAA', target_type='gene'):
    if 'N' in kmer:
        return output_kmer_coord_on_feature_with_gap(kmer_table_dir, feature_list, kmer_size, kmer, target_type)
    if __name__ == "__main__":
        signal(SIGPIPE, SIG_DFL) # Supress broken pipe error
    
    target_kmer_coord_table = rw_table.read_kmer_table(kmer_table_dir, kmer)
    target_kmer_revcomp_coord_table = rw_table.read_kmer_table(kmer_table_dir, reverse_complement(kmer))
    
    for feature in feature_list:
        type = feature['type']
        if target_type and type != target_type:
            continue
        seqid = feature['seqid']
        source = feature['source']
        start = feature['start']
        end = feature['end']
        score = feature['score']
        strand = feature['strand']
        phase = feature['phase']
        id = feature['id']
        parent = feature['parent']
        
        if seqid not in target_kmer_coord_table.keys():
            if verbose:
                print(f'Warning: {seqid} not found in table.', file=sys.stderr)
                print('Skipping ...', seqid, source, type, start, end, score, strand, phase, f'{f"ID={id}" if id else f"Parent={parent}"}', sep='\t', file=sys.stderr)
            continue

        if strand == '-':
            revcomp_coord_list = target_kmer_revcomp_coord_table[seqid]
            #result_list = [c + kmer_size - 1 for c in revcomp_coord_list if start <= c and c  <= end - kmer_size + 1]
            result_list = get_kmer_coord_list_revcomp(revcomp_coord_list, start, end, kmer_size)
        else:
            coord_list = target_kmer_coord_table[seqid]
            #result_list = [c for c in coord_list if start <= c and c  <= end - kmer_size + 1]
            result_list = get_kmer_coord_list(coord_list, start, end, kmer_size)
        if len(result_list) == 0:
            continue
        feature['kmer_'+kmer] = result_list
        kmer_pos_list = ','.join(map(str, result_list))
        if __name__ == "__main__":
            print(seqid, source, type, start, end, score, strand, phase, f'{f"ID={id}" if id else f"Parent={parent}"};pos_{kmer}={kmer_pos_list}', sep='\t')
    
    return feature_list


def create_kmer_combination_with_gap(kmer='AANAAA'):
    gap_char = 'N'
    kmer_list = ['']
    for base in kmer:
        if base == gap_char:
            kmer_list = [kmer + 'A' for kmer in kmer_list] + [kmer + 'T' for kmer in kmer_list] + [kmer + 'G' for kmer in kmer_list] + [kmer + 'C' for kmer in kmer_list]
        else:
            kmer_list = [kmer + base for kmer in kmer_list]
    return kmer_list


def output_kmer_coord_on_feature_with_gap(kmer_table_dir, feature_list, kmer_size=6, kmer='AANAAA', target_type='gene'):
    if __name__ == "__main__":
        signal(SIGPIPE, SIG_DFL) # Supress broken pipe error
    
    target_kmer_list = create_kmer_combination_with_gap(kmer)
    target_kmer_coord_table_list = [rw_table.read_kmer_table(kmer_table_dir, tmp_kmer) for tmp_kmer in target_kmer_list]
    target_kmer_revcomp_coord_table_list = [rw_table.read_kmer_table(kmer_table_dir, reverse_complement(tmp_kmer)) for tmp_kmer in target_kmer_list]
    
    for feature in feature_list:
        type = feature['type']
        if target_type and type != target_type:
            continue
        seqid = feature['seqid']
        source = feature['source']
        start = feature['start']
        end = feature['end']
        score = feature['score']
        strand = feature['strand']
        phase = feature['phase']
        id = feature['id']
        parent = feature['parent']
        
        result_all_list = []
        for target_kmer_coord_table, target_kmer_revcomp_coord_table in zip(target_kmer_coord_table_list, target_kmer_revcomp_coord_table_list):
            if seqid not in target_kmer_coord_table.keys():
                if verbose:
                    print(f'Warning: {seqid} not found in table.', file=sys.stderr)
                    print('Skipping ...', seqid, source, type, start, end, score, strand, phase, f'{f"ID={id}" if id else f"Parent={parent}"}', sep='\t', file=sys.stderr)
                continue

            if strand == '-':
                revcomp_coord_list = target_kmer_revcomp_coord_table[seqid]
                #result_list = [c + kmer_size - 1 for c in revcomp_coord_list if start <= c and c  <= end - kmer_size + 1]
                result_list = get_kmer_coord_list_revcomp(revcomp_coord_list, start, end, kmer_size)
            else:
                coord_list = target_kmer_coord_table[seqid]
                #result_list = [c for c in coord_list if start <= c and c  <= end - kmer_size + 1]
                result_list = get_kmer_coord_list(coord_list, start, end, kmer_size)
            
            result_all_list += result_list

        if len(result_all_list) == 0:
            continue
        result_all_list = sorted(list(set(result_all_list)))
        feature['kmer_'+kmer] = result_all_list
        kmer_pos_list = ','.join(map(str, result_all_list))

        if __name__ == "__main__":
            print(seqid, source, type, start, end, score, strand, phase, f'{f"ID={id}" if id else f"Parent={parent}"};pos_{kmer}={kmer_pos_list}', sep='\t')
    
    return feature_list


def count_kmer(kmer_table_dir, strand='both', graph_flag=False, yticks_space_size=None):
    file_list = [os.path.basename(file_path) for file_path in glob.glob(kmer_table_dir+'/*')]
    exist_kmer_list = []
    for file in file_list:
        match = re.search(r'([ATGC]{1,}).(csv|pkl).gz$', file)
        if match:
            exist_kmer_list.append(match[1])
    exist_kmer_list = sorted(list(set(exist_kmer_list)))
    if len(exist_kmer_list) == 0:
        print('Error: table not found.', file=sys.stderr)
        return
    
    cnt_list = []
    cnt_table = {}
    print('Mer\tCount')    
    kmer_table_dir_list = [kmer_table_dir for i in range(len(exist_kmer_list))]
    strand_list = [strand for i in range(len(exist_kmer_list))]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for kmer, cnt in zip(exist_kmer_list, executor.map(task_count_kmer, kmer_table_dir_list, exist_kmer_list, strand_list)):
            cnt_list.append(cnt)
            cnt_table[kmer] = cnt
            print(f'{kmer}\t{cnt}')

    if(graph_flag):
        import matplotlib.pyplot as plt
        import numpy as np

        fig = plt.figure(figsize=(len(exist_kmer_list)/ 4 / 6, 10), tight_layout=True)
        strand_to_symbol = {'both': '+/-', 'plus': '+', 'minus': '-'}
        fig.suptitle(f"k-mer count (strand {strand_to_symbol[strand]})")
        #fig.supxlabel("k-mer")
        #fig.supylabel("count")

        sort_base = 'ACGT'
        max_cnt = max(cnt_list)
        if(yticks_space_size):
            ticks_num = max_cnt // yticks_space_size + 1
            tmp_max =  ticks_num * yticks_space_size
            yticks = np.asarray(range(0, tmp_max+1, yticks_space_size), dtype=int)
        else:
            yticks = np.asarray(np.linspace(0, max_cnt if max_cnt % 2 == 0 else max_cnt + 1, 4), dtype=int)
        
        for i in range(len(sort_base)):
            sub_kmer_list = [kmer for kmer in exist_kmer_list if kmer[0] == sort_base[i]]
            sub_cnt_list = [cnt for kmer, cnt in cnt_table.items() if kmer[0] == sort_base[i]]

            ax = fig.add_subplot(4, 1, i+1)
            ax.set_xmargin(0)
            ax.set_ymargin(0)
            ax.grid(axis='y', alpha=0.3)
            ax.bar(sub_kmer_list, sub_cnt_list, 0.5)
            ax.set_xticks(range(len(sub_kmer_list)))
            ax.set_xticklabels(sub_kmer_list, rotation='vertical', fontsize=10, fontname='monospace')
            ax.set_yticks(yticks)
            ax.ticklabel_format(style='plain',axis='y')
        
        fig.savefig('graph.png')

def task_count_kmer(kmer_table_dir, kmer, strand):
    cnt = 0
    if strand in ('both', 'plus'):
        kmer_table = rw_table.read_kmer_table(kmer_table_dir, kmer)
        for kmer_coord_list in kmer_table.values():
            cnt += len(kmer_coord_list)
    if strand in ('both', 'minus'):
        kmer_table_rev = rw_table.read_kmer_table(kmer_table_dir, reverse_complement(kmer))
        for kmer_coord_list in kmer_table_rev.values():
            cnt += len(kmer_coord_list)
    return cnt

def get_table_path_list(root_path):
    p = re.compile(r'.*/[ATGC]{1,}\.(csv|pkl)\.gz$')
    table_path_list = [dir for dir in glob.glob(f'{root_path}/**/', recursive=True) if any([ p.match(file) for file in glob.glob(dir+'/*')])]
    return table_path_list

def read_kmer_table_all(table_dir_path):
    p = re.compile(r'.*/([ATGC]{1,})\.(csv|pkl)\.gz$')
    kmer_file_list = glob.glob(table_dir_path+'/*.csv.gz') + glob.glob(table_dir_path+'/*.pkl.gz')
    kmer_exist_list = sorted(list(set([p.match(kmer_file)[1] for kmer_file in kmer_file_list])))
    
    table_dir_path_list = [table_dir_path for i in range(len(kmer_exist_list))]
    with concurrent.futures.ProcessPoolExecutor() as executer:
            kmer_table_list = executer.map(rw_table.read_kmer_table, table_dir_path_list, kmer_exist_list)
    return {kmer: kt for kmer, kt in zip(kmer_exist_list, kmer_table_list)}
    #return { kmer: read_kmer_table(table_dir_path, kmer) for kmer in kmer_exist_list }


def load_kmer_tables(tables_dir_root_path):
    table_path_list = get_table_path_list(tables_dir_root_path)
    

def get_gff_path_list(root_path):
    p = re.compile(r'.*\.(gff|gff3)(\.gz)?$')
    gff_file_list = [file_path for file_path in glob.glob(root_path+'/**/*', recursive=True) if p.match(file_path)]
    return gff_file_list


def output_kmer_coord_in_table(table_dir_path, target_kmer, seqid_keyword="*"):
    if 'N' in target_kmer:
        return output_kmer_coord_in_table_with_gap(table_dir_path, target_kmer, seqid_keyword)
    if __name__ == "__main__":
        signal(SIGPIPE, SIG_DFL) # Supress broken pipe error
    target_kmer_coord_table = rw_table.read_kmer_table(table_dir_path, target_kmer)
    escape_table = {'\\':'\\\\', '.':'\\.', '*':'.*', '+':'\\+', '?':'\\?', '{':'\\{', '}':'\\}', '(':'\\(', ')':'\\)', '[':'\\[', ']':'\\]', '^':'\\^', '$':'\\$', '|':'\\|'}
    seqid_keyword_re_str = seqid_keyword
    for befor, after in escape_table.items():
        seqid_keyword_re_str = seqid_keyword_re_str.replace(befor, after)
    seqid_keyword_re_str = '^' + seqid_keyword_re_str + '$'
    p = re.compile(seqid_keyword_re_str)
    picked_kmer_coord_table = {}
    for seqid, coord_list in target_kmer_coord_table.items():
        if len(coord_list) > 0 and p.match(seqid):
            picked_kmer_coord_table[seqid] = coord_list
            if __name__ == "__main__":
                print(seqid, ','.join(map(str,coord_list)), sep=',')
    return picked_kmer_coord_table


def output_kmer_coord_in_table_with_gap(table_dir_path, target_kmer, seqid_keyword="*"):
    if __name__ == "__main__":
        signal(SIGPIPE, SIG_DFL) # Supress broken pipe error
    
    target_kmer_list = create_kmer_combination_with_gap(target_kmer)
    target_kmer_coord_table_list = [rw_table.read_kmer_table(table_dir_path, tmp_kmer) for tmp_kmer in target_kmer_list]
    escape_table = {'\\':'\\\\', '.':'\\.', '*':'.*', '+':'\\+', '?':'\\?', '{':'\\{', '}':'\\}', '(':'\\(', ')':'\\)', '[':'\\[', ']':'\\]', '^':'\\^', '$':'\\$', '|':'\\|'}
    seqid_keyword_re_str = seqid_keyword
    for befor, after in escape_table.items():
        seqid_keyword_re_str = seqid_keyword_re_str.replace(befor, after)
    seqid_keyword_re_str = '^' + seqid_keyword_re_str + '$'
    p = re.compile(seqid_keyword_re_str)
    picked_kmer_coord_table = {}
    for seqid in target_kmer_coord_table_list[0].keys():
        coord_list = []
        for target_kmer_coord_table in target_kmer_coord_table_list:
            coord_list += target_kmer_coord_table[seqid]
        coord_list = sorted(set(coord_list))
        if len(coord_list) > 0 and p.match(seqid):
            picked_kmer_coord_table[seqid] = coord_list
            if __name__ == "__main__":
                print(seqid, ','.join(map(str,coord_list)), sep=',')
    return picked_kmer_coord_table


def command_create(args):
    generate_kmer_coord_table(args.fasta_file, args.kmer_size, args.table_dir, args.format)


def command_lookup(args):
    output_kmer_coord_on_feature(args.table_dir, read_gff(args.gff_file), len(args.target_kmer), args.target_kmer, args.type)

def command_stat(args):
    count_kmer(args.table_dir, args.strand, args.graph, args.yticks_space_size)

def command_pickup(args):
    output_kmer_coord_in_table(args.table_dir, args.target_kmer, args.seqid)

def main():
    parser = argparse.ArgumentParser(description='k-mer table')
    
    parser.add_argument('-v', '--verbose', action='store_true', help='output warning')

    subparsers = parser.add_subparsers()

    parser_create = subparsers.add_parser('create')
    parser_create.add_argument('table_dir')
    parser_create.add_argument('fasta_file')
    parser_create.add_argument('kmer_size', type=int)
    parser_create.add_argument('--format', choices=['csv', 'pkl'])
    parser_create.set_defaults(handler=command_create)

    parser_lookup = subparsers.add_parser('lookup')
    parser_lookup.add_argument('table_dir')
    parser_lookup.add_argument('gff_file')
    parser_lookup.add_argument('target_kmer')
    parser_lookup.add_argument('--type')
    parser_lookup.set_defaults(handler=command_lookup)

    parser_pickup = subparsers.add_parser('pickup')
    parser_pickup.add_argument('table_dir')
    parser_pickup.add_argument('target_kmer')
    parser_pickup.add_argument('--seqid', default="*")
    parser_pickup.set_defaults(handler=command_pickup)

    parser_stat = subparsers.add_parser('stat')
    parser_stat.add_argument('table_dir')
    parser_stat.add_argument('--strand', choices=['both', 'plus', 'minus'], default='both', help='default: both')
    parser_stat.add_argument('--graph', action='store_true', help="output bar graph option. filename=graph.png")
    parser_stat.add_argument('--yticks_space_size', default=None, type=int, help="specify graph yticks space size. e.g. 1000000")
    parser_stat.set_defaults(handler=command_stat)

    args = parser.parse_args()
    
    global verbose
    verbose = args.verbose
    if hasattr(args, 'handler'):
        args.handler(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()