import os, array, gzip, pickle, glob, re
import concurrent.futures
import itertools

store_kmer_table = {}
store_kmer_table_status = {}

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

def get_table_path_list(root_path):
    p = re.compile(r'.*/[ATGC]{1,}\.(csv|pkl)\.gz$')
    table_path_list = [dir for dir in glob.glob(f'{root_path}/**/', recursive=True) if any([ p.match(file) for file in glob.glob(dir+'/*')])]
    return table_path_list

def get_exist_kmer_file_path_list(table_dir_path):
    p = re.compile(r'.*/([ATGC]{1,})\.(csv|pkl)\.gz$')
    exist_kmer_list = [p.match(file)[1] for file in glob.glob(table_dir_path+'/*') if p.match(file)]
    return sorted(set(exist_kmer_list))

def is_exist_kmer_file_path(table_dir_path, kmer):
    if len(glob.glob(f'{table_dir_path}/{kmer}.*')) > 0:
        return True
    else:
        return False

def get_exist_kmer_file_max_length(table_dir_path):
    exist_kmer_list = get_exist_kmer_file_path_list(table_dir_path)
    max_length = len(max(exist_kmer_list, key=len))
    return max_length


def generate_longer_kmer_table(table_dir_path, kmer):
    sub_kmer_max_length = get_exist_kmer_file_max_length(table_dir_path)
    offset_and_kmer_list = []
    offset_pos = 0

    while True:
        q, mod = divmod(offset_pos+sub_kmer_max_length, len(kmer))
        if q == 1 and mod > 0:
            offset_pos -= mod
            continue
        sub_kmer = kmer[offset_pos:offset_pos+sub_kmer_max_length]
        offset_and_kmer_list.append( (offset_pos, sub_kmer) )
        if q == 1 and mod == 0:
            break
        offset_pos += sub_kmer_max_length
    
    target_kmer_coord_table = None
    for offset_and_kmer in offset_and_kmer_list:
        offset_pos, sub_kmer = offset_and_kmer
        target_sub_kmer_coord_table = read_kmer_table(table_dir_path, sub_kmer)
        if offset_pos == 0:
            target_kmer_coord_table = target_sub_kmer_coord_table
            continue

        for seqid, coord_list in target_kmer_coord_table.items():
            downstream_coord_list = set(target_sub_kmer_coord_table[seqid])
            filterd_coord_list = [coord for coord in coord_list if coord+offset_pos in downstream_coord_list]
            target_kmer_coord_table[seqid] = filterd_coord_list

    return target_kmer_coord_table


def read_kmer_table(table_dir_path, kmer):
    if table_dir_path in  store_kmer_table and kmer in store_kmer_table[table_dir_path]:
        return store_kmer_table[table_dir_path][kmer]
    if len(kmer) > get_exist_kmer_file_max_length(table_dir_path):
        return generate_longer_kmer_table(table_dir_path, kmer)
    
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


def load_kmer_table_all(table_dir_path):
    p = re.compile(r'.*/([ATGC]{1,})\.(csv|pkl)\.gz$')
    exist_kmer_list = [p.match(file)[1] for file in glob.glob(table_dir_path+'/*') if p.match(file)]
    exist_kmer_list = sorted(set(exist_kmer_list))
    global store_kmer_table
    store_kmer_table[table_dir_path] = {}
    global store_kmer_table_status
    store_kmer_table_status[table_dir_path] = "loading"
    for kmer in exist_kmer_list:
        store_kmer_table[table_dir_path][kmer] = read_kmer_table(table_dir_path, kmer)
    
    store_kmer_table_status[table_dir_path] = "loaded"
    return store_kmer_table


def load_kmer_table_all_async(table_dir_path):
    executor = concurrent.futures.ThreadPoolExecutor()
    executor.submit(load_kmer_table_all, table_dir_path)
    executor.shutdown(wait=False)

def check_status_table_loaded(table_dir_path):
    status = "not loaded"
    if table_dir_path in store_kmer_table_status:
        status = store_kmer_table_status[table_dir_path]
    return status


def generate_kmer_coord_table_async(fasta_file, kmer_size, out_dir, format):
    executor = concurrent.futures.ThreadPoolExecutor()
    executor.submit(generate_kmer_coord_table, fasta_file, kmer_size, out_dir, format)
    executor.shutdown(wait=False)


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
