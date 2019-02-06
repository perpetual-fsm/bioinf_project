import numpy as np
import argparse

def get_subs_score(char_a, char_b):
    if char_a == char_b:
        return 1
    return -1

def semiglobal_align(seq_a, seq_b, gap_penalty=-1):
    rows_len = len(seq_b) + 1
    col_len = len(seq_a) + 1

    score_matrix = np.zeros((rows_len, col_len))
    path_matrix  = np.zeros((rows_len, col_len))
    for i in range(1, rows_len):
        for j in range(1, col_len):
            match = score_matrix[i-1][j-1] + get_subs_score(seq_a[j-1], seq_b[i-1])
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            t_max = max(match, delete, insert)
            score_matrix[i][j] = t_max
            if t_max == match:
                path_matrix[i][j] = 0 # D
            elif t_max == delete:
                path_matrix[i][j] = 1 # U
            elif t_max == insert:
                path_matrix[i][j] = 2 # L

    max_row = max(range(col_len), key=lambda x: score_matrix[rows_len-1][x])
    max_col = max(range(rows_len), key=lambda x: score_matrix[x][col_len-1])
    if score_matrix[rows_len-1][max_row] >= score_matrix[max_col][col_len-1]:
        i = rows_len-1
        j = max_row
    else:
        i = max_col
        j = col_len - 1
    max_score = int(score_matrix[i][j])
    return score_matrix, path_matrix, i, j, max_score

def print_one_optimal(path_mat, seq_a, seq_b, row, col):
    str_a = seq_a
    str_b = seq_b
    for _ in range(len(seq_b) - row):
        str_a += '-'
    for _ in range(len(seq_a) - col):
        str_b += '-'
    insert_gap = lambda word, i: word[:i] + '-' + word[i:]
    i = row
    j = col
    while i > 1 and j > 1:
        if path_mat[i][j] == 1:
            i -= 1
            str_a = insert_gap(str_a, j)
        elif path_mat[i][j] == 2:
            j -= 1
            str_b = insert_gap(str_b, i)
        else:
            i -= 1
            j -= 1
    for _ in range(i-1):
        str_a = insert_gap(str_a, 0)
    for _ in range(j-1):
        str_b = insert_gap(str_b, 0)
    print(str_a)
    print(str_b)
    return str_a, str_b

def parse_input(t_file):
    fasta_file = open(t_file)
    t_seq = ''
    sequences = []
    for line in fasta_file.readlines():
        if line[0] == '>':
            if t_seq != '':
                sequences.append(t_seq.replace('\n', '').strip())
                t_seq = ''
        else:
            t_seq += line
    sequences.append(t_seq.replace('\n', '').strip())
    return sequences

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Arguments for semiglobal alignment')
    parser.add_argument(
        '-file',
        dest='file',
        type=str,
        nargs=1,
        default=None,
        help='Specify filename with input'
    )
    parser.add_argument(
        '-input',
        dest='input',
        default=None,
        type=str,
        nargs=2,
        help='Get input via cmd args.'
    )
    parser.add_argument(
        '-example',
        action='store_true',
        dest='example',
        default=False,
        help='Run semiglobal alignment with basic example'
    )
    parser.add_argument(
        '-verbose',
        action='store_true',
        dest='verbose',
        default=False,
        help='Be verbose'
    )
    args = parser.parse_args()
    if args.example:
        seqa = 'CAGCACTTGGATTCTCGG'
        seqb = 'CAGCGTGG'
    if args.file:
        seqa, seqb = parse_input(args.file)
    if args.input:
        seqa, seqb = args.input
    score_matrix, path_matrix, row, col, max_score = semiglobal_align(seqa, seqb, -1)
    if args.verbose:
        print('Sekwencje:')
        print(seqa)
        print(seqb)
        print('Macierz wyników:')
        print(score_matrix)
        print('Macierz ścieżek:')
        print(path_matrix)
        print('Lokalizacja najlepszego wyniku: {}-{}'.format(row,col))
        print('Wartość najlepszego wyniku: {}'.format(int(max_score)))
    print(max_score)
    print_one_optimal(path_matrix, seqa, seqb, row, col)