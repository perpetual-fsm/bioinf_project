import numpy as np
import argparse

def maximal_gap(seq_a, seq_b, gap_penalty=-1):
    rows_len = len(seq_b) + 1
    col_len = len(seq_a) + 1

    score_matrix = np.zeros((rows_len, col_len))
    path_matrix  = np.zeros((rows_len, col_len))
    for i in range(1, rows_len):
        for j in range(1, col_len):
            if seq_a[j-1] == seq_b[i-1]:
                match = score_matrix[i-1][j-1] + 1
                score_matrix[i][j] = match
                path_matrix[i][j] = 0 # D
            else:
                delete = score_matrix[i-1][j]
                insert = score_matrix[i][j-1]
                t_max = max(delete, insert)
                score_matrix[i][j] = t_max
                if score_matrix[i][j] == insert:
                    path_matrix[i][j] = 2 # L
                elif score_matrix[i][j] == match:
                    path_matrix[i][j] = 0 # D
            
    return int(len(seq_a) + len(seq_b) - 2* score_matrix[-1][-1]), score_matrix, path_matrix

def print_one_optimal(path_mat, seq_a, seq_b):
    str_a = seq_a
    str_b = seq_b
    insert_gap = lambda word, i: word[:i] + '-' + word[i:]
    i = len(seq_b)
    j = len(seq_a)
    while i > 0 and j > 0:
        if path_mat[i][j] == 1:
            i -= 1
            str_a = insert_gap(str_a, j)
        elif path_mat[i][j] == 2:
            j -= 1
            str_b = insert_gap(str_b, i)
        else:
            i -= 1
            j -= 1
    print(str_a)
    print(str_b)


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
    fasta_file.close()
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
    )
    parser.add_argument(
        '-verbose',
        action='store_true',
        dest='verbose',
        default=False,
    )
    args = parser.parse_args()
    if args.example:
        seqa = 'AACGTA'
        seqb = 'ACACCTA'
    if args.file:
        seqa, seqb = parse_input(args.file[0])
    if args.input:
        seqa, seqb = args.input
    score, score_matrix, path_matrix = maximal_gap(seqa, seqb, -1)
    if args.verbose:
        print('Sekwencje')
        print(seqa)
        print(seqb)
        print('Macierz wyników')
        print(score_matrix)
        print('Macierz ścieżek')
        print(path_matrix)
        print('Wynik')
    print(score)