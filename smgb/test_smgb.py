from smgb import semiglobal_align, print_one_optimal

def test_smgb():
	seqa = 'CAGCACTTGGATTCTCGG'
	seqb = 'CAGCGTGG'
	score_matrix, path_matrix, row, col, max_score = semiglobal_align(seqa, seqb, -1)
	align_a, align_b = print_one_optimal(path_matrix, seqa, seqb, row, col)
	assert max_score == 4
	assert align_a == 'CAGCA-CTTGGATTCTCGG'
	assert align_b == '---CAGCGTGG--------'
