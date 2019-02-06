from mgap import maximal_gap

def test_smgb():
	seqa = 'AACGTA'
	seqb = 'ACACCTA'
	score, score_matrix, path_matrix = maximal_gap(seqa, seqb, -1)
	assert score == 3