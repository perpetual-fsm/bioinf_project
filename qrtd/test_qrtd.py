from qrtd import qrtd

def test_qrtd():
	taxa = 'A B C D E'
	a = '(A,C,((B,D),E));'
	b = '(C,(B,D),(A,E));'
	out = qrtd(taxa, a, b)
	assert 2 == out