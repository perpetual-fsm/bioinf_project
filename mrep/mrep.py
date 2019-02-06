import argparse
import itertools

class Node:
	def __init__(self, no, parent):
		self.parent = parent
		self.no = no
		self.children = []

	def add_child(self, child):
		self.children.append(child)

	def remove_child(self, child):
		self.children.remove(child)

	def update_parent(self, parent):
		self.parent = parent

class SuffixTree:
	def __init__(self, word=""):
		self.nodes = [Node(0, None)]
		self.edges = {}
		self.children = {}
		if type(word) == str:
			self.add_word(word)


	def find_suffix_parent(self, substring_index, root):
		for child in root.children:
			source_edge, target_edge = self.edges[(root.no, child.no)]
			if self.word[substring_index : substring_index - source_edge + target_edge] == self.word[source_edge : target_edge]:
				return self.find_suffix_parent(substring_index - source_edge + target_edge, child)
			elif self.word[substring_index] == self.word[source_edge]:
				return child, substring_index, True
		return root, substring_index, False

	def add_node(self, parent, source_edge, target_edge, child=None):
		if child is None:
			child = Node(len(self.nodes), parent)
		self.nodes.append(child)
		parent.add_child(child)
		self.edges[(parent.no, child.no)] = [source_edge, target_edge]


	def add_word(self, word):
		if word[-1] != '$':
			word += '$'
		self.word = word
		self.length = len(word)
		for i in range(self.length):
			parent, source_edge, overlap = self.find_suffix_parent(i, self.nodes[0])
			if overlap:
				p_edge_source, p_edge_target = self.edges[(parent.parent.no, parent.no)]
				insert = 0	
				while word[source_edge:source_edge + insert] == word[p_edge_source:p_edge_source + insert]:
					insert += 1
				new_node = Node(len(self.nodes), parent.parent)
				new_node.add_child(parent)
				self.add_node(parent.parent, p_edge_source, p_edge_source + insert - 1, new_node)
				del self.edges[(parent.parent.no, parent.no)]
				parent.parent.remove_child(parent)
				parent.update_parent(new_node)
				self.edges[(parent.parent.no, parent.no)] = [p_edge_source + insert - 1, p_edge_target]
				self.add_node(new_node, source_edge + insert - 1, self.length)
			else:
				self.add_node(parent, source_edge, self.length)

	def get_word(self, leaf):
		word = ''
		while leaf.no != 0:
			parent_index = self.edges[(leaf.parent.no, leaf.no)]
			word = self.word[parent_index[0]:parent_index[1]]+ word
			leaf = leaf.parent
		return word.rstrip("$")

	def get_children_no(self, node):
		if node not in self.children:
			self.children[node] = len(node.children) + sum([self.get_children_no(child) for child in node.children])
		return self.children[node]

def mrep(dna_string, disable_len=False):
	repeats = {}
	dna_suffix_tree = SuffixTree(dna_string)
	for node in dna_suffix_tree.nodes[1:]:
		children_no = dna_suffix_tree.get_children_no(node)
		node_word = dna_suffix_tree.get_word(node)
		if children_no >= 2 and (len(node_word) >= 20 or disable_len):
			if children_no not in repeats:
				repeats[children_no] = [node_word]
			else:
				repeats[children_no].append(node_word)
	max_repeats = []
	for value in repeats.values():
		if len(value) == 1:
			max_repeats.append(value)
		else:
			t_max = []
			for ta_word in value:
				t_value = []
				for tb_word in value:
					if ta_word != tb_word:
						t_value.append(ta_word not in tb_word)
				if all(t_value) and t_value:
					t_max.append(ta_word)

			max_repeats.append(t_max)
	final_result = list(itertools.chain(*max_repeats))
	return final_result

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Arguments for Identifying Maximal Repeats')
	parser.add_argument(
		'-file',
		dest='file',
		type=str,
		nargs=1,
		default=None,
		help='Specify filename with input'
	)
	parser.add_argument(
		'-input'
		'-nr', '--no_recov',
		dest='input',
		default=None,
		type=str,
		nargs=1,
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
		seq = 'TAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTATTATATAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTAT'
	if args.file:
		with open(args.file[0]) as tf:
			seq = tf.readlines()[0].strip()
	if args.input:
		seq = args.input[0]
	if args.verbose:
		print('Sekwencja:')
		print(seq)
		print('Wynik:')
	print('\n'.join(mrep(seq)))
