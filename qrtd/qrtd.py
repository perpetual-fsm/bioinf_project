import numpy as np
from newick import loads

def k_subsets_i(n, k):
	if k == 0 or n < k:
		yield []
	elif n == k:
		yield list(range(n))
	else:
		# Use recursive formula based on binomial coeffecients:
		# choose(n, k) = choose(n - 1, k - 1) + choose(n - 1, k)
		for s in k_subsets_i(n - 1, k - 1):
			s.append(n - 1)
			yield s
		for s in k_subsets_i(n - 1, k):
			yield s

def gen_triplets(n):
	return k_subsets_i(n, 3)


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

class Tree:
	def __init__(self, newick):
		self.string = newick
		self.parse_newick(newick)
		# self.parents = np.unique([node.no for node in self.tree])
		self.ntaxa.sort()
		# print('<<< Print Tree')
		# for key, node in self.tree.items():
		# 	print('node : {}'.format(key))
		# 	try:
		# 		print('parent: {}'.format(node.parent.no))
		# 	except AttributeError:
		# 		print('parent: None')
		# 	t = [x.no for x in node.children if x]
		# 	print('children: {}'.format(t))
		# 	try:
		# 		t = self.labels[node.no]
		# 	except KeyError:
		# 		t = 'None'
		# 	print('Label: {}'.format(t))
		# print('>>> End Tree')

	def create_tree(self, node, p_id):
		for child in node.descendants:
			if not child.descendants:
				tmp = len(self.tree)
				self.tree[tmp] = Node(tmp, self.tree[p_id])
				self.tree[p_id].add_child(self.tree[tmp])
				self.label_type[child.name] = tmp
				self.labels[tmp] = child.name
				self.ntaxa.append(child.name)
			else:
				tmp = len(self.tree)
				self.tree[tmp] = Node(tmp, self.tree[p_id])
				self.tree[p_id].add_child(self.tree[tmp])
				self.create_tree(child, tmp)



	def parse_newick(self, newick):

		tree = loads(newick)
		self.tree = {}
		self.label_type = {}
		self.labels = {}
		self.ntaxa = []
		self.tree[0] = Node(0, None)
		self.create_tree(tree[0], 0)

		# index = 0

		# symbols = ['(', ')', ':', ';', ',']

		# self.tree = {}
		# self.tree[0] = Node(0, None) # artificial root
		# self.labels = {}
		# self.weights = {}
		# self.label_type = {}
		# self.ntaxa = []
		# parent = -1
		# current = 0
		# save = 0
		# prev_symbol = None
		# while index < len(newick):
		# 	if newick[index] == '(':
		# 		# new node
		# 		parent = current
		# 		current = len(self.tree)
		# 		self.tree[current] =  Node(current, self.tree[parent])
		# 		self.tree[parent].add_child(self.tree[current])
		# 		prev_symbol = '('
		# 		index += 1
		# 		continue
		# 		# done open
		# 	elif newick[index] == ')':
		# 		if prev_symbol == ',':
		# 			temp = len(self.tree)
		# 			self.tree[temp] = Node(temp, self.tree[current])
		# 			self.tree[current].add_child(self.tree[temp])
		# 			temp = None
		# 		save = current
		# 		current = parent
		# 		parent = self.tree[current].parent.no
		# 		index += 1
		# 		prev_symbol = ')'
		# 		continue
		# 		# done close
		# 	elif newick[index] == ':':
		# 		if prev_symbol in ['(', ',']:
		# 			temp = len(self.tree)
		# 			self.tree[temp] = Node(temp, self.tree[current])
		# 			self.tree[current].add_child(self.tree[temp])
		# 			temp = None
		# 		index += 1
		# 		prev_symbol = ':'
		# 		continue
		# 	elif newick[index] == ';':
		# 		index = len(newick)
		# 		break
		# 		# done - break
		# 	elif newick[index] == ',':
		# 		if prev_symbol == ',':
		# 			temp = len(self.tree)
		# 			self.tree[temp] = Node(temp, self.tree[current])
		# 			temp = None
		# 		index += 1
		# 		prev_symbol = ','
		# 		# done comma
		# 		continue
		# 	else:
		# 		taxa = ''
		# 		while index < len(newick) and newick[index] not in symbols:
		# 			taxa += newick[index]
		# 			index += 1
		# 		if prev_symbol == ')':
		# 			self.labels[taxa] = save
		# 			self.label_type[save] = taxa
		# 		elif prev_symbol == ':':
		# 			self.weights[save] = int(taxa)
		# 		else:
		# 			temp = len(self.tree)
		# 			self.labels[taxa] = temp
		# 			self.label_type[temp] = taxa
		# 			self.tree[temp] = Node(temp, self.tree[current])
		# 			self.tree[current].add_child(self.tree[temp])
		# 			save = temp
		# 			self.ntaxa.append(taxa)
		# 			temp = None
		# 		index += 1
		# 		prev_symbol = 'L'
		# 		continue
		# 		# done label

	def path_to_root(self, taxa):
		tree=self.tree
		start_node=self.label_type[taxa]
		path = [tree[start_node].no]
		parent = tree[start_node].parent.no
		while parent != 0:
			path.append(parent)
			parent = tree[parent].parent.no
		path.append(parent)
		return path

	def find_path(self, x, y):
		path_x = self.path_to_root(x)
		path_y = self.path_to_root(y)

		rev_path_x = list(reversed(path_x))
		rev_path_y = list(reversed(path_y))

		while len(rev_path_y) != len(rev_path_x):
			if len(rev_path_y) < len(rev_path_x):
				rev_path_y += [-1]
			else:
				rev_path_x += [-1]

		assert rev_path_y[0] == rev_path_x[0]

		prev = rev_path_x[0]
		for ind, tupla in enumerate(zip(rev_path_x, rev_path_y)):
			a, b = tupla
			if a == b:
				prev = a
				continue
			else:
				pathA = [x for x in rev_path_x[ind:] if x != -1] 
				pathB = [x for x in rev_path_y[ind:] if x != -1] 
				whole_path = list(reversed(pathA)) + [prev] +pathB
				break
		return whole_path

def find_center_from_three_paths(t, path1, path2, path3):

	internals = {}
	for i in range(len(t.tree)):
		internals[i] = 0

	for path in (path1, path2, path3):
		for node in path:
			internals[node] += 1
			if internals[node] == 3:
				center = node
				break
	prevA = path1[0]
	for node in path1[1:]:
		if node == center:
			break
		else:
			prevA = node
	prevB = path2[0]
	for node in path2[1:]:
		if node == center:
			break
		else:
			prevB = node
	prevC = path3[0]
	for node in path3[1:]:
		if node == center:
			break
		else:
			prevC = node

	# print('prevA: {}, prevB: {}, prevC: {}'.format(prevA, prevB, prevC))

	return center, prevA, prevB, prevC


def find_center(tree, a,b,c):
	pathAB = tree.find_path(a, b)
	pathBC = tree.find_path(b, c)
	pathCA = tree.find_path(c, a)
	center, prevA, prevB, prevC = find_center_from_three_paths(tree, pathAB, pathBC, pathCA)

	parent = None
	try:
		if prevA == tree.tree[center].parent.no:
			parent = 0
	except AttributeError:
		pass
	try:
		if prevB == tree.tree[center].parent.no:
			parent = 1
	except AttributeError:
		pass
	try:
		if prevC == tree.tree[center].parent.no:
			parent = 2
	except AttributeError:
		pass
	return center, prevA, prevB, prevC, parent

def get_leaves(tree, forbidden, start_node, center):
	if tree.tree[start_node].children and start_node not in [forbidden, center]:
		result = []
		for node in tree.tree[start_node].children:
			ans = get_leaves(tree, start_node, node.no, center)
			if ans:
				try:
					result.extend(ans)
				except TypeError:
					result.append(ans)
		return result
	elif start_node in [forbidden, center]:
		return None
	else:
		return tree.tree[start_node].no


def get_leaves_from_subtree(tree, forbidden, start_node, from_child=False):
	if not from_child:
		leaves = []
		if tree.tree[start_node].children:
			tnodes = tree.tree[start_node].children + [tree.tree[start_node].parent]
			nodes = [node.no for node in tnodes]
			for node in nodes:
				if node != forbidden:
					rec = get_leaves(tree, forbidden, node, forbidden)
					if rec:
						try:
							leaves.extend(rec)
						except TypeError:
							leaves.append(rec)
				else:
					continue
			return leaves
		else:
			return [tree.tree[start_node].no]
	else:
		leaves = []
		root = tree.tree[0].no
		tnodes = tree.tree[root].children
		nodes = [node.no for node in tnodes]
		for node in nodes:
			if node != forbidden:
				rec = get_leaves(tree, forbidden, node, forbidden)
				if rec:
					try:
						leaves.extend(rec)
					except TypeError:
						leaves.append(rec)
			else:
				continue
		return leaves

def count_diff_topologies(tree_a, tree_b, labels, a,b,c):

	# print()
	# print()
	# print('Start counting diff')
	# print('Chosen leaves')
	# print(a,b,c)

	nleaves = len(labels)
	label_dict = {}

	for ind, label in enumerate(labels):
		label_dict[label] = ind

	center_a = find_center(tree_a, a,b,c)
	center_b = find_center(tree_b, a,b,c)

	prev_nodes_a = center_a[1:4]
	prev_parent_a = center_a[4]

	leaves_a_subtree_a = get_leaves_from_subtree(tree_a, center_a[0], prev_nodes_a[0], from_child=prev_parent_a==0)
	leaves_a_subtree_b = get_leaves_from_subtree(tree_a, center_a[0], prev_nodes_a[1], from_child=prev_parent_a==1)
	leaves_a_subtree_c = get_leaves_from_subtree(tree_a, center_a[0], prev_nodes_a[2], from_child=prev_parent_a==2)
	# print('a subtree_a')
	# print(leaves_a_subtree_a)
	# print('a subtree_b')
	# print(leaves_a_subtree_b)
	# print('a subtree_c')
	# print(leaves_a_subtree_c)

	prev_nodes_b = center_b[1:4]
	prev_parent_b = center_b[4]

	leaves_b_subtree_a = get_leaves_from_subtree(tree_b, center_b[0], prev_nodes_b[0], from_child=prev_parent_b==0)
	leaves_b_subtree_b = get_leaves_from_subtree(tree_b, center_b[0], prev_nodes_b[1], from_child=prev_parent_b==1)
	leaves_b_subtree_c = get_leaves_from_subtree(tree_b, center_b[0], prev_nodes_b[2], from_child=prev_parent_b==2)

	# print('b subtree_a')
	# print(leaves_b_subtree_a)
	# print('b subtree_b')
	# print(leaves_b_subtree_b)
	# print('b subtree_c')
	# print(leaves_b_subtree_c)

	topology_a = [3] * nleaves
	topology_b = [3] * nleaves

	for x in leaves_a_subtree_a:
		t_taxa = tree_a.labels[x]
		t_key = label_dict[t_taxa]
		topology_a[t_key] = 0
	for x in leaves_a_subtree_b:
		t_taxa = tree_a.labels[x]
		t_key = label_dict[t_taxa]
		topology_a[t_key] = 1
	for x in leaves_a_subtree_c:
		t_taxa = tree_a.labels[x]
		t_key = label_dict[t_taxa]
		topology_a[t_key] = 2

	# print(leaves_b_subtree_c)

	for x in leaves_b_subtree_a:
		t_taxa = tree_b.labels[x]
		t_key = label_dict[t_taxa]
		topology_b[t_key] = 0
	for x in leaves_b_subtree_b:
		t_taxa = tree_b.labels[x]
		t_key = label_dict[t_taxa]
		topology_b[t_key] = 1
	for x in leaves_b_subtree_c:
		t_taxa = tree_b.labels[x]
		t_key = label_dict[t_taxa]
		topology_b[t_key] = 2

	diff = 0
	for i in range(nleaves):
		# print('topology a : topology b')
		# print('{} : {}'.format(topology_a[i], topology_b[i]))
		if topology_a[i] != topology_b[i]:
			diff += 1

	# print('Stop counting diff')
	# print()
	# print()
	return diff

def qrtd(taxa, a, b):
	labels = taxa.split()
	nleaves = len(labels)

	tree_a = Tree(a)
	tree_b = Tree(b)

	triplets = gen_triplets(nleaves)
	diff = 0
	for triplet in triplets:
		taxa_triplet = labels[triplet[0]], labels[triplet[1]], labels[triplet[2]]
		diff += count_diff_topologies(tree_a, tree_b, labels, *taxa_triplet)
	return int(diff / 4)

if __name__ == '__main__':
	taxa = 'A B C D E'
	a = '(A,C,((B,D),E));'
	b = '(C,(B,D),(A,E));'
	# a = '(A,B,(C,D),E);'
	# b = '(A,C,(B,D),E);'

	# file = 'input.txt'
	# with open(file) as tf:
	# 	t_input = tf.readlines()
	# 	taxa = t_input[0].strip()
	# 	a = t_input[1].strip()
	# 	b = t_input[2].strip()
	# print(taxa)
	# print(a)
	# print(b)
	print(qrtd(taxa, a, b))	