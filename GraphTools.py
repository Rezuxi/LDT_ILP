import networkx as nx
import asymmetree.tools as tools
import asymmetree.cograph as cg
from asymmetree.tools.GraphTools import disturb_graph, symmetric_diff
import random
import copy
import itertools


#####################################################################################
#																					#
#						Induced paths on 3 vertices (P3)							#
#																					#
#####################################################################################


def find_all_P3(G, get_triples=False, colored_G=None):
	'''
		Finds all connected triples that don't form a triangle
		if triples = True (the graph needs to be colored), we also check that the colors are distinct for each node in a triple (a, b, c)

		returns the P3s as a list of lists
		returns the triples as a list of tuples
	'''
	leaves = set()
	triples = []

	if colored_G:
		nodes = colored_G.nodes()	
	else:
		nodes = G.nodes()
	
	for e in G.edges():
		u_ID = e[0]
		v_ID = e[1]
		
		u_adj = G.adj[u_ID]
		v_adj = G.adj[v_ID]


		for w_ID in u_adj:
			# if they're the same node, skip
			if w_ID == v_ID:
				#print("{} and {} is the same node".format(w_ID, v_ID))
				continue
			# if they form a triangle, skip
			if w_ID in v_adj:
				#print("Triangle detected for the nodes: {}".format((v_ID, u_ID, w_ID)))
				continue
			# if we want a set of triples
			if get_triples == True:
				# if the nodes' colors are not pairwise distinct, skip
				if nodes[w_ID]['color'] == nodes[v_ID]['color']:
					#print("Not pairwise distinct colors for: {}".format((v_ID, w_ID, u_ID)))
					continue
				else:
					# if the triples has been counted already
					if (w_ID, v_ID, u_ID) in triples or (v_ID, w_ID, u_ID) in triples:
						#print("This triple is already included: {}".format((v_ID, w_ID, u_ID)))
						continue
					else:
						#print("Adding the triple: {}".format((w_ID, v_ID, u_ID)))
						triples.append((w_ID, v_ID, u_ID))
						leaves.add(u_ID)
						leaves.add(v_ID)
						leaves.add(w_ID)
			# if we want all P3
			# all P3 will be of the form (a, b, c) or (c, b, a). i.e., the edges are (a,b) and (b, c)
			else:
				if [w_ID, u_ID, v_ID] in triples or [v_ID, u_ID, w_ID] in triples:
					#print("This P3 is already included: {}".format((v_ID, u_ID, w_ID)))
					continue
				else:
					#print("Adding the P3: {}".format((w_ID, u_ID, v_ID)))
					triples.append([w_ID, u_ID, v_ID])


		for w_ID in v_adj:
			if w_ID == u_ID:
				#print("{} and {} is the same node".format(w_ID, u_ID))
				continue
			if w_ID in u_adj:
				#print("Triangle detected for the nodes: {}".format((v_ID, u_ID, w_ID)))
				continue
			if get_triples == True:
				if nodes[w_ID]['color'] == nodes[u_ID]['color']:
					#print("This P3 is already included: {}".format((u_ID, w_ID, v_ID)))
					continue
				else:
					if (w_ID, u_ID, v_ID) in triples or (u_ID, w_ID, v_ID) in triples:
						#print("This triple is already included: {}".format((u_ID, w_ID, v_ID)))
						continue
					else:
						#print("Adding the triple: {}".format((w_ID, u_ID, v_ID)))
						triples.append((w_ID, u_ID, v_ID))
						leaves.add(u_ID)
						leaves.add(v_ID)
						leaves.add(w_ID)
			else:
				if [w_ID, v_ID, u_ID] in triples or [u_ID, v_ID, w_ID] in triples:
					#print("This P3 is already included: {}".format((u_ID, v_ID, w_ID)))
					continue
				else:
					#print("Adding the P3: {}".format((w_ID, v_ID, u_ID)))
					triples.append([w_ID, v_ID, u_ID])
	return triples, leaves


def P3_regions(l, a = 1):
	'''
		l is a list of P3s (lists)
		a is the minimum number of elements regions need to have in common for them to be
		considered the same.

		returns a list of sets (regions) and a list with amount of P3s per region.
		the idx of the amount list matches the idx of the list of regions.
	'''

	regions = []
	amounts = []
	#curr_key = 0
	curr_value = 1
	while len(l) > 0:
		head, *tail = l
		head = set(head)
		n = -1
		while len(head) > n:
			n = len(head)
			rest = []
			for t in tail:
				t = set(t)
				if len(head.intersection(t)) >= a:
					curr_value += 1
					head |= t
				else:
					rest.append(t)
			tail = rest
		regions.append(head)
		amounts.append(curr_value)
		#curr_key += 1
		curr_value = 1
		l = tail
	return regions, amounts


def P3_distance(lengths, p1, p2):
	'''
		returns the min distance between the path p1 and
		assumes the shortest path between all pairs has been calculated
		no need to run nx.all_pairs_shortest_path_length() here
	'''
	min_dist = float('inf')
	for a in p1:
		d1 = lengths[a][p2[0]]
		d2 = lengths[a][p2[1]]
		d3 = lengths[a][p2[2]]
		min_value = min((d1, d2, d3))
		if min_value < min_dist:
			min_dist = min_value
	return min_dist



def unique_combinations(*l):
	'''
		get all combinations from the lists in l of length n.
	'''
	for c in itertools.combinations(l, 2):
		for pair in itertools.product(*c):
			yield pair

def regions_distance(G, regions, lengths=None):
	'''
		regions is a list of sets containing nodes. each set is its own region.
		G is the graph containing the nodes in the regions.

		return
			dictionary with key mapping to the index of the region in the regions list, the value of which is
			a dictionary with keys mapping to the index of other regions the value of which is the min distance 
			between the regions.
	'''

	if not lengths:
		lengths = dict(nx.all_pairs_shortest_path_length(G))
	region_distances = {}
	k = len(regions)
	for i in range(k):
		region_distances[i] = dict()
		for j in range(i+1, k):
			min_dist = float('inf')

			combinations = unique_combinations(regions[i], regions[j])		
			for c in combinations:
				# check if keys exist
				if c[0] in lengths:
					if c[1] in lengths[c[0]]:
						if lengths[c[0]][c[1]] < min_dist:
							min_dist = lengths[c[0]][c[1]]
			region_distances[i][j] = min_dist
	return region_distances 


def get_P3_data(G, colored_G=None):
	if colored_G:
		P3s, _ = find_all_P3(G, get_triples=True, colored_G=colored_G)
	else:
		P3s, _ = find_all_P3(G, get_triples=True)

	print("P3s: \n{}".format(P3s))
	print("length of P3s: {}".format(len(P3s)))
	regions, amounts = P3_regions(P3s)
	regions_distances = regions_distance(G, regions)

	return regions, amounts, regions_distances


#####################################################################################
#																					#
#									Graph Editing Tools 							#
#																					#
#####################################################################################


def is_compatible(G, colored_G = None):
	if colored_G:
		triples, leaves = find_all_P3(G, get_triples=True, colored_G=colored_G)
	else:
		triples, leaves = find_all_P3(G, get_triples=True)

	# if triples , leaves are empty then it's an LDT graph (if cograph)
	'''
	if triples == []:
		return True
	'''
	B = tools.Build(triples, leaves)
	tree_triples = B.build_tree()

	return True if tree_triples else False
	

def is_cograph(G):
	cotree = cg.Cotree.cotree(G)
	if cotree:
		return True
	return False


def is_properly_colored(G, colored_G = None):
	'''
		Checks if a graph is properly colored. (no adjacent vertices have the same color)
	'''
	if colored_G:
		nodes = colored_G.nodes()
	else:
		nodes = G.nodes()

	for e in G.edges():
		if nodes[e[0]]['color'] == nodes[e[1]]['color']:
			return False
	return True



class InvestigateGraph:

	def __init__(self, G):
		'''
			G is an LDT graph
		'''
		self._G = G
		self._G_perturbed = None

		self._is_cograph = True
		self._is_compatible = True

		self._triplesEdit_to_LDT = False
		self._cographEdit_to_LDT = False

		'''
			Count the P3s and regions of the graphs that are/become LDT-graphs and compare them to
			non LDT-graphs. See if there is anything standing out.
		'''

		# count how many times the perturbed graph remains properly colored after cograph editing and also for triples editing when adding edges
		self._count_cographEdit_remain_properly_colored = 0

		# count how many times any of the edits results in an LDT-graph.
		self._count_cographEdit_to_LDT = 0
		self._count_triplesEdit_to_LDT = 0
		
		# count how many times a graph went from being consistent to inconsistent by doing cograph editing
		# and the same for triple editing		
		self._count_triplesEdit_broke_cograph = 0
		self._count_cographEdit_broke_consistency = 0

		# count how many times a graph went from broken consistency to fixed by doing cograph editing
		# and the same for triple editing
		self._count_triplesEdit_fixed_cograph = 0
		self._count_cographEdit_fixed_consistency = 0

		# count how many times a graph went from broken consistency to fixed by doing cograph editing
		# and the same for triple editing
		self._count_triplesEdit_remained_cograph = 0
		self._count_cographEdit_remained_consistent = 0


		# count how many times the disturbed graph is not a cograph and not consistent
		# so that we can compare the success rate of each heuristic.
		self._count_dG_not_cograph = 0
		self._count_dG_not_consistent = 0
		
		# count how many times the disturbed graph remains a cograph or consistent so that 
		# we can compare the frequency at which each heuristic breaks the other property of 
		# LDT-graphs. i.e. cograph editing breaking consistency and vice versa.
		self._count_dG_cograph = 0
		self._count_dG_consistent = 0

		self._count_dG_notCograph_notConsistent = 0
		self._count_dG_cograph_notConsistent = 0
		self._count_dG_notCograph_consistent = 0


		self._count_triplesEdit_success = 0

	def perturb_graph(self, i_rate = None, d_rate = None):
		'''
			perturbs a graph until it is not an LDT-graph
			by default random values for deletion/insertion rate
		'''
		if i_rate == None:
			i_rate = round(random.random(), 1)
		if d_rate == None:
			d_rate = round(random.random(), 1)
		if i_rate == 0.0 and d_rate == 0.0:
			i_rate = round(random.uniform(0.1, 1.0), 1)
			d_rate = round(random.uniform(0.1, 1.0), 1)
		self._G_perturbed = disturb_graph(self._G, insertion_prob=i_rate, deletion_prob=d_rate)
		if is_cograph(self._G_perturbed):
			self._is_cograph = True
		else:
			self._is_cograph = False
		if is_compatible(self._G_perturbed):
			self._is_compatible = True
		else:
			self._is_compatible = False
		# make sure the disturbed graph is not an LDT
		if self._is_cograph and self._is_compatible:
			#print("adding noise again!")
			self.perturb_graph()




	def compare_graphs(self):
		"""
			compare the P3s of the LDT graph and the perturbed graph
		"""
		pass





	#########################################################################################################
	#																										#
	#										EDITING HEURISTICS												#
	#																										#
	#########################################################################################################


	# only deletes edges with highest weight for every triple to be removed
	def triples_editing(self):
		copy_G = self._G_perturbed.copy()
		#copy_G = copy.deepcopy(self.disturbed_G)
		# copy_G is colored (?) so dont need to pass in a different colored graph here
		triples, leaves = find_all_P3(copy_G, get_triples=True)
		
		if len(triples) == 0:
			#print("The set of triples is empty!\n")
			return None, None
		# NOTE: If G has less than two nodes an error will occur in the BUILD alg, specifically in stoer_wagner alg.

		B = tools.Build(triples, leaves, mincut=True)
		tree_triples = B.build_tree()
		#print("Cut list: {}".format(B.cut_list))
		# set weights to 0
		for a, b, c in triples:
			copy_G[a][c]['weight'] = 0
			copy_G[b][c]['weight'] = 0

		# set weights to the edges based on how often they appear in the set of triples
		if len(B.cut_list) > 0:
			for a, b, c in triples:
				copy_G[a][c]['weight'] += 1
				copy_G[b][c]['weight'] += 1

			# decide which edge(s) to remove
			for a, b, c in triples:
				if (a, b) in B.cut_list or (b, a) in B.cut_list:
					'''
						check if both edges exist. If not then we dont need to do anything since
						1 of the edges has already been cut.
					'''
					if not (copy_G.has_edge(a, c) and copy_G.has_edge(b, c)):
						continue
					# check which edge of (a, c) and (b, c) has the highest weight
					if copy_G[a][c]['weight'] >= copy_G[b][c]['weight']:
						# remove edge (a, c)
						copy_G.remove_edge(a, c)
						#print("Removing edge ({}, {})".format(a, c))
					else:
						copy_G.remove_edge(b, c)
						#print("Removing edge ({}, {})".format(b, c))
		return copy_G, tree_triples


	def cograph_editing(self, G=None):
		return cg.edit_to_cograph(G) if G else cg.edit_to_cograph(self._G_perturbed)
		

	def color_editing(self, G=None):
		if G:
			nodes = G.nodes()
		else:
			nodes = self._G.nodes()
		copy_G = self._G.copy()
		for e in G.edges():
			if nodes[e[0]]['color'] == nodes[e[1]]['color']:
				copy_G.remove_edge(e)
		return copy_G



	#########################################################################################################
	#																										#
	#											PRINT DATA													#
	#																										#
	#########################################################################################################


	def print_P3_data(self):
		
		pass

	def print_symmetric_diff(self, G):
		# edit distance between edited graph and perturbed graph
		edit_dist = symmetric_diff(self._G_perturbed, G)
		print("The edit distance is: {}".format(edit_dist))
		

	def print_perturbed_G_data(self):
		print("The amount of nodes and edges in the perturbed graph is: {}, {}".format(len(self._G_perturbed.nodes()), len(self._G_perturbed.edges())))

	def print_data(self):

		# check that the variable being divided by is > 0
		'''
			Need to count the number of times the perturbed graph is neither a cograph nor "consistent" to properly count the frequency at which they are both "fixed" by each edit heuristic
			also how many times the perturbed graph is a cograph but not "consistent" and vice versa.
		'''
		print("\n\t\t------------------------------------Cograph editing data------------------------------------")
		if self._count_dG_notCograph_notConsistent > 0:
			print("\nFrequency of cograph editing making the set of triples compatible: {}".format(self._count_cographEdit_fixed_consistency/self._count_dG_notCograph_notConsistent))
		else:
			print("\nNo perturbed graph ended up not being a cograph and having a non-compatible set of triples.")
		if self._count_dG_notCograph_consistent > 0:
			print("\nFrequency of cograph editing not making the set of triples incompatible: {}".format(self._count_cographEdit_remained_consistent/self._count_dG_notCograph_consistent))
		else:
			print("\nNo perturbed graph ended up not being a cograph and having a set of compatible triples.")
		if self._count_dG_notCograph_consistent > 0:
			print("\nFrequency of cograph editing making the set of triples incompatible: {}".format(self._count_cographEdit_broke_consistency/self._count_dG_notCograph_consistent))
		if self._count_dG_not_cograph > 0:
			print("\nFrequency of cograph editing turning an arbitrary properly colored graph into an LDT-graph: {}".format(self._count_cographEdit_to_LDT/self._count_dG_not_cograph))
		else:
			print("\nAll perturbed graphs remained cographs.")

		#print("The amount of times cograph editing made it not properly colored is: {}".format(100-self._count_cographEdit_remain_properly_colored))



		print("\n\t\t------------------------------------Triples editing data------------------------------------")
		if self._count_dG_notCograph_notConsistent > 0:
			print("\nFrequency of triples editing turning the graph into a cograph: {}".format(self._count_triplesEdit_fixed_cograph/self._count_dG_notCograph_notConsistent))
		if self._count_dG_cograph_notConsistent > 0:
			print("\nFrequency of triples editing not breaking the cograph: {}".format(self._count_triplesEdit_remained_cograph/self._count_dG_cograph_notConsistent))
		if self._count_dG_cograph_notConsistent > 0:
			print("\nFrequency of triples editing breaking the cograph: {}".format(self._count_triplesEdit_broke_cograph/self._count_dG_cograph_notConsistent))
		if self._count_dG_not_consistent > 0:
			print("\nFrequency of triples editing turning an arbitrary properly colored graph into an LDT-graph: {}".format(self._count_triplesEdit_to_LDT/self._count_dG_not_consistent))
			print("\nFrequency of triples editing making the graphs set of triples compatible (not introducing new inconsistent P3s): {}.".format(self._count_triplesEdit_success/self._count_dG_not_consistent))
		
		#print("The edit distance (symmetric difference) between the perturbed graph and the triples edited graph is: {}".format())			
		#print("The edit distance (symmetric difference) between the perturbed graph and the cograph edited graph is: {}".format())

		print("\nOf 100 perturbed graphs, the amount of times it wasn't a cograph is: {}".format(self._count_dG_not_cograph))
		print("Of 100 perturbed graphs, the amount of times its set of triples wasn't compatible is: {}".format(self._count_dG_not_consistent))

if __name__ == "__main__":
	pass