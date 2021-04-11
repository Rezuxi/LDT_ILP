"""
	ILP formulation for LDT editing

	The constants:
		E(x, y) = 1 if there is an edge (x, y) in G, else 0


	The variables:
		a, b, c are the colors of the nodes (x, y, z)
		T(a, b, c) = 1 if there is a P3 on (x, y, z) with colors (a, b, c) in G*
		
		e(x, y) = 1 if there is an edge (x, y) in G* and a/=b , else 0
		e(x, y)=e(y, x)
		


	The objective function:

	The constraints:
		(1)	for all ordered (w, x, y, z):
			e(w, x) + e(x, y) + e(y, z) - e(x, z) - e(w, z) - e(w, y) <= 2
			(cograph)
		(2)	for all ordered (a, b, c, d):
			2T(a, b, c) + 2T(a, d, b) - T(b, d, c) - T(a, d, c) <= 2
			(consistency)
		(3) for all (a, b, c):
			T(a, b, c) + T(a, c, b) + T(b, c, a) = 1
			(strictly dense)
		(4) for all ordered (x, y, z):
			e(x, y) + e(y, z) + (1 - e(x, z)) - T(a, c, b) <= 2
			(a, b, c) is the coloring for (x, y, z).
			(triples)



"""


import itertools
import gurobipy as gp
from gurobipy import GRB

import networkx as nx
import asymmetree.tools.GraphTools as gt

import os
import json
from networkx.readwrite import json_graph
import re





class LDTEditor:




	def __init__(self, G, only_add=False, only_delete=False):
		self.graph = G
		self.color_dict = gt.sort_by_colors(G)
		self._only_add = only_add
		self._only_delete = only_delete
		#print("The colors are:")
		#print([*self.color_dict])



	def build_model(self):

		# set constants
		self.E = {(x, y): int(self.graph.has_edge(x, y))
					for x in self.graph.nodes()
					for y in self.graph.nodes()
					if x != y}
		# build model with variables and objective function
		self.model = gp.Model()

		# add variables for edges of edited graph
		self.e = self.model.addVars(self.E, vtype=GRB.BINARY, name='e')
		# add variables for triples
		self._additional_variables()
		#self.model.update()


		self._objective()

		# set constraints
		if self._only_add:
			self._only_insertion_constraint()
		elif self._only_delete:
			self._only_deletion_constraint()
		self._undirected_constraint()
		self._properly_colored_constraint()
		self._cograph_constraint()
		self._strictly_dense_constraint()
		self._consitency_constraint()
		self._triples_constraint()

	def _objective(self):
		self.obj = gp.LinExpr()
		"""

		# the obj func (1 - e_xy) * E_xy + (1 - E_xy) * e_xy for all x, y
		# for directed graphs
		for x, y in self.E:
			if self.E[x, y]:
				self.obj += 1 - self.e[x, y]
			else:
				self.obj += self.e[x, y]
		"""
		# because the graph is undirected we only count each edge once. e_xy=e_yx
		for x, y in itertools.combinations(self.graph.nodes(), 2):
			if self.E[x, y]:
				self.obj += 1 - self.e[x, y]
			else:
				self.obj += self.e[x, y]

		self.model.setObjective(self.obj, GRB.MINIMIZE)

	

	def _additional_variables(self):
		
		self.triples = ((a, b, c)
						for a, b, c in itertools.permutations([*self.color_dict], 3)
						if a < b)
		self.t = self.model.addVars(self.triples, vtype=GRB.BINARY, name='t')


		self.model.update()


	######################################################################################################################
	#																													 #
	#													CONSTRAINTS														 #
	#																													 #
	######################################################################################################################

	# constraint 0
	def _undirected_constraint(self):
		for x, y in self.E:
			self.model.addConstr((self.e[x, y] == self.e[y, x]), name='undir[({}, {}) = ({}, {})]'.format(x, y, y, x))


	# constraint 1
	def _properly_colored_constraint(self):
		"""
			set forbidden edges to preserve the graph being properly colored
		"""
		for a in [*self.color_dict]:
			# get all nodes with color a
			# enforce e variable between all x,y in a to be 0. i.e. there can be no edge between any nodes in the same color set.
			nodes = self.color_dict[a]

			if len(nodes) >= 2:
				for x, y in itertools.combinations(nodes, 2):
					self.model.addConstr((self.e[x, y] == 0), 
										name='pc[{}, {}]'.format(x, y))

	# constraint 2
	def _cograph_constraint(self):
		for w, x, y, z in itertools.permutations(self.graph.nodes(), 4):
			self.model.addConstr((self.e[w, x] + self.e[x, y] + self.e[y, z] - 
								self.e[x, z] - self.e[w, z] - self.e[w, y] <= 2),
								name='c[{}, {}, {}, {}]'.format(w, x, y, z))

		
	# constraint 3
	def _strictly_dense_constraint(self):
		for a, b, c in itertools.combinations(sorted([*self.color_dict]), 3):

			ab = sorted([a, b])
			ac = sorted([a, c])
			bc = sorted([b, c])

			self.model.addConstr((self.t[(*ab, c)] + self.t[(*ac, b)] + self.t[(*bc, a)] == 1), 
								name='sd[{}, {}, {}]'.format(a, b, c))
		

	# constraint 4
	def _consitency_constraint(self):
		for a, b, c, d in itertools.permutations([*self.color_dict], 4):

			ab = sorted([a, b])
			ad = sorted([a, d])
			bd = sorted([b, d])

			self.model.addConstr((2*self.t[(*ab, c)] + 2*self.t[(*ad, b)] - self.t[(*bd, c)] - self.t[(*ad, c)] <= 2), 
								name='tc[{}, {}, {}, {}]'.format(a, b, c, d))


	# constraint 5
	def _triples_constraint(self):
		for x, y, z in itertools.permutations(self.graph.nodes(), 3):
			
			# this constraint should be set for x, y, z with distinct colors
			# so if any of these 3 colors are the same then we dont add a constraint for this x, y, z.
			x_color = self.graph.nodes[x]['color']
			y_color = self.graph.nodes[y]['color']
			z_color = self.graph.nodes[z]['color']
			#print("a, b, c: {}, {}, {}".format(x_color, y_color, z_color))
			if (x_color != y_color and x_color != z_color) and (y_color != z_color):
				xz_color  = sorted([x_color, z_color])
				self.model.addConstr((self.e[x, y] + self.e[y, z] + (1 - self.e[x, z]) - self.t[(*xz_color, y_color)] <= 2), 
									name='triple[({}, {}, {}), ({}, {}, {})]'.format(x, y, z, x_color, y_color, z_color))


	def _only_insertion_constraint(self):
		for x, y in itertools.combinations(self.graph.nodes(), 2):
			self.model.addConstr((self.e[x, y] >= self.E[x, y]), name='oi[{}, {}]'.format(x, y))
		
	def _only_deletion_constraint(self):
		for x, y in itertools.combinations(self.graph.nodes(), 2):
			self.model.addConstr((self.e[x, y] <= self.E[x, y]), name='od[{}, {}]'.format(x, y))


	def optimize(self, time_limit=False):
		if time_limit:
			original_time_limit = self.model.Params.timeLimit
			self.model.Params.timeLimit = time_limit
		
		self.model.optimize()
		
		if time_limit:
			self.model.Params.timeLimit = original_time_limit


	def get_solution(self):

		sol_graph = nx.Graph()
		sol_graph.add_nodes_from(self.graph.nodes.data())
		# set edges (only need combinations of x, y since e[x, y]=e[y, x])
		for x, y in itertools.combinations(sol_graph.nodes(), 2):
			#print("e_{},{}: {}".format(x, y, self.e[x, y].X))
			if self.e[x, y].X:
				sol_graph.add_edge(x, y)
		return sol_graph, self.model.objVal


	def get_solve_time(self):
		return self.model.Runtime


	def _save_ILP_data(self, G, edited_G, solve_time, min_edit_dist, only_add=False, only_delete=False, filename='LDTEdit_exact_solution'):
		n = G.order()
		ILP_data = {'solve_time': solve_time, 'n': n, 'only_add': only_add, 'only_delete': only_delete, 'min_edit_dist': min_edit_dist,
					'G': json_graph.node_link_data(G), 'edited_G': json_graph.node_link_data(edited_G)
					}

		existing_files = []

		for _, _, files in os.walk('./exact_results'):
			for file in files:
				existing_files.append(file)

		def find_ID(f):
			s = re.findall('(?<=\_)([0-9]*?)(?=\.)', f)
			return int(s[0]) if s else -1

		next_ID = 0
		for f in existing_files:
			ID = find_ID(f)
			if ID >= next_ID:
				next_ID = ID + 1

		new_name = 'exact_results/' + filename + "_{}_{}_{}_{}.json".format(n, int(only_add), int(only_delete), next_ID)
		
		# TODO: handle errors
		with open(new_name, 'w', encoding='utf-8') as f:
			json.dump(ILP_data, f, ensure_ascii=False, indent=6)

		print("The results of the ILP were successfully saved to {}".format(new_name))

	@staticmethod
	def get_ILP_data(filename):
		# TODO: handle errors
		with open(filename) as f:
			ILP_data = json.load(f)

		G = json_graph.node_link_graph(ILP_data['G'])
		edited_G = json_graph.node_link_graph(ILP_data['edited_G'])
		
		only_add = ILP_data['only_add']
		only_delete = ILP_data['only_delete']
		min_edit_dist = ILP_data['min_edit_dist']
		#solve_time = ILP_data['solve_time']

		return G, edited_G, only_add, only_delete, min_edit_dist





if __name__ == "__main__":

	# non ldt graph
	G = nx.Graph()
	G_nodes = [(0, {"color": 0}), (1, {"color": 1}), (2, {"color": 2}), (3, {"color": 0}),
			   (4, {"color": 0}) 
			  ]
	G_edges = [(0, 1), (1, 2), (2, 3), (2, 4)]

	G.add_nodes_from(G_nodes)
	G.add_edges_from(G_edges)

	solver = LDTEditor(G)
	solver.build_model()
	solver.optimize(time_limit=None)

	sol_graph, sol_distance = solver.get_solution()
	# sol_distance will be 2*symmetric_diff if using obj func for directed graph since we're counting xy edges and yx edges in the objective function. the graph is undirected so only xy or yx is needed


	edit_dist = gt.symmetric_diff(G, sol_graph)
	print("------------------------------------------------------------------------------------------------------------------------")
	print("The symmetric difference between G and G* is: {}".format(edit_dist))
	print("The value for the objective function is: {}".format(sol_distance))
	#print("------------------------------------------------------------------------------------------------------------------------")
	#print("The edges of the solution graph: {}".format(sol_graph.edges()))
	print("------------------------------------------------------------------------------------------------------------------------")
	print("Runtime: {}".format(solver.get_solve_time()))
	print("Saving data...")
	solver._save_ILP_data(G, sol_graph, solver.get_solve_time(), edit_dist, only_add=False, only_delete=False)