from LDT_ILP import LDTEditor
import asymmetree.treeevolve as te
from asymmetree.datastructures import PhyloTree
from asymmetree.hgt import ldt_graph
import asymmetree.tools as tools
from asymmetree.tools.GraphTools import disturb_graph
import asymmetree.tools.GraphTools as gt
from GraphTools import *

S = te.simulate_species_tree(10, model='innovation')
TGT = te.simulate_dated_gene_tree(S, dupl_rate=0.5, loss_rate=0.5, hgt_rate=0.5, prohibit_extinction='per_family', replace_prob=0.0)
OGT = te.observable_tree(TGT)
ldt = ldt_graph(OGT, S)

print("Amount of nodes in the graph: {}".format(len(ldt.nodes())))

if __name__ == '__main__':
	'''
		given and ldt graph, add noise and use ILP solver on perturbed graph
	'''
	G = InvestigateGraph(ldt)

	# this makes sure we get a perturbed graph that is not an LDT graph.
	G.perturb_graph()

	solver = LDTEditor(G._G_perturbed)
	solver.build_model()
	solver.optimize(time_limit=None)

	sol_graph, sol_distance = solver.get_solution()

	edit_dist = gt.symmetric_diff(G._G_perturbed, sol_graph)
	print("The value of the ILP: {}".format(sol_distance))
	print("The value of the symmetric distance: {}".format(edit_dist))
	print("Runtime: {}".format(solver.get_solve_time()))
	print("Saving data...")
	solver._save_ILP_data(G._G, sol_graph, solver.get_solve_time(), edit_dist, only_add=False, only_delete=False)

	properly_colored = is_properly_colored(sol_graph)
	cograph = is_cograph(sol_graph)
	compatible = is_compatible(sol_graph)

	if properly_colored and cograph and compatible:
		print("Successfully edited into an LDT-graph!")
	else:
		print("The editing was unsuccessful!")

	
