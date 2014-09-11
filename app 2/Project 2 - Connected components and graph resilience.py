"""
Project 2 - Connected components and graph resilience
"""
import poc_queue

EX_GRAPH0 = {0: set([1,2]),
             1: set([]),
             2: set([])}
EX_GRAPH1 = {0: set([1,4,5]),
             1: set([2,6]),
             2: set([3]),
             3: set([]),
             4: set([1]),
             5: set([2]),
             6: set([])}
EX_GRAPH2 = {0: set([1,4,5]),
             1: set([2,6]),
             2: set([3,7]),
             3: set([7]),
             4: set([1]),
             5: set([2]),
             6: set([]),
             7: set([3]),
             8: set([1,2]),
             9: set([0,3,4,5,6,7])}

def bfs_visited(ugraph, start_node):
    """
    Takes the undirected graph ugraph and the node start_node
    and returns the set consisting of all nodes that are 
    visited by a breadth-first search that starts at start_node.
    """
    bfs_queue = poc_queue.Queue()
    visited = set([])
    visited.add(start_node)
    bfs_queue.enqueue(start_node)
    while len(bfs_queue) != 0:
        node_j = bfs_queue.dequeue()
        neighbors = ugraph[node_j]
        for nbr in neighbors:
            if nbr not in visited:
                visited.add(nbr)
                bfs_queue.enqueue(nbr)
    return visited

def cc_visited(ugraph):
    """
    Takes the undirected graph ugraph and returns a list of 
    sets, where each set consists of all the nodes (and 
    nothing else) in a connected component, and there 
    is exactly one set in the list for each connected 
    component in ugraph and nothing else.
    """
    remaining_nodes = ugraph.keys()
    connected_comp_list = []
    while len(remaining_nodes) != 0:
        working_set = bfs_visited(ugraph, remaining_nodes[0])
        connected_comp_list.append(working_set)
        for node in working_set:
            remaining_nodes.remove(node)
    return connected_comp_list
        
def largest_cc_size(ugraph):
    """
    Takes the undirected graph ugraph and returns the size
    (an integer) of the largest connected component in 
    ugraph.
    """
    largest = 0;
    connected_comp_list = cc_visited(ugraph)
    for connected_comp in connected_comp_list:
        size = len(connected_comp)
        if size > largest:
            largest = size
    return largest
    
def compute_resilience(ugraph, attack_order):
    """
    Takes the undirected graph ugraph, a list of nodes 
    attack_order and iterates through the nodes in 
    attack_order. 
    """
    resilience_list = []
    size = largest_cc_size(ugraph)
    resilience_list.append(size)
    for attacked_node in attack_order:
        ugraph.pop(attacked_node)
        for dummy_node in ugraph:
            if attacked_node in ugraph[dummy_node]:
                ugraph[dummy_node].remove(attacked_node)
        size = largest_cc_size(ugraph)
        resilience_list.append(size)
    return resilience_list
            