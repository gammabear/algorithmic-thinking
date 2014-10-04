"""
Template for Project 3
Student will implement four functions:

slow_closest_pairs(cluster_list)
fast_closest_pair(cluster_list) - implement fast_helper()
hierarchical_clustering(cluster_list, num_clusters)
kmeans_clustering(cluster_list, num_clusters, num_iterations)

where cluster_list is a list of clusters in the plane
"""

import math
import alg_cluster
import urllib2

#########################################################################
# Helper function for converting list of clusters to set of county tuples

def set_of_county_tuples(cluster_list):
    """
    Input: A list of Cluster objects
    Output: Set of sorted tuple of counties corresponds to counties in each cluster
    """
    set_of_clusters = set([])
    for cluster in cluster_list:
        ##print cluster
        #print "cluster in cluster list", cluster
        
        counties_in_cluster = cluster.fips_codes()
        
        # convert to immutable representation before adding to set
        county_tuple = tuple(sorted(list(counties_in_cluster)))
        set_of_clusters.add(county_tuple)
    return set_of_clusters


#############################################################################

# Assets
DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
DATA_3108_URL = DIRECTORY + "data_clustering/unifiedCancerData_3108.csv"
DATA_111_URL = DIRECTORY + "data_clustering/unifiedCancerData_111.csv"

DATA_24_URL = DIRECTORY + "data_clustering/unifiedCancerData_24.csv"


############################################################
# Load data tables

def load_data_table(data_url):
    """
    Import a table of county-based cancer risk data
    from a csv format file
    """
    data_file = urllib2.urlopen(data_url)
    data = data_file.read()
    data_lines = data.split('\n')
    print "Loaded", len(data_lines), "data points"
    data_tokens = [line.split(',') for line in data_lines]
    return [[tokens[0], float(tokens[1]), float(tokens[2]), int(tokens[3]), float(tokens[4])] 
            for tokens in data_tokens]



def pair_distance(cluster_list, idx1, idx2):
    """
    Helper function to compute Euclidean distance between two clusters
    in cluster_list with indices idx1 and idx2
    
    Returns tuple (dist, idx1, idx2) with idx1 < idx2 where dist is distance between
    cluster_list[idx1] and cluster_list[idx2]
    """
    return (cluster_list[idx1].distance(cluster_list[idx2]), idx1, idx2)


def slow_closest_pairs(cluster_list):
    """
    Compute the set of closest pairs of cluster in list of clusters
    using O(n^2) all pairs algorithm
    
    Returns the set of all tuples of the form (dist, idx1, idx2) 
    where the cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.   
    
    """
    result = set([(float("inf"), -1, -1)])
    min_distance = float("inf")
    for idx1 in range(0, len(cluster_list)):
        for idx2 in range(idx1+1, len(cluster_list)):
            pair_dist = pair_distance(cluster_list, idx1, idx2)[0]
            if pair_dist < min_distance:
                min_distance = pair_dist
                result = set([(min_distance, idx1, idx2)])
            if pair_dist == min_distance:
                result.add((min_distance, idx1, idx2))
    return result


def fast_closest_pair(cluster_list):
    """
    Compute a closest pair of clusters in cluster_list
    using O(n log(n)) divide and conquer algorithm
    
    Returns a tuple (distance, idx1, idx2) with idx1 < idx 2 where
    cluster_list[idx1] and cluster_list[idx2]
    have the smallest distance dist of any pair of clusters
    """
        
    def fast_helper(cluster_list, horiz_order, vert_order):
        """
        Divide and conquer method for computing distance between closest pair of points
        Running time is O(n * log(n))
        
        horiz_order and vert_order are lists of indices for clusters
        ordered horizontally and vertically
        
        Returns a tuple (distance, idx1, idx2) with idx1 < idx 2 where
        cluster_list[idx1] and cluster_list[idx2]
        have the smallest distance dist of any pair of clusters
    
        """
        # base case
        size_h = len(horiz_order)
        if size_h <= 3:          
            s_list = []
            for idx in horiz_order:
                s_list.append(cluster_list[idx])
            (distance, idx, jdx) = slow_closest_pairs(s_list).pop()
            return (distance, horiz_order[idx], horiz_order[jdx])
            
        # divide
        mid_idx = int(math.ceil(size_h/2)) 
        mid_hcoord =  0.5 * (cluster_list[horiz_order[mid_idx-1]].horiz_center() + cluster_list[horiz_order[mid_idx]].horiz_center())  # bug
        horiz_order_left = horiz_order[0:mid_idx]
        horiz_order_right = horiz_order[mid_idx:size_h]

        
        vcoord_and_index_left = [(cluster_list[idx].vert_center(), idx) 
                        for idx in horiz_order_left]
        vcoord_and_index_left.sort()
        vert_order_left = [vcoord_and_index_left[idx][1] for idx in range(len(vcoord_and_index_left))]

        vcoord_and_index_right = [(cluster_list[idx].vert_center(), idx) 
                        for idx in horiz_order_right] 
        vcoord_and_index_right.sort()
        vert_order_right = [vcoord_and_index_right[idx][1] for idx in range(len(vcoord_and_index_right))]

        closest_pair_left = fast_helper(cluster_list, horiz_order_left, vert_order_left)
        closest_pair_right = fast_helper(cluster_list, horiz_order_right, vert_order_right)
        if closest_pair_left[0] < closest_pair_right[0]:
            closest_pair = closest_pair_left
        else:
            closest_pair = closest_pair_right

        # conquer
        s_list = []
        for idx in vert_order:
            if cluster_list[idx].horiz_center() - mid_hcoord < closest_pair[0]:
                s_list.append(idx)

        for idx in range(0, len(s_list) - 1):
            for jdx in range(idx + 1, min(idx + 4, len(s_list))):
                distance = cluster_list[s_list[idx]].distance(cluster_list[s_list[jdx]])
                if distance < closest_pair[0]:
                    closest_pair = (distance, s_list[idx], s_list[jdx])
  
        return closest_pair
            
    # compute list of indices for the clusters ordered in the horizontal direction
    hcoord_and_index = [(cluster_list[idx].horiz_center(), idx) 
                        for idx in range(len(cluster_list))] 
    
    hcoord_and_index.sort()
    
    horiz_order = [hcoord_and_index[idx][1] for idx in range(len(hcoord_and_index))]
    
    # compute list of indices for the clusters ordered in vertical direction
    vcoord_and_index = [(cluster_list[idx].vert_center(), idx) 
                        for idx in range(len(cluster_list))]    
    vcoord_and_index.sort()
    
    vert_order = [vcoord_and_index[idx][1] for idx in range(len(vcoord_and_index))]
    
    # compute answer recursively
    answer = fast_helper(cluster_list, horiz_order, vert_order) 
    
    return (answer[0], min(answer[1:]), max(answer[1:]))

    

def hierarchical_clustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function mutates cluster_list
    
    Input: List of clusters, number of clusters
    Output: List of clusters whose length is num_clusters
    """
    while len(cluster_list) > num_clusters:
        closest_pair = fast_closest_pair(cluster_list)
        cluster_i = cluster_list[closest_pair[1]]
        cluster_j = cluster_list[closest_pair[2]]
        cluster_i_union_j = cluster_i.merge_clusters(cluster_j)
        cluster_list.append(cluster_i_union_j)
        cluster_list.remove(cluster_i)
        cluster_list.remove(cluster_j)
    return cluster_list


    
def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters
    Note: the function mutates cluster_list
    
    Input: List of clusters, number of clusters, number of iterations
    Output: List of clusters whose length is num_clusters
    """
    
    # 2 initialize k-means clusters to be initial clusters with largest populations
    population_and_index = [(cluster_list[idx].total_population(), idx) 
                        for idx in range(len(cluster_list))] 
    #print "population", population_and_index
    
    population_and_index.sort()
    
    #print "population", population_and_index
    
    population_order = [population_and_index[idx][1] for idx in range(len(population_and_index))]
    
    #print "pop order", population_order
    
    #print "number of clusters", len(cluster_list)
    #print "num clusters we want", num_clusters
    #print "number of iterations", num_iterations
    
    centers_list = []
    for idx in range(len(cluster_list), len(cluster_list)-num_clusters, -1):
        #print cluster_list[population_order[idx-1]]  #clusters sorted by pop
        (hcoord, vcoord) = cluster_list[population_order[idx-1]].horiz_center(), cluster_list[population_order[idx-1]].vert_center()
        #print "(x,y)", (hcoord, vcoord)
        centers_list.append((hcoord, vcoord))
    
    #print "centers", centers_list    
 
    #main for loop
    for idx in range(0, num_iterations):
        # initialize k empty clusters
        k_list = [alg_cluster.Cluster(set([]),centers_list[idx][0], centers_list[idx][1], 0, 0.0) for idx in range(0, num_clusters)] 
        answer_list = [alg_cluster.Cluster(set([]),centers_list[idx][0], centers_list[idx][1], 0, 0.0) for idx in range(0, num_clusters)]
        #print "OG answer_list", answer_list
        #print "len of answer", len(answer_list)

        for jdx in range(0, len(cluster_list)):
            min_distance = float("inf")
            min_kdx = -1
            for kdx in range(0, num_clusters):
                #print "jdx, kdx", jdx, kdx
                #print "anwer_list, cluster_list", answer_list[kdx], cluster_list[jdx]
                #if jdx == 11 and kdx == 5:
                    #print "answer list", answer_list
                    #print "answer_list5", kdx, answer_list[kdx]
                    #print "clusterlist11", cluster_list[jdx]
                distance = k_list[kdx].distance(cluster_list[jdx])
                #print "distance from", jdx, "to", kdx, distance
                if distance < min_distance:
                    min_distance = distance
                    min_kdx = kdx
            
            #print "for jdx=", jdx, "min kdx=", min_kdx
            answer_list[min_kdx].merge_clusters(cluster_list[jdx])
            #print "merged answer list", answer_list
        
        # recompute its center
        for kdx in range(0, num_clusters):
            #print "current center", centers_list[kdx]
            #print "new horiz", answer_list[kdx].horiz_center()
            #print "new vert", answer_list[kdx].vert_center()
            
            (new_hcoord, new_vcoord) = answer_list[kdx].horiz_center(), answer_list[kdx].vert_center()
            centers_list[kdx] = (new_hcoord, new_vcoord)
       
    return answer_list

def test_fast():
    """
    Test for fast_closest_pair
    kmeans_clustering should not mutate cluster_list, but make a new copy of each test anyways
    """
    data_24_table = load_data_table(DATA_24_URL)
    cluster_list = []
    for idx in range(len(data_24_table)):
        line = data_24_table[idx]
        cluster_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))  
    
    #print "OG cluster", cluster_list
    print "distance, 1,2", cluster_list[1].distance(cluster_list[2])
    
    print "222slow closest", slow_closest_pairs(cluster_list)
    print "222fast closest", fast_closest_pair(cluster_list)

def test_owl():
    """
    Test for fast_closest_pair
    kmeans_clustering should not mutate cluster_list, but make a new copy of each test anyways
    """
    #cluster_list =  [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0), alg_cluster.Cluster(set([]), 24.1715119657, -31.4764880007, 1, 0), alg_cluster.Cluster(set([]), 55.9347036476, -48.5578980541, 1, 0), alg_cluster.Cluster(set([]), -45.668221279, -40.2204165013, 1, 0), alg_cluster.Cluster(set([]), -92.1133508388, 41.1264397769, 1, 0), alg_cluster.Cluster(set([]), -61.335334931, -39.9238885371, 1, 0), alg_cluster.Cluster(set([]), -41.3587295901, 40.3245981796, 1, 0), alg_cluster.Cluster(set([]), -12.4948622123, 98.3102571423, 1, 0), alg_cluster.Cluster(set([]), 84.331269134, 87.4096437485, 1, 0), alg_cluster.Cluster(set([]), -64.3793374021, -29.3604139823, 1, 0), alg_cluster.Cluster(set([]), 88.1863663091, 87.1319063356, 1, 0), alg_cluster.Cluster(set([]), -17.1454765601, 63.574041831, 1, 0), alg_cluster.Cluster(set([]), 94.065223686, -24.2896714813, 1, 0), alg_cluster.Cluster(set([]), 41.3803043889, 33.56171226, 1, 0), alg_cluster.Cluster(set([]), 85.2641824191, 52.6931063404, 1, 0), alg_cluster.Cluster(set([]), -56.3981249563, -83.5611722575, 1, 0), alg_cluster.Cluster(set([]), -67.2363423428, -56.91272576, 1, 0), alg_cluster.Cluster(set([]), 28.8357363013, 14.416522777, 1, 0)]
    #cluster_base =  [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0)]
    #cluster_four =  [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0), alg_cluster.Cluster(set([]), 24.1715119657, -31.4764880007, 1, 0)]
    #cluster_five =  [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0), alg_cluster.Cluster(set([]), 24.1715119657, -31.4764880007, 1, 0), alg_cluster.Cluster(set([]), 55.9347036476, -48.5578980541, 1, 0)]
    #cluster_six =   [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 28.8357363013, 14.416522777, 1, 0)]
    #cluster_seven = [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0),alg_cluster.Cluster(set([]), 28.8357363013, 14.416522777, 1, 0)]
    #cluster_eight = [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0), alg_cluster.Cluster(set([]), 24.1715119657, -31.4764880007, 1, 0), alg_cluster.Cluster(set([]), 28.8357363013, 14.416522777, 1, 0)]
    #cluster_nine =  [alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0), alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0), alg_cluster.Cluster(set([]), 24.1715119657, -31.4764880007, 1, 0), alg_cluster.Cluster(set([]), 55.9347036476, -48.5578980541, 1, 0),alg_cluster.Cluster(set([]), 28.8357363013, 14.416522777, 1, 0)]
    #cluster_nine_left = [alg_cluster.Cluster(set([]), -67.5757896921, -6.52776165362, 1, 0), alg_cluster.Cluster(set([]), 24.1715119657, -31.4764880007, 1, 0), alg_cluster.Cluster(set([]), 27.9189799338, 17.6102324623, 1, 0)]
    #cluster_nine_right = [alg_cluster.Cluster(set([]), 28.8357363013, 14.416522777, 1, 0), alg_cluster.Cluster(set([]), 55.9347036476, -48.5578980541, 1, 0), alg_cluster.Cluster(set([]), 86.0830468948, -59.1167021002, 1, 0)]
    #print "slow closest" , slow_closest_pairs(cluster_list)
    #print "slow base", slow_closest_pairs(cluster_base)
    #print "fast base", fast_closest_pair(cluster_base)
    #print "slow five", slow_closest_pairs(cluster_list)
    #print "fast five", fast_closest_pair(cluster_list)

    
#test_fast()

test_owl()

def test_kmeans():
    """
    Test for hierarchial
 
    """
    data_24_table = load_data_table(DATA_24_URL)
    cluster_list = []
    for idx in range(len(data_24_table)):
        line = data_24_table[idx]
        cluster_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
        
    print "OG cluster", cluster_list
    #print "h clustering answer", hierarchical_clustering(cluster_list, 20)
    student_clustering = kmeans_clustering(cluster_list, 15, 1)
    student_county_tuple = set_of_county_tuples(student_clustering)
    print student_county_tuple

#test_kmeans()



