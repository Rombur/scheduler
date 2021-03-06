# Python code
# Author: Bruno Turcksin
# Date: 2015-01-25 16:31:13.044256

#----------------------------------------------------------------------------#
## Classes BRANCH_AND_BOUND and NODE                                        ##
#----------------------------------------------------------------------------#

"""This module contains the class Branch_and_bound that solve the SOP using
banch-and-bound algorithm and the class Node that represent the different node
of the graph created by the branch-and-bound algorithm. (This algorithm is
actually closer to A* than Branch-And-Bound)"""

import itertools
import numpy as np
from sklearn.cluster import MiniBatchKMeans
import random
import TASK


class NODE(object):
    """Class representing a node in the tree created by the branch-and-bound
    algorithm."""

    def __init__(self, graph, tasks_done, tasks_ready, cost, min_bound):

        super().__init__()
        self.graph = graph
        self.tasks_done = tasks_done
        self.tasks_ready = tasks_ready
        self.cost = cost
        self.min_bound = min_bound

#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#


class BRANCH_AND_BOUND(object):
    """Class solving the SOP using branch-and-bound algorithm."""

    def __init__(self, tasks, n_processors, best_time, output, n_clusters,
            bound_estimator_type):

        super().__init__()
        self.tasks = tasks
        self.n_tasks = len(tasks)
        self.n_processors = n_processors
        self.best_time = best_time
        self.output = output
        self.n_clusters = n_clusters
        self.bound_estimator_type = bound_estimator_type
        self.current_directions = []
# Create a map whose the key is the task_id and the value is the position of the
# task in self.tasks
        self.task_id_to_pos = [0 for i in range(self.n_tasks)]
        for i in range(self.n_tasks):
            self.task_id_to_pos[self.tasks[i].task_id] = i

#----------------------------------------------------------------------------#

    def Compute_b_level_recursive(self, task):
        """Recursive function that computes the b-level of the task."""

        b_level = task.delta_t
# The size of the list needs to be at least one for the max function to work.
        downstream_b_level = [0 for i in range(len(task.waiting_tasks)+1)]
        pos = 0
        for waiting_task_id in task.waiting_tasks:
            waiting_task = self.tasks[self.task_id_to_pos[waiting_task_id]]
            if waiting_task.b_level != 0:
                downstream_b_level[pos] = waiting_task.b_level
            else:
                downstream_b_level[pos] = self.Compute_b_level_recursive(waiting_task)
            pos += 1

# Tuple are immutable so we need to replace task
        b_level += max(downstream_b_level)
        pos = self.task_id_to_pos[task.task_id]
        self.tasks[pos] = self.tasks[pos]._replace(b_level=b_level)

        return b_level 

#----------------------------------------------------------------------------#

    def Compute_b_level(self):
        """Compute the b-level of all the tasks."""

        pos = 0
        for task in self.tasks:
            if len(task.required_tasks) == 0:
                self.Compute_b_level_recursive(task)

#----------------------------------------------------------------------------#

    def Compute_bound_estimator(self, tasks_ready, tasks_done):
        """Compute a minimum bound estimator. The 'adams' bound estimator uses
        the optimal sweep as a bound. The 'b-level' bound estimator uses the
        b-level of the ready tasks."""
        
        if self.bound_estimator_type == 'adams':
#min_bound = max n_tasks left for a processor + 
#            max(2*sqrt(procs that haven't executed)-4,0) +
#            max(min(distance to reach proc boundary)
            tasks_left = [0 for i in range(self.n_processors)]
            started_procs = set()
            for task in self.tasks:
                if task.task_id not in tasks_done:
                    tasks_left[task.subdomain_id] += task.delta_t
                else:
                    started_procs.add(task.subdomain_id)
            sqrt_n_idle_procs = np.sqrt(self.n_processors-len(started_procs))
            min_dist_to_new_proc = self.Compute_min_distance_to_new_proc(current_level)
            bound_estimator = max(tasks_left) + max(2*sqrt_n_idle_procs-4.0)+\
                    max(min_dist_to_new_proc)
        else:
# The bound is computed as follows: for each quadrant compute the minimum
# b_level and add the largest number of ready tasks on a same processors. The
# bound is the maximum of these "effective" b-level.
            min_b_level = [1e6, 1e6, 1e6, 1e6]
            serial_tasks = [ [0 for i in range(self.n_processors)] for j in range(4)]
            first_serial_tasks = [ [True for i in range(self.n_processors)] for j in range(4)]
            for task in tasks_ready:
                dir = task.dir%4
                if task.b_level < min_b_level[dir]:
                    min_b_level[dir] = task.b_level
                if first_serial_tasks[dir][task.subdomain_id] is True:
                    first_serial_tasks[dir][task.subdomain_id] = False
                else:
                    serial_tasks[dir][task.subdomain_id] += task.delta_t

            bound_estimator_list = [0 for i in range(4)]
            for dir in range(4):
                bound_estimator_list[dir] = min_b_level[dir] +\
                        max(serial_tasks[dir])
            bound_estimator = max(bound_estimator_list)

        return bound_estimator

#----------------------------------------------------------------------------#

    def Compute_subdomains_list(self,node):
        """Compute the minimal subdomains list."""

        subdomains_set = set()
        for task in node.tasks_ready:
            subdomains_set.add(task.subdomain_id)

        return list(subdomains_set)

#----------------------------------------------------------------------------#

    def Compute_min_distance_recursive(self,task):
        """This is a recursive function that computes the minimum number of
        stages (distance) before a task on a new processor can bexecuted."""

        distance = task.delta_t
        new_distance = [0 for i in range(len(task.waiting_tasks))]
        pos = 0
        for waiting_task_id in task.waiting_tasks:
            waiting_task = self.tasks[self.task_id_to_pos[waiting_task_id]]
            if waiting_task.subdomain_id == task.subdomain_id:
                new_distance[pos] = self.Compute_min_distance_recursive(waiting_task)
            pos += 1

        return distance + min(new_distance)

#----------------------------------------------------------------------------#

    def Compute_min_distance_to_new_proc(self,current_level):
        """For each processor, compute the minimum number of stages before a
        task on a new processor can be executed."""

        min_distance = [0 for i in range(self.n_processors)]
        for task in current_level:
            for waiting_task_id in task.waiting_tasks:
                waiting_taks = self.tasks[self.task_id_to_pos[waiting_task_id]]
                if waiting_task.subdomain_id == task.subdomain_id:
                    distance = self.Compute_min_distance_recursive(waiting_task)
                    min_distance[task.subdomain_id] = distance

        return min_distance

#----------------------------------------------------------------------------#

    def Add_tasks_to_new_tasks_ready(self, given_tasks, new_tasks_done, new_tasks_ready):
        """Add new taks to ready_tasks."""

        for task in given_tasks:
            for waiting_task_id in task.waiting_tasks:
                waiting_task = self.tasks[self.task_id_to_pos[waiting_task_id]]
                ready = True
                for required_task_id in waiting_task.required_tasks:
                    if required_task_id not in new_tasks_done:
                        ready = False
                        break
                if ready is True:
                    if waiting_task not in new_tasks_ready:
                        new_tasks_ready.append(waiting_task)

#----------------------------------------------------------------------------#

    def Create_node(self, given_procs, given_tasks, parent_node):
        """Create a new node given a processor combination, a task combination,
        and a node."""

        max_delta_t = 0
        new_graph = parent_node.graph[:]
        current_level = []
        new_tasks_done = parent_node.tasks_done.copy()
        new_tasks_ready = parent_node.tasks_ready[:]
        for task in given_tasks:
            current_level.append(task)
            if task.task_id not in new_tasks_done:
                new_tasks_done.append(task.task_id)
            new_tasks_ready.remove(task)
            if max_delta_t < task.delta_t:
                max_delta_t = task.delta_t
        new_graph.append(current_level)
# Add new tasks to the tasks_ready list
        self.Add_tasks_to_new_tasks_ready(given_tasks, new_tasks_done, new_tasks_ready)
# The new cost is the cost of parent_node plus max_delta_t
        new_cost = parent_node.cost + max_delta_t
# Compute new_min_bound
        new_min_bound = new_cost + self.Compute_bound_estimator(new_tasks_ready,
                new_tasks_done)
# Create the new node
        node = NODE(new_graph, new_tasks_done, new_tasks_ready, new_cost, 
                new_min_bound)

        return node

#----------------------------------------------------------------------------#

    def Compute_nodes_centroids(self, nodes_to_explore):
        """Compute clusters of nodes_to_explors using k-means, compute their
        centroids, and then finally return the closest nodes to the
        centroids."""
        
        n_nodes_to_explore = len(nodes_to_explore)
        if n_nodes_to_explore > self.n_clusters and\
                len(self.nodes[nodes_to_explore[0]].graph) != 0 and\
                len(self.nodes[nodes_to_explore[0]].graph[-1]) > 1:

# Build the observation matrix
            obs = np.zeros([n_nodes_to_explore,\
                len(self.nodes[nodes_to_explore[0]].graph[-1])])
            for i in range(n_nodes_to_explore):
                node = self.nodes[nodes_to_explore[i]]
                for j in range(len(node.graph[-1])):
                    obs[i][j] = node.graph[-1][j].task_id

            minibatch_kmeans = MiniBatchKMeans(self.n_clusters)
            minibatch_kmeans.fit(obs)
            distances = minibatch_kmeans.transform(obs)
            reduced_nodes_to_explore = []
            for i in range(self.n_clusters):
                pos = distances[:,i].argmin()
                reduced_nodes_to_explore.append(nodes_to_explore[pos])

            return reduced_nodes_to_explore

        else:
            return nodes_to_explore

#----------------------------------------------------------------------------#

    def Reduce_tasks_ready(self, tasks_ready):
        """Group the tasks that are ready by subdomain id and then, for each
        subdomain id suppress from tasks_ready the equivalent tasks. Tasks are
        equivalent if they corresponds to the same cell and their direction is
        in the same quadrant."""

# Group cells by subdomain id.
        clustered_tasks = []
        for task in tasks_ready:        
            task_added = False
            for cluster in clustered_tasks:
                if task.subdomain_id == cluster[0].subdomain_id:
                    cluster.append(task)
                    task_added = True
            if task_added is False:
                clustered_tasks.append([task])

        reduced_clustered_tasks = []
        for cluster in clustered_tasks:
            equivalent_nodes = set()
            accepted_tasks = []
            for task in cluster:
                i, j = task.pos
# Lists are unhashable, so use a tuple instead
                node_geom_tuple = tuple((i, j, task.dir%4))
                if node_geom_tuple not in equivalent_nodes:
# Set are unhashable, so use a tuple instead. The tuple needs to be sorted
# first.
                    accepted_tasks.append(task)
                    equivalent_nodes.add(node_geom_tuple)
            reduced_clustered_tasks.append(accepted_tasks)

        return reduced_clustered_tasks

#----------------------------------------------------------------------------#

    def Compute_used_tasks_centroids(self, full_used_tasks):
        """Compute clusters of used_tasks using k-means, compute their
        centroids, and then finally return the used_tasks to the
        centroids."""
        
        n_clusters = 10*self.n_clusters
        if len(full_used_tasks) > n_clusters:
# Build the observation matrix
            obs = np.zeros([len(full_used_tasks), len(full_used_tasks[0])])
            for i in range(len(full_used_tasks)):
                used_tasks = full_used_tasks[i]
                for j in range(len(used_tasks)):
                    obs[i][j] = used_tasks[j].task_id

            minibatch_kmeans = MiniBatchKMeans(n_clusters)
            minibatch_kmeans.fit(obs)
            distances = minibatch_kmeans.transform(obs)
            reduced_used_tasks = []
            for i in range(n_clusters):
                pos = distances[:,i].argmin()
                reduced_used_tasks.append(full_used_tasks[pos])

            return reduced_used_tasks

        else:
            return full_used_tasks

#----------------------------------------------------------------------------#

    def Run(self):
        """Solve the SOP using branch-and-bound algorithm."""
# TODO: using the clustering reveals that some nodes and used_tasks are
# identical. This should be improved to improved speed.

# Create the first node
        self.nodes = []
        graph = []
        tasks_done = []
        tasks_ready = []
        cost = 0
        if self.bound_estimator_type == 'adams':
# When the grid is regular min_bound is given by P_x + P_y + N_tasks - 4 A good
# approximation is to replace P_x + P_y by 2*sqrt(P) and N_tasks by 
# max_proc(sum_i delta_t_i)
            tasks_per_proc = [0 for i in range(self.n_processors)]
            for task in self.tasks:
                tasks_per_proc[task.subdomain_id] += task.delta_t
            min_bound = 2*np.sqrt(self.n_processors) + max(tasks_per_proc) - 4
        else:
            self.Compute_b_level()
            min_bound = self.Compute_bound_estimator(tasks_ready, tasks_done)
# Build the tasks_ready list. This list constains all the tasks that are ready
# to be executed. At first, these tasks are the one that do not require
# another task.
        for task in self.tasks:
            if len(task.required_tasks) == 0:
                tasks_ready.append(task)
        first_node = NODE(graph, tasks_done, tasks_ready, cost, min_bound)
        self.nodes.append(first_node)

        done = False
        while not done:
            counter = 0
            max_graph_length = -1
            for node in self.nodes:
                if len(node.graph) > max_graph_length:
                    max_graph_length = len(node.graph)
                    min_bound_deepest = node.min_bound
                    nodes_to_explore = [counter]
                    proc_used = [len(self.Compute_subdomains_list(self.nodes[counter]))]
                elif len(node.graph) == max_graph_length:
                    if node.min_bound < min_bound_deepest:
                        min_bound_deepest = node.min_bound
                        nodes_to_explore = [counter]
                        proc_used = [len(self.Compute_subdomains_list(self.nodes[counter]))]
                    elif node.min_bound == min_bound_deepest:
                        nodes_to_explore.append(counter)
                        proc_used.append(len(self.Compute_subdomains_list(self.nodes[counter])))
                counter += 1

# Only keep the nodes that are using the most processors. max_proc_used can be
# choosen lower than max(proc_used) if there are convergence problems.
            max_proc_used = max(proc_used)
            nodes_to_pop = []
            for pos in range(len(nodes_to_explore)):
                if proc_used[pos] < max_proc_used:
                    nodes_to_pop.append(pos)
            for pos in reversed(nodes_to_pop):
                nodes_to_explore.pop(pos)

# Cluster the tasks and then, only explore the centroids of theses clusters.
            print()
            print()
            print("number of nodes to explore",len(nodes_to_explore))
            nodes_to_explore = self.Compute_nodes_centroids(nodes_to_explore)
            print("number of reduced nodes to explore",len(nodes_to_explore))

            for pos in nodes_to_explore:
                subdomains_list = self.Compute_subdomains_list(self.nodes[pos])
                print()
                print('graph length', len(self.nodes[pos].graph))
                print('subdomain list', subdomains_list)
                print('number of nodes', len(self.nodes))
                print('min_bound', self.nodes[pos].min_bound)

# The combination increases too fast so we first try to run as many processors
# as possible. If we can add a node using as many processors as possible, we
# exit the loop.
                for i in range(len(subdomains_list), 0, -1):
# itertools create a generator, i.e., a list of iterators that can be used only
# once. If we want to reuse the list, we need to recall itertools
                    proc_comb = itertools.combinations(subdomains_list, i)
                    node_added = False
                    print('processors used', i)
                    for used_procs in proc_comb:
# Get ride of equivalent tasks in tasks_ready
                        reduced_tasks_ready = self.Reduce_tasks_ready(
                                self.nodes[pos].tasks_ready)
# Compute the number of nodes to create
                        n_tasks_ready_per_proc = []
                        for tmp in reduced_tasks_ready:
                            n_tasks_ready_per_proc.append(len(tmp))
                        n_tasks_ready = 1
                        for tmp in n_tasks_ready_per_proc:
                            n_tasks_ready *= tmp
                        print('number of tasks ready', n_tasks_ready)
                        current_pos = [0 for i in range(len(n_tasks_ready_per_proc))]
# Create all the used_tasks possible so that they can be clustered later on.
                        full_used_tasks = []
                        for j in range(n_tasks_ready):
                            used_tasks = []
                            for proc in range(len(n_tasks_ready_per_proc)):
                                used_tasks.append(reduced_tasks_ready[proc][current_pos[proc]])
                            full_used_tasks.append(used_tasks)
                            current_pos[0] += 1
                            for proc in range(len(n_tasks_ready_per_proc)):
                                if current_pos[proc] == n_tasks_ready_per_proc[proc] and\
                                        proc+1 < len(n_tasks_ready_per_proc):
                                    current_pos[proc] = 0
                                    current_pos[proc+1] += 1

# Reduce the number of full_used_tasks by clustering them
                        reduced_used_tasks = self.Compute_used_tasks_centroids(
                                full_used_tasks)
                        print('reduced number of tasks ready',
                                len(reduced_used_tasks))
                        for used_tasks in reduced_used_tasks:
                            node = self.Create_node(used_procs, used_tasks,
                                                self.nodes[pos])
                            if node.min_bound < self.best_time:
                                self.nodes.append(node)
                                node_added = True
                    if node_added is True:
                        break

# The parent nodes can be deleted. We deleted the elements from the end first so
# the other indices stay valid
            for pos in reversed(nodes_to_explore):
                self.nodes.pop(pos)

# Removes all the nodes but the ones at the deepest level
            nodes_to_remove = []
            for pos in range(len(self.nodes)):
                if len(self.nodes[pos].graph) <= max_graph_length:
                    nodes_to_remove.append(pos)
            for pos in reversed(nodes_to_remove):
                self.nodes.pop(pos)

            done = False
            if (len(self.nodes) == 0):
                print('All the possibilities have been exhausted')
                done = True
            for node in self.nodes:
                if len(node.tasks_ready) == 0:
# Stop if the solution is better than the best_time
                    if node.cost <= self.best_time:
                        done = True
                        self.nodes[0] = node
                        print('Output result.')
                        print('Node cost: ', node.cost)
                        if self.output is True:
                            self.Output_results()

                        break
# Print the temporary solution
                    else:
                        if tmp_node is not node:
                            print('Output intermediate result because BB did',
                                  'not converge.')
                            print('Node cost: ', node.cost)

#----------------------------------------------------------------------------#

    def Output_results(self):
        """Output intermediate results when Branch-and Bound does not
        converge."""
        print('Minimum cost: '+str(self.nodes[0].cost))
        for level in self.nodes[0].graph:
            print('-------------')
            for task in level:
                print('subdomain_id: '+str(task.subdomain_id))
                print('Id: '+str(task.task_id))
