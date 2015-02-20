# Python code
# Author: Bruno Turcksin
# Date: 2015-01-25 16:31:13.044256

#----------------------------------------------------------------------------#
## Classes BRANCH_AND_BOUND and NODE                                        ##
#----------------------------------------------------------------------------#

"""This module contains the class Branch_and_bound that solve the SOP using
banch-and-bound algorithm and the class Node that represent the different node
of the graph created by the branch-and-bound algorithm."""

import itertools
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

    def __init__(self, tasks, n_processors, best_time):

        super().__init__()
        self.tasks = tasks
        self.n_tasks = len(tasks)
        self.n_processors = n_processors
        self.best_time = best_time
        self.current_directions = []

#----------------------------------------------------------------------------#

    def Create_node(self, given_procs, given_tasks, parent_node):
        """Create a new node given a processor combination, a task combination,
        and a node."""

        n_tasks = len(given_tasks)
        used_procs = set()

        for task in given_tasks:
            if task.subdomain_id in given_procs and\
               task.subdomain_id not in used_procs:
                n_tasks -= 1
                used_procs.add(task.subdomain_id)
# Check that the processor and the task combination exist
        if n_tasks == 0:
            max_delta_t = 0
            new_graph = parent_node.graph[:]
            current_level = []
            new_tasks_done = parent_node.tasks_done.copy()
            new_tasks_ready = parent_node.tasks_ready[:]
            for task in given_tasks:
                current_level.append(task)
                new_tasks_done.add(task)
                new_tasks_ready.remove(task)
                if max_delta_t < task.delta_t:
                    max_delta_t = task.delta_t
            new_graph.append(current_level)
# Add new tasks to the tasks_ready list
            for task in self.tasks:
                if task not in new_tasks_done:
                    ready = True
                    for required_task in task.required_tasks:
                        if required_task not in new_tasks_done:
                            ready = False
                            break
                    if ready is True:
                        if task not in new_tasks_ready:
                            new_tasks_ready.append(task)
# The new cost is the cost of parent_node plus max_delta_t
            new_cost = parent_node.cost + max_delta_t
#TODO min_bound = max n_tasks left for a processor + sqrt(procs that have
#     executed a task but have tasks left)/2 + sqrt(procs that haven't executed
#     a task) - 4
# The min_bound is given by the largest number of tasks left for one
# processors
            n_tasks_proc = [0 for i in range(self.n_processors)]
            for task in self.tasks:
                if task not in new_tasks_done:
                    n_tasks_proc[task.subdomain_id] += 1
            new_min_bound = new_cost + max(n_tasks_proc)
# Create the new node
            node = Node(new_graph, new_tasks_done, new_tasks_ready, new_cost,
                        new_min_bound)
        else:
# The node given by the processor combination and the task combination does
# not exist. cost and min_bound are set to 2*len(self.tasks)
            node = NODE(parent_node.graph, parent_node.tasks_done,
                        parent_node.tasks_ready, 2*len(self.tasks),
                        2*len(self.tasks))

        return node

#----------------------------------------------------------------------------#

    def Solve(self):
        """Solve the SOP using branch-and-bound algorithm."""

# Create the first node
        self.nodes = []
        graph = []
        tasks_done = set()
        tasks_ready = []
        cost = 0
# When the grid is regular min_bound is given by P_x + P_y + N_tasks - 4 A good
# approximation is to replace P_x + P_y by 2*sqrt(P) and N_tasks by 
# max_proc(sum_i delta_t_i)
        tasks_per_proc = [0 for i in range(self.n_processors)]
        for task in self.tasks:
          tasks_per_proc[task.subdomain_id] += task.delta_t
        min_bound = 2*int(np.sqrt(self.n_processors)) + max(tasks_per_proc) - 4
# Build the tasks_ready list. This list constains all the tasks that are ready
# to be executed. At first, these tasks are the one that do not require
# another task.
        for task in self.tasks:
            if task.n_required_tasks == 0:
                tasks_ready.append(task)
        first_node = NODE(graph, tasks_done, tasks_ready, cost, min_bound)
        self.nodes.append(first_node)

        subdomains_list = list(range(0, self.n_processors))
        done = False
        best_first = False
        while not done:
            lowest_bound = self.best_time
            n_tasks_done = 0
            pos = 0
            counter = 0

            max_level = 0
            for node in self.nodes:
# Best-first search
                if node.min_bound < lowest_bound:
                    pos = counter
                    lowest_bound = node.min_bound
                counter += 1

            for i in range(1, self.n_processors+1):
# itertools create a generator, i.e., a list of iterators that can be used only
# once. If we want to reuse the list, we need to recall itertools
                proc_comb = itertools.combinations(subdomains_list, i)
                for used_procs in proc_comb:
                    task_comb = itertools.combinations(self.nodes[pos].tasks_ready,
                                                       i)
                    for used_tasks in task_comb:
                        node = self.Create_node(used_procs, used_tasks,
                                                self.nodes[pos])
                        self.nodes.append(node)

# Every possible children from the node was created, thus the node itself can
# be deleted.
            if self.nodes[pos].cost != self.nodes[pos].min_bound:
                self.nodes.pop(pos)
# Delete all the nodes that have min_bound larger than the best_time
            self.nodes = [node for node in self.nodes
                          if node.min_bound <= self.best_time]

            done = False
            for node in self.nodes:
                if len(node.tasks_ready) == 0:
# Stop if the solution is better than the best_time
                    if node.cost <= self.best_time:
                        done = True
                        self.nodes[0] = node
                        print('Output result.')
                        print('Node cost: ', node.cost)
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
                print('Id: '+str(task.ID))
