# Python code
# Author: Bruno Turcksin
# Date: 2014-09-03 16:51:38.125245

#----------------------------------------------------------------------------#
## Module CAP_PFB                                                           ##
#----------------------------------------------------------------------------#

"""In this module, we implement all the functions and the class needed to
compute a CAP-PFB ordering."""

import itertools
import copy
import LOCAL_SEARCH
import TASK


def build_tasks(n_x, n_y, n_p, d_t, sn, cells_to_refine):
    """Build the tasks used in the rest of the program. It creates n_x x n_y x 4
    tasks, corresponding to an uniform mesh with n_x by n_y cells and a S2
    discretization. n_p is map between the position of the cells and the
    processor ids. d_t is a map between the position of the cells and the time
    the tasks associated will need to be executed. sn is the order of the
    discrete ordinate methods. cells_to_refine contains the position (i,j) of
    the cells that must be refined."""

    tasks = []
    n_dir = {2: 1, 4: 3, 8: 10, 16: 36}
    task_id_offset = 0
    tasks_to_refine = []
    idir = 0
    for dir in range(n_dir[sn]):
# First direction (x: +, y: +)
        for pos in itertools.product(range(n_x), range(n_y)):
            i, j = pos
            task_id = i + j*n_x + task_id_offset
            subdomain_id = n_p[i][j]
            delta_t = d_t[i][j]
            required_tasks = []
            waiting_tasks = []
            if i >= 1:
                required_tasks.append((i-1)+j*n_x+task_id_offset)
            if j >= 1:
                required_tasks.append(i+(j-1)*n_x+task_id_offset)
            if i < (n_x-1):
                waiting_tasks.append((i+1)+j*n_x+task_id_offset)
            if j < (n_y-1):
                waiting_tasks.append(i+(j+1)*n_x+task_id_offset)
            task = TASK.TASK(subdomain_id, task_id, required_tasks,
                             waiting_tasks, delta_t, pos, idir)
            tasks.append(task)
            if [i, j] in cells_to_refine:
                tasks_to_refine.append([task, 0])
        idir += 1
        task_id_offset += n_x*n_y

# Second direction (x: +, y: -)
        for pos in itertools.product(range(n_x), range(n_y-1, -1, -1)):
            i, j = pos
            task_id = i + j*n_x + task_id_offset
            subdomain_id = n_p[i][j]
            delta_t = d_t[i][j]
            required_tasks = []
            waiting_tasks = []
            if i >= 1:
                required_tasks.append((i-1)+j*n_x+task_id_offset)
            if j < (n_y-1):
                required_tasks.append(i+(j+1)*n_x+task_id_offset)
            if i < (n_x-1):
                waiting_tasks.append((i+1)+j*n_x+task_id_offset)
            if j >= 1:
                waiting_tasks.append(i+(j-1)*n_x+task_id_offset)
            task = TASK.TASK(subdomain_id, task_id, required_tasks,
                             waiting_tasks, delta_t, pos, idir)
            tasks.append(task)
            if [i, j] in cells_to_refine:
                tasks_to_refine.append([task, 1])
        idir += 1
        task_id_offset += n_x*n_y

# Third direction (x: -, y: +)
        for pos in itertools.product(range(n_x-1, -1, -1), range(n_y)):
            i, j = pos
            task_id = i + j*n_x + task_id_offset
            subdomain_id = n_p[i][j]
            delta_t = d_t[i][j]
            required_tasks = []
            waiting_tasks = []
            if i < (n_x-1):
                required_tasks.append((i+1)+j*n_x+task_id_offset)
            if j >= 1:
                required_tasks.append(i+(j-1)*n_x+task_id_offset)
            if i >= 1:
                waiting_tasks.append((i-1)+j*n_x+task_id_offset)
            if j < (n_y-1):
                waiting_tasks.append(i+(j+1)*n_x+task_id_offset)
            task = TASK.TASK(subdomain_id, task_id, required_tasks,
                             waiting_tasks, delta_t, pos, idir)
            tasks.append(task)
            if [i, j] in cells_to_refine:
                tasks_to_refine.append([task, 2])
        idir += 1
        task_id_offset += n_x*n_y


# Fourth direction (x: -, y: -)
        for pos in itertools.product(range(n_x-1, -1, -1), range(n_y-1, -1, -1)):
            i, j = pos
            task_id = i + j*n_x + task_id_offset
            subdomain_id = n_p[i][j]
            delta_t = d_t[i][j]
            required_tasks = []
            waiting_tasks = []
            if i < (n_x-1):
                required_tasks.append((i+1)+j*n_x+task_id_offset)
            if j < (n_y-1):
                required_tasks.append(i+(j+1)*n_x+task_id_offset)
            if i >= 1:
                waiting_tasks.append((i-1)+j*n_x+task_id_offset)
            if j >= 1:
                waiting_tasks.append(i+(j-1)*n_x+task_id_offset)
            task = TASK.TASK(subdomain_id, task_id, required_tasks,
                             waiting_tasks, delta_t, pos, idir)
            tasks.append(task)
            if [i, j] in cells_to_refine:
                tasks_to_refine.append([task, 3])
        idir += 1
        task_id_offset += n_x*n_y

# Refine the cells necessary.
    for task_to_refine in tasks_to_refine:
        task = task_to_refine[0]
        dir = task_to_refine[1]
        subdomain_id = task.subdomain_id
        task_id = task.task_id
        if dir == 0:
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                   task.required_tasks, [task_id_offset+1,
                                     task_id_offset+2], int(task.delta_t/4),
                                   task.pos, task.dir))
            # On the corner no task may be required or waiting
            if task.required_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+1,
                                       [task_id_offset+0], [task_id_offset+3,
                                         min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                       [task_id_offset+0], [task_id_offset+3,
                                         max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
            elif task.waiting_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+1,
                                      [task_id_offset+0,
                                        min(task.required_tasks)],
                                      [task_id_offset+3], int(task.delta_t/4),
                                      task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                       [task_id_offset+0,
                                         max(task.required_tasks)],
                                       [task_id_offset+3], int(task.delta_t/4),
                                       task.pos, task.dir))
            else:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+1,
                                       [task_id_offset+0,
                                         min(task.required_tasks)],
                                       [task_id_offset+3,
                                         min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos,
                                       task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                       [task_id_offset+0,
                                         max(task.required_tasks)],
                                       [task_id_offset+3,
                                         max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos,
                                       task.dir))
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                   [task_id_offset+1, task_id_offset+2],
                                   task.waiting_tasks, int(task.delta_t/4),
                                   task.pos, task.dir))
            for other_task in tasks:
                if task_id in other_task.required_tasks:
                    other_task.required_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.required_tasks.append(task_id_offset+2)
                        other_task.required_tasks.append(task_id_offset+3)
                    else:
                        other_task.required_tasks.append(task_id_offset+1)
                        other_task.required_tasks.append(task_id_offset+3)
                if task_id in other_task.waiting_tasks:
                    other_task.waiting_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.waiting_tasks.append(task_id_offset+0)
                        other_task.waiting_tasks.append(task_id_offset+1)
                    else:
                        other_task.waiting_tasks.append(task_id_offset+0)
                        other_task.waiting_tasks.append(task_id_offset+2)
        if dir == 1:
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                   task.required_tasks, [task_id_offset+0,
                                     task_id_offset+3],
                                   int(task.delta_t/4), task.pos, task.dir))
            # On the corner no task may be required or waiting
            if task.required_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                       [task_id_offset+2], [task_id_offset+1,
                                         min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                       [task_id_offset+2], [task_id_offset+1,
                                         max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
            elif task.waiting_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                       [task_id_offset+2,
                                         min(task.required_tasks)],
                                       [task_id_offset+1], int(task.delta_t/4),
                                       task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                       [task_id_offset+2,
                                         max(task.required_tasks)],
                                       [task_id_offset+1], int(task.delta_t/4),
                                       task.pos, task.dir))
            else:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                       [task_id_offset+2,
                                         min(task.required_tasks)],
                                       [task_id_offset+1,
                                         min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                       [task_id_offset+2,
                                         max(task.required_tasks)],
                                       [task_id_offset+1,
                                         max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+1,
                                   [task_id_offset+0, task_id_offset+3],
                                   task.waiting_tasks, int(task.delta_t/4),
                                   task.pos, task.dir))
            for other_task in tasks:
                if task_id in other_task.required_tasks:
                    other_task.required_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.required_tasks.append(task_id_offset+0)
                        other_task.required_tasks.append(task_id_offset+1)
                    else:
                        other_task.required_tasks.append(task_id_offset+1)
                        other_task.required_tasks.append(task_id_offset+3)
                if task_id in other_task.waiting_tasks:
                    other_task.waiting_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.waiting_tasks.append(task_id_offset+2)
                        other_task.waiting_tasks.append(task_id_offset+3)
                    else:
                        other_task.waiting_tasks.append(task_id_offset+0)
                        other_task.waiting_tasks.append(task_id_offset+2)
        if dir == 2:
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+1, task.required_tasks,
                                   [task_id_offset+0, task_id_offset+3],
                                   int(task.delta_t/4), task.pos, task.dir))
            # On the corner no task may be required or waiting
            if task.required_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                       [task_id_offset+1], [task_id_offset+2,
                                         min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                       [task_id_offset+1], [task_id_offset+2,
                                         max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
            elif task.waiting_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                       [task_id_offset+1, min(task.required_tasks)],
                                       [task_id_offset+2], int(task.delta_t/4),
                                       task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                       [task_id_offset+1,
                                         max(task.required_tasks)],
                                       [task_id_offset+2], int(task.delta_t/4),
                                       task.pos, task.dir))
            else:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                       [task_id_offset+1,
                                         min(task.required_tasks)],
                                       [task_id_offset+2,
                                         min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos,
                                       task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                       [task_id_offset+1,
                                         max(task.required_tasks)],
                                       [task_id_offset+2,
                                         max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos, task.dir))
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                   [task_id_offset+0, task_id_offset+3],
                                   task.waiting_tasks, int(task.delta_t/4),
                                   task.pos, task.dir))
            for other_task in tasks:
                if task_id in other_task.required_tasks:
                    other_task.required_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.required_tasks.append(task_id_offset+2)
                        other_task.required_tasks.append(task_id_offset+3)
                    else:
                        other_task.required_tasks.append(task_id_offset+0)
                        other_task.required_tasks.append(task_id_offset+2)
                if task_id in other_task.waiting_tasks:
                    other_task.waiting_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.waiting_tasks.append(task_id_offset+0)
                        other_task.waiting_tasks.append(task_id_offset+1)
                    else:
                        other_task.waiting_tasks.append(task_id_offset+1)
                        other_task.waiting_tasks.append(task_id_offset+3)
        if dir == 3:
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+3,
                                   task.required_tasks, [task_id_offset+1,
                                     task_id_offset+2], int(task.delta_t/4),
                                   task.pos, task.dir))
            # On the corner no task may be required or waiting
            if task.required_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+1,
                                       [task_id_offset+3], [task_id_offset+0,
                                         min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos,
                                       task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                       [task_id_offset+3], [task_id_offset+0,
                                         max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos,
                                       task.dir))
            elif task.waiting_tasks == []:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+1,
                                       [task_id_offset+3,
                                        min(task.required_tasks)],
                                       [task_id_offset+0], int(task.delta_t/4),
                                       task.pos, task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                       [task_id_offset+3,
                                        max(task.required_tasks)],
                                       [task_id_offset+0], int(task.delta_t/4),
                                       task.pos, task.dir))
            else:
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+1,
                                       [task_id_offset+3,
                                        min(task.required_tasks)],
                                       [task_id_offset+0,
                                        min(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos,
                                       task.dir))
                tasks.append(TASK.TASK(subdomain_id, task_id_offset+2,
                                       [task_id_offset+3,
                                        max(task.required_tasks)],
                                       [task_id_offset+0,
                                        max(task.waiting_tasks)],
                                       int(task.delta_t/4), task.pos,
                                       task.dir))
            tasks.append(TASK.TASK(subdomain_id, task_id_offset+0,
                                   [task_id_offset+1, task_id_offset+2],
                                   task.waiting_tasks, int(task.delta_t/4),
                                   task.pos, task.dir))
            for other_task in tasks:
                if task_id in other_task.required_tasks:
                    other_task.required_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.required_tasks.append(task_id_offset+0)
                        other_task.required_tasks.append(task_id_offset+1)
                    else:
                        other_task.required_tasks.append(task_id_offset+0)
                        other_task.required_tasks.append(task_id_offset+2)
                if task_id in other_task.waiting_tasks:
                    other_task.waiting_tasks.remove(task_id)
                    if task.pos[0] == other_task.pos[0]:
                        other_task.waiting_tasks.append(task_id_offset+2)
                        other_task.waiting_tasks.append(task_id_offset+3)
                    else:
                        other_task.waiting_tasks.append(task_id_offset+1)
                        other_task.waiting_tasks.append(task_id_offset+3)
        task_id_offset += 4
        tasks.remove(task)

    return tasks

#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#


class CAP_PFB(object):
    """This class creates the CAP-PFB sweep ordering."""

    def __init__(self, tasks, max_it, tol, n_p, forbid_oscillations,
                 do_local_search):
        """The constructor requires a list of tasks to order, a maximum number of
        iterations, a tolerance, and the processors arrangement."""

        super().__init__()
        self.tasks = tasks
        self.max_it = max_it
        self.tol = tol
        self.n_p = n_p
        self.forbid_oscillations = forbid_oscillations
        self.do_local_search = do_local_search
        self.task_id_map = dict()
        for task in self.tasks:
            self.task_id_map[task.task_id] = task

#----------------------------------------------------------------------------#

    def Create_initial_scheduling(self):
        """Create a random initial scheduling."""

        self.schedule = []
        self.schedule_id_map = dict()
        tasks_done = set()
        proc_end_time = [0]*(max(max(self.n_p))+1)
        pos = 0
        while len(tasks_done) != len(self.tasks):
            for task in self.tasks:
# Check that the task is ready to be executed
                if task.task_id not in tasks_done:
                    ready = True
                    for required_task in task.required_tasks:
                        if self.task_id_map[required_task].task_id not in tasks_done:
                            ready = False
                            break
                    if ready == True:
                        starting_time = 0
# Search for the time when the latest required task is done
                        for scheduled_task in self.schedule:
                            scheduled_task_id = scheduled_task[0].task_id
                            if scheduled_task_id in task.required_tasks:
                                if scheduled_task[2] > starting_time:
                                    starting_time = scheduled_task[2]
                        starting_time = max(starting_time,
                                            proc_end_time[task.subdomain_id])
                        self.schedule.append([task, starting_time,
                                              starting_time+task.delta_t])
                        self.schedule_id_map[task.task_id] = pos
                        tasks_done.add(task.task_id)
                        proc_end_time[task.subdomain_id] = starting_time+task.delta_t
                        pos += 1

        self.start_time = 100000
        self.end_time = 0
        for task in self.schedule:
            if task[1] < self.start_time:
                self.start_time = task[1]
            if task[2] > self.end_time:
                self.end_time = task[2]
        print("Initial scheduling")
        print("Start time", self.start_time)
        print("End time", self.end_time)
        print("Execution time", self.end_time-self.start_time)

#----------------------------------------------------------------------------#

    def Create_serial_initial_scheduling(self):
        """Compute an initial scheduling. This scheduling is similar to a serial
        scheduling, i.e. only one task is executed at any given time."""

        self.schedule = []
        self.schedule_id_map = dict()
        tasks_done = set()
        pos = 0
        current_dir = 0
        starting_time = 0
        while len(tasks_done) != len(self.tasks):
            update_dir = True
            for task in self.tasks:
# Check that the task is ready to be executed
                if task.task_id not in tasks_done and task.dir == current_dir:
                    update_dir = False
                    ready = True
                    for required_task in task.required_tasks:
                        if self.task_id_map[required_task].task_id not in tasks_done:
                            ready = False
                            break
                    if ready is True:
                        self.schedule.append([task, starting_time,
                                              starting_time+task.delta_t])
                        self.schedule_id_map[task.task_id] = pos
                        tasks_done.add(task.task_id)
                        starting_time += task.delta_t
                        pos += 1
            if update_dir is True:
                current_dir += 1

        self.start_time = 100000
        self.end_time = 0
        for task in self.schedule:
            if task[1] < self.start_time:
                self.start_time = task[1]
            if task[2] > self.end_time:
                self.end_time = task[2]
        print("Initial scheduling")
        print("Start time", self.start_time)
        print("End time", self.end_time)
        print("Execution time", self.end_time-self.start_time)

#----------------------------------------------------------------------------#

    def Create_dfds_local_initial_scheduling(self):
        """Compute an initial scheduling using interleaved DFDS. To compute the
        b-level of each task, a serial initial scheduling is computed first."""


# Use the serial scheduling to sort the task
        self.Create_serial_initial_scheduling()

# Compute the ranks
        ranks = [0]*len(self.tasks)
        b_level = [0]*len(self.tasks)
# Constant must be choosen such that it is larger than the largest b-level but
# it is unknwon before the computation. Thus, take an arbitrary large number.
        constant = 10000
        done = False
        for i in range(len(self.schedule)-1, -1, -1):
            task = self.schedule[i][0]
            for waiting_task_id in task.waiting_tasks:
                wait_sch_task = self.schedule_id_map[waiting_task_id]
                if self.schedule[wait_sch_task][0].subdomain_id != task.subdomain_id:
                    candidate_rank = b_level[wait_sch_task] + constant
                else:
                    candidate_rank = ranks[wait_sch_task] - 1
                if candidate_rank > ranks[i]:
                    ranks[i] = candidate_rank
                if (b_level[wait_sch_task]+1) > b_level[i]:
                    b_level[i] = b_level[wait_sch_task]+1

# Compute the new scheduling
        new_schedule = []
# Set is a mutable so cannot use [set()]*n_procs
        n_procs = max(max(self.n_p))+1
        proc_start_time = [set() for i in range(n_procs)]
        end_time = int(self.end_time/(n_procs/10))
        for i in range(n_procs):
            for j in range(0, end_time+1):
                proc_start_time[i].add(j)
        tasks_done = set()
        pos_list = []
        for task in self.schedule:
            if task[0].required_tasks == []:
                pos_list.append(self.schedule_id_map[task[0].task_id])
        for i in range(len(self.schedule)):
            rejected_tasks = set()
            found = False
            while found is False:
                found = True
                max_rank = -1
                pos = -1
                k = 0
                for j in pos_list:
                    if ranks[j] > max_rank and j not in rejected_tasks:
                        max_rank = ranks[j]
                        pos = j
                        pos_to_pop = k
                    k += 1
                tmp_task = self.schedule[pos][0]
                for req_task_id in tmp_task.required_tasks:
                    if req_task_id not in tasks_done:
                        found = False
                        rejected_tasks.add(pos)
                        break

            task = self.schedule[pos][0]
            pos_list.pop(pos_to_pop)
            for waiting_task_id in task.waiting_tasks:
                wait_sch_task = self.schedule_id_map[waiting_task_id]
                if waiting_task_id not in tasks_done and\
                   wait_sch_task not in pos_list:
                    pos_list.append(wait_sch_task)
            starting_time = 0
            for scheduled_task in new_schedule:
                scheduled_task_id = scheduled_task[0].task_id
                if scheduled_task_id in task.required_tasks:
                    if scheduled_task[2] > starting_time:
                        starting_time = scheduled_task[2]
            candidate_t = max(starting_time,
                              min(proc_start_time[task.subdomain_id]))
            enough_free_t = False
            while enough_free_t is False:
                enough_free_t = True
                for j in range(task.delta_t):
                    if (candidate_t+j) not in proc_start_time[task.subdomain_id]:
                        enough_free_t = False
                        candidate_t += 1
                        break
            end_time = candidate_t+task.delta_t
            new_schedule.append([task, candidate_t, end_time])
            for j in range(task.delta_t):
                proc_start_time[task.subdomain_id].remove(candidate_t+j)
            tasks_done.add(task.task_id)
        self.schedule = new_schedule

        self.Update_schedule_id_map()

        self.start_time = 100000
        self.end_time = 0
        for task in self.schedule:
            if task[1] < self.start_time:
                self.start_time = task[1]
            if task[2] > self.end_time:
                self.end_time = task[2]
        print("Initial scheduling")
        print("Start time", self.start_time)
        print("End time", self.end_time)
        print("Execution time", self.end_time-self.start_time)

#----------------------------------------------------------------------------#

    def Run(self, init_sched):
        """Driver for CAP-PFB."""

        if init_sched == 'serial':
            print('Create serial initial scheduling')
            self.Create_serial_initial_scheduling()
        elif init_sched == 'DFDS':
            print('Create DFDS initial scheduling')
            self.Create_dfds_local_initial_scheduling()
        else:
            print('Create initial scheduling')
            self.Create_initial_scheduling()
        tmp_schedule = copy.deepcopy(self.schedule)
        tmp_start = self.start_time
        tmp_end = self.end_time
        for iter in range(self.max_it):
            old_start = self.start_time
            old_end = self.end_time
            old_execution = old_end-old_start
            self.Backward_iteration()
            if self.forbid_oscillations is True:
                execution = self.end_time-self.start_time
                if execution > old_execution:
                    self.schedule = copy.deepcopy(tmp_schedule)
                    self.schedule.reverse()
                    execution = old_execution
                    self.start_time = tmp_start
                    self.end_time = tmp_end
                    self.Update_schedule_id_map()
                else:
                    tmp_schedule = copy.deepcopy(self.schedule)
                    tmp_schedule.reverse()
                    tmp_start = self.start_time
                    tmp_end = self.end_time

            self.Forward_iteration()
            if self.forbid_oscillations is True:
                new_execution = self.end_time-self.start_time
                if new_execution > execution:
                    self.schedule = copy.deepcopy(tmp_schedule)
                    new_execution = execution
                    self.start_time = tmp_start
                    self.end_time = tmp_end
                    self.Update_schedule_id_map()
                else:
                    tmp_schedule = copy.deepcopy(self.schedule)
                    tmp_start = self.start_time
                    tmp_end = self.end_time

                if (old_execution-new_execution) < self.tol:
                    break

        if self.do_local_search is True:
            n_samples = 500
            n_procs = max(max(self.n_p))+1
            local_search = LOCAL_SEARCH.LOCAL_SEARCH(self.schedule, n_procs,
                                                     self.start_time,
                                                     self.end_time, n_samples)
            local_search.Execute_search()
            self.schedule = []
            for proc in range(n_procs):
                self.schedule.append(local_search.schedule[proc])

#----------------------------------------------------------------------------#

    def Update_schedule_id_map(self):
        """"Update the schedule_id_map."""

# Update schedule_id_map
        pos = 0
        for task in self.schedule:
            self.schedule_id_map[task[0].task_id] = pos
            pos += 1

#----------------------------------------------------------------------------#

    def Backward_iteration(self):
        """Execute backward iteration."""

        ranks = self.Compute_backward_ranks()
        self.Backward_sort(ranks)
        self.Backward_scheduling()
        print("Backward scheduling")
        print("Start time", self.start_time)
        print("End time", self.end_time)
        print("Execution time", self.end_time-self.start_time)

#----------------------------------------------------------------------------#

    def Compute_backward_ranks(self):
        """Compute the backward rank of each vertex (task)."""

        ranks = [-1]*len(self.schedule)
        done = False
        pos = 0
        for task in self.schedule:
            for required_task_id in task[0].required_tasks:
                req_sch_task = self.schedule_id_map[required_task_id]
                if self.schedule[req_sch_task][0].subdomain_id != task[0].subdomain_id:
                    candidate_rank = self.schedule[req_sch_task][2]
                else:
                    candidate_rank = ranks[req_sch_task]
                if ranks[pos] < candidate_rank:
                    ranks[pos] = candidate_rank
            pos += 1

        return ranks

#----------------------------------------------------------------------------#

    def Backward_sort(self, ranks):
        """Sort the tasks in the backward iteration."""

        for new_pos in range(len(ranks)):
            max_rank = max(ranks[new_pos:])
            old_pos = new_pos+ranks[new_pos:].index(max_rank)
# Put the task earlier in the list
            ranks[new_pos], ranks[old_pos] = ranks[old_pos], ranks[new_pos]
            self.schedule[new_pos], self.schedule[old_pos] =\
                self.schedule[old_pos], self.schedule[new_pos]

# Sort tasks with the same rank
        pos = 0
        while pos < len(ranks):
            end_pos = pos + 1
            while ranks[pos] == ranks[end_pos]:
                end_pos += 1
                if end_pos == len(ranks):
                    break
            for new_pos in range(pos, end_pos):
                latest_time = -1
                for i in range(new_pos, end_pos):
                    if self.schedule[i][2] > latest_time:
                        latest_time = self.schedule[i][2]
                for old_pos in range(new_pos, end_pos):
                    if self.schedule[old_pos][2] == latest_time:
                        break
                self.schedule[new_pos], self.schedule[old_pos] =\
                    self.schedule[old_pos], self.schedule[new_pos]
            pos = end_pos

# Update schedule_id_map
        self.Update_schedule_id_map()

#----------------------------------------------------------------------------#

    def Backward_scheduling(self):
        """Compute a scheduling that starts the tasks as late as possible."""

        n_procs = max(max(self.n_p))+1
# Simulate multiple processsors
        part_schedule = []
        for i in range(n_procs):
            part_schedule.append([])
        for task in self.schedule:
            part_schedule[task[0].subdomain_id].append(task)
# Set is a mutable so cannot use [set()]*n_procs
        proc_end_time = [set() for i in range(n_procs)]
        end_time = self.end_time+1
        start_time = max(self.start_time-100,
                         int(self.end_time-len(self.schedule)/(n_procs/10)))
        if start_time < 0:
            end_time -= start_time
            start_time = 0
            self.end_time = end_time-1
        for i in range(n_procs):
            for j in range(start_time, end_time):
                proc_end_time[i].add(j)
# Flag to know if the time of the waiting task has been update
        updated_time = [False]*len(self.schedule)

        n_tasks = len(self.schedule)
        n_tasks_done = 0
        old_n_tasks_done = -1
        while n_tasks_done != n_tasks:
            assert old_n_tasks_done != n_tasks_done,\
                'No improvement in Backward scheduling.'
            old_n_tasks_done = n_tasks_done
            for i in range(n_procs):
                if len(part_schedule[i]) != 0:
                    task = part_schedule[i][0]
                    ready = True
                    t = self.end_time
                    if task[0].waiting_tasks != []:
                        for waiting_task_id in task[0].waiting_tasks:
                            waiting_task = self.schedule[
                                self.schedule_id_map[waiting_task_id]]
                            ready = updated_time[self.schedule_id_map[
                                waiting_task_id]]
                            if ready is False:
                                break
                            if t > waiting_task[1]:
                                t = waiting_task[1]
                    if ready is True:
                        candidate_t = min(t, max(proc_end_time[i]))
                        enough_free_t = False
                        while enough_free_t is False:
                            enough_free_t = True
                            for j in range(task[0].delta_t):
                                if (candidate_t-j) not in proc_end_time[i]:
                                    enough_free_t = False
                                    candidate_t -= 1
                                    break
                        task[2] = candidate_t
                        task[1] = task[2]-task[0].delta_t
                        for j in range(task[0].delta_t):
                            proc_end_time[i].remove(task[2]-j)
                        self.schedule[self.schedule_id_map[task[0].task_id]] = task
                        updated_time[self.schedule_id_map[task[0].task_id]] = True
                        part_schedule[i].pop(0)
                        n_tasks_done += 1

# Update starting time
        self.start_time = self.end_time
        for task in self.schedule:
            if self.start_time > task[1]:
                self.start_time = task[1]

#----------------------------------------------------------------------------#

    def Forward_iteration(self):
        """Execute forkward iteration."""

        ranks = self.Compute_forward_ranks()
        self.Forward_sort(ranks)
        self.Forward_scheduling()
        print("Forward scheduling")
        print("Start time", self.start_time)
        print("End time", self.end_time)
        print("Execution time", self.end_time-self.start_time)

#----------------------------------------------------------------------------#

    def Compute_forward_ranks(self):
        """Compute the forward rank of each vertex (task)."""

        ranks = [10*len(self.schedule)]*len(self.schedule)
        pos = 0
        for task in self.schedule:
            for waiting_task_id in task[0].waiting_tasks:
                wait_sch_task = self.schedule_id_map[waiting_task_id]
                if self.schedule[wait_sch_task][0].subdomain_id != task[0].subdomain_id:
                    candidate_rank = self.schedule[wait_sch_task][1]
                else:
                    candidate_rank = ranks[wait_sch_task]
                if ranks[pos] > candidate_rank:
                    ranks[pos] = candidate_rank
            pos += 1

        return ranks

#----------------------------------------------------------------------------#

    def Forward_sort(self, ranks):
        """Sort the tasks in the backward iteration."""

        for new_pos in range(len(ranks)):
            min_rank = min(ranks[new_pos:])
            old_pos = new_pos+ranks[new_pos:].index(min_rank)
# Put the task earlier in the list
            ranks[new_pos], ranks[old_pos] = ranks[old_pos], ranks[new_pos]
            self.schedule[new_pos], self.schedule[old_pos] =\
                self.schedule[old_pos], self.schedule[new_pos]

# Sort tasks with the same rank
# This is very slow because there are many tasks with the same ranks
        pos = 0
        size = len(ranks)
        while pos < size:
            end_pos = pos + 1
            rank = ranks[pos]
            while rank == ranks[end_pos]:
                end_pos += 1
                if end_pos == size:
                    break
            for new_pos in range(pos, end_pos):
                earliest_time = 1e100
                for i in range(new_pos, end_pos):
                    if self.schedule[i][1] < earliest_time:
                        earliest_time = self.schedule[i][1]
                        old_pos = i
                self.schedule[new_pos], self.schedule[old_pos] =\
                    self.schedule[old_pos], self.schedule[new_pos]
            pos = end_pos

# Update schedule_id_map
        self.Update_schedule_id_map()

#----------------------------------------------------------------------------#

    def Forward_scheduling(self):
        """Compute a scheduling that starts the tasks as late as possible."""

        n_procs = max(max(self.n_p))+1
# Simulate multiple processors
        part_schedule = []
        for i in range(n_procs):
            part_schedule.append([])
        for task in self.schedule:
            part_schedule[task[0].subdomain_id].append(task)
# Set is a mutable so cannot use [set()]*n_procs
        proc_start_time = [set() for i in range(n_procs)]
        for i in range(n_procs):
            for j in range(self.start_time, self.end_time+100):
                proc_start_time[i].add(j)
# Flag to know if the time of the required task has been updated
        updated_time = [False]*len(self.schedule)

        n_tasks = len(self.schedule)
        n_tasks_done = 0
        old_n_tasks_done = -1
        while n_tasks_done != n_tasks:
            assert old_n_tasks_done != n_tasks_done,\
                'No improvement in Forward scheduling.'
            old_n_tasks_done = n_tasks_done
            for i in range(n_procs):
                if len(part_schedule[i]) != 0:
                    task = part_schedule[i][0]
                    ready = True
                    t = self.start_time
                    if task[0].required_tasks != []:
                        for required_task_id in task[0].required_tasks:
                            local_schedule_id = self.schedule_id_map[required_task_id]
                            required_task = self.schedule[local_schedule_id]
                            ready = updated_time[local_schedule_id]
                            if ready is False:
                                break
                            if t < required_task[2]:
                                t = required_task[2]
                    if ready is True:
                        candidate_t = max(t, min(proc_start_time[i]))
                        enough_free_t = False
                        while enough_free_t is False:
                            enough_free_t = True
                            for j in range(task[0].delta_t):
                                if (candidate_t+j) not in proc_start_time[i]:
                                    enough_free_t = False
                                    candidate_t += 1
                                    break
                        task[1] = candidate_t
                        task[2] = task[1] + task[0].delta_t
                        for j in range(task[0].delta_t):
                            proc_start_time[i].remove(task[1]+j)
                        self.schedule[self.schedule_id_map[task[0].task_id]] = task
                        updated_time[self.schedule_id_map[task[0].task_id]] = True
                        part_schedule[i].pop(0)
                        n_tasks_done += 1

# Update ending time
        self.end_time = self.start_time
        for task in self.schedule:
            if self.end_time < task[2]:
                self.end_time = task[2]
