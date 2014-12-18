# Python code
# Author: Bruno Turcksin
# Date: 2014-12-08 14:31:36.414075

#----------------------------------------------------------------------------#
## Class LOCAL_SEARCH                                                       ##
#----------------------------------------------------------------------------#

"""This class does a local search using a pure downhill method: a task is
randomly chosen and moved to be executed as early as possible."""

import copy
import random
import TASK


class LOCAL_SEARCH(object):
    """This class does a local search using a downhill method."""

    def __init__(self, serial_schedule, n_procs, start_time, end_time, n_samples):
        """The input are: the current best schedule in a 1D list, the number of
        processors, the start time of the schedule, the end time of the schedukle,
        the number of samples to try."""

        super().__init__()
        self.n_procs = n_procs
        self.start_time = start_time
        self.end_time = end_time
        self.n_samples = n_samples
        self.n_tasks = len(serial_schedule)
        self.schedule = []
        for i in range(n_procs):
            self.schedule.append([])
        for task in serial_schedule:
            self.schedule[task[0].subdomain_id].append(task)
        for proc in range(n_procs):
            self.schedule[proc].sort(key=lambda tup: tup[1])

#----------------------------------------------------------------------------#

    def Compute_schedule_id_map(self, schedule):
        """Compute a schedule id map for the given schedule. The key of the map
        is the task ID and the value is the processor and the position in the
        schedule."""

        schedule_id_map = dict()
        for proc in range(self.n_procs):
            pos = 0
            for task in schedule[proc]:
                schedule_id_map[task[0].task_id] = [proc, pos]
                pos += 1

        return schedule_id_map

#----------------------------------------------------------------------------#

    def Compute_sweep(self):
        """Compute a new schedule by moving a random task as early as possible.
        This function returns true if the new schedule is faster than the
        previous one."""

        new_schedule = copy.deepcopy(self.schedule)

        randomize = True
        while randomize is True:
            rand_proc = random.randrange(self.n_procs)
            local_id = random.randrange(len(new_schedule[rand_proc]))
            task_to_move = new_schedule[rand_proc].pop(local_id)
            insert_done = False

            schedule_id_map = self.Compute_schedule_id_map(new_schedule)

            t = self.start_time
            for required_task_id in task_to_move[0].required_tasks:
                local_schedule_id = schedule_id_map[required_task_id]
                required_task = new_schedule[local_schedule_id[0]][local_schedule_id[1]]
                if required_task[2] > t:
                    t = required_task[2]

            for i in range(len(new_schedule[rand_proc])):
                if new_schedule[rand_proc][i][1] > t:
                    insert_done = True
                    new_schedule[rand_proc].insert(i, task_to_move)
                    if i != local_id:
                        print('task_to_move', rand_proc, local_id)
                        print ('new pos', i, 'delta t', i-local_id)
                        randomize = False
                    break
            if insert_done is False:
                new_schedule[rand_proc].insert(local_id, task_to_move)

        schedule_id_map = self.Compute_schedule_id_map(new_schedule)

        candidate_schedule = [[] for i in range(self.n_procs)]
# Set is a mutable so cannot use [set()]*n_procs
        proc_start_time = [set() for i in range(self.n_procs)]
        n_tasks = 0
        for i in range(self.n_procs):
            n_tasks += len(new_schedule[i])
            for j in range(self.start_time, 10*self.end_time):
                proc_start_time[i].add(j)
        assert self.n_tasks == n_tasks, str(self.n_tasks)+' '+str(n_tasks)
        tasks_done = set()

        next_tasks = [0]*self.n_procs
        n_tasks_done_old = -1
        n_tasks_done = 0
        while n_tasks_done != self.n_tasks:
            assert n_tasks_done_old != n_tasks_done, 'Something went wrong. ' +\
                    str(n_tasks_done) + '/' + str(self.n_tasks) + ' ' +\
                    str(next_tasks[rand_proc]) + '/' + str(len(new_schedule[rand_proc]))
            n_tasks_done_old = n_tasks_done
            for proc in range(self.n_procs):
                if len(new_schedule[proc]) != next_tasks[proc]:
                    task = new_schedule[proc][next_tasks[proc]]
                    ready = True
                    t = self.start_time
                    if task[0].required_tasks != []:
                        for required_task_id in task[0].required_tasks:
                            local_schedule_id = schedule_id_map[required_task_id]
                            required_task = new_schedule[local_schedule_id[0]][local_schedule_id[1]]
                            if required_task[0].task_id in tasks_done:
                                ready = True
                                if t < required_task[2]:
                                    t = required_task[2]
                            else:
                                ready = False
                                break
                    if ready is True:
                        candidate_t = max(t, min(proc_start_time[proc]))
                        enough_free_t = False
                        while enough_free_t is False:
                            enough_free_t = True
                            for j in range(task[0].delta_t):
                                if (candidate_t+j) not in proc_start_time[proc]:
                                    enough_free_t = False
                                    candidate_t += 1
                                    break
                        task[1] = candidate_t
                        task[2] = task[1] + task[0].delta_t
                        for j in range(task[0].delta_t):
                            proc_start_time[proc].remove(task[1]+j)
                        candidate_schedule[proc].append(task)
                        local_schedule_id = schedule_id_map[task[0].task_id]
                        tasks_done.add(task[0].task_id)
                        n_tasks_done += 1
                        next_tasks[proc] += 1


# Compute new execution time
        end_time = self.start_time
        for proc in range(self.n_procs):
            for task in candidate_schedule[proc]:
                if end_time < task[2]:
                    end_time = task[2]

        print ('execution time', end_time-self.start_time)
        if end_time < self.end_time:
            self.end_time = end_time
            self.schedule = candidate_schedule
            for proc in range(self.n_procs):
                self.schedule[proc].sort(key=lambda tup: tup[1])
            return True
        else:
            return False

#----------------------------------------------------------------------------#

    def Execute_search(self):
        """Driver for the local search. Try a given number of new schedule."""

        print('Execution time before local search: ',
              self.end_time-self.start_time)
        n_improvements = 0
        for i in range(self.n_samples):
            if self.Compute_sweep() is True:
                n_improvements += 1
            print('Executed local search', i, 'execution time',
                  self.end_time-self.start_time)
        print('Execution time after local search:',
              self.end_time-self.start_time)

        print('Number of steps that improve the scheduling:', n_improvements,
              'out of', self.n_samples)
