# Python code
# Author: Bruno Turcksin
# Date: 2014-12-09 09:27:41.683262

#----------------------------------------------------------------------------#
## MODULE TASK                                                              ##
#----------------------------------------------------------------------------#

"""This module encapsulated all the data necessary to define a task from
Ouranos."""


import collections


TASK = collections.namedtuple('TASK', ['subdomain_id', 'task_id',
                                       'required_tasks', 'waiting_tasks',
                                       'delta_t', 'pos', 'dir', 'b_level'])
TASK.__doc__ = """Represent a task in Ouranos.

Each task is associated to a processor (subdomain_id), has an id (task_id), has
set of tasks that must be executed before it can be executed, has set of tasks
that are waiting for this task to be executed, and needs a certain time to be
executed (delta_t)."""
